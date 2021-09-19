#![allow(dead_code)]

use std::ops::{Index, IndexMut};
use std::boxed::Box;

use rug::Integer;

mod matrix;
use matrix::Matrix;

mod vector;
use vector::Vector;

mod diophantine;

#[derive(Debug)]
struct Valuation {
    /// The key value pairs are stored as a Vector
    /// because I doubt a hashmap/tree would be faster
    /// when there are so few variables.
    vals: Vec<(char, Integer)>,
}

impl Valuation {
    /// Initializes a valuation from a list of variables
    /// each of which will be Initialized to 0.
    pub fn zero(vars: &Vec<char>) -> Self {
        let vals = vars.iter()
            .map(|c| (*c, 0.into()))
            .collect();

        Self { vals }
    }
}

impl Index<char> for Valuation {
    type Output = Integer;
    fn index(&self, index: char) -> &Self::Output {
        &self.vals.iter()
            .find(|(name, _)| *name == index)
            .unwrap().1
    }
}

impl IndexMut<char> for Valuation {
    fn index_mut(&mut self, index: char) -> &mut Self::Output {
        &mut self.vals.iter_mut()
            .find(|(name, _)| *name == index)
            .unwrap().1
    }
}

#[derive(Clone, Debug)]
enum Expr {
    Var(char),
    Add(Vec<Expr>),
    Mul(Integer, Box<Expr>),
    And(Box<Expr>, Box<Expr>),
    Or(Box<Expr>, Box<Expr>),
    Xor(Box<Expr>, Box<Expr>),
    Not(Box<Expr>),
}

impl Expr {
    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self, v: &mut Vec<char>) {
        match self {
            Self::Var(c)        => v.push(*c),
            Self::Add(es)       => es.iter().for_each(|e| e.vars(v)),
            Self::Mul(_, e)     => e.vars(v),
            Self::And(e1, e2)   => { e1.vars(v); e2.vars(v) },
            Self::Or(e1, e2)    => { e1.vars(v); e2.vars(v) },
            Self::Xor(e1, e2)   => { e1.vars(v); e2.vars(v) },
            Self::Not(e)        => e.vars(v),
        }
    }

    /// Evaluate an expression with a valuation for the occuring variables.
    pub fn eval(&self, v: &Valuation) -> Integer {
        match self {
            Self::Var(c)        => v[*c].clone(),
            Self::Add(es)       => es.iter().map(|e| e.eval(v)).sum(),
            Self::Mul(i, e)     => i * e.eval(v),
            Self::And(e1, e2)   => (e1.eval(v) & e2.eval(v)) & 1,
            Self::Or(e1, e2)    => (e1.eval(v) | e2.eval(v)) & 1,
            Self::Xor(e1, e2)   => (e1.eval(v) ^ e2.eval(v)) & 1,
            Self::Not(e)        => (!e.eval(v)) & 1,
        }
    }

    fn write_safe(e1: &Expr, e2: &Expr, op: char, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Self::Var(c) = e1 {
            write!(f, "{} {}", c, op)?;
        } else {
            write!(f, "({}) {}", e1, op)?;
        }

        if let Self::Var(c) = e2 {
            write!(f, " {}", c)
        } else {
            write!(f, " ({})", e2)
        }
    }

    /// Parse a string to an expression.
    /// Note that this function is extremly limited
    /// and expects very specific syntax.
    pub fn from_string<T: ToString>(s: T) -> Option<Self> {
        let mut s = s.to_string();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it)
    }

    fn parse(s: &mut std::iter::Peekable<std::str::Chars>) -> Option<Self> {
        let lookahead = match s.next() {
            None => return None,
            Some(c) => c,
        };

        let lhs = if lookahead == '(' {
            let e = Self::parse(s)?;
            s.next();
            e
        } else if lookahead == '~' {
            Self::Not(Box::new(Self::parse(s)?))
        } else if lookahead.is_ascii_digit() {
            let mut num: Integer = lookahead.to_digit(10).unwrap().into();

            let mut lookahead = s.peek();
            while let Some(c) = lookahead {
                if !c.is_ascii_digit() {
                    break;
                }

                num *= 10;
                num += c.to_digit(10).unwrap();
                s.next();
                lookahead = s.peek();
            }

            match lookahead {
                Some('*') => {
                    s.next();
                    let e = Self::parse(s)?;
                    return Some(Self::Mul(num, Box::new(e)));
                },
                _ => return None,
            };
        } else if lookahead.is_alphabetic() {
            let v = lookahead;
            Self::Var(v)
        } else {
            return None;
        };

        let c = match s.peek() {
            None | Some(')') => return Some(lhs),
            Some(c) => *c,
        };

        s.next();

        let rhs = Self::parse(s)?;

        Some(match c {
            '+' => Expr::Add(vec![lhs, rhs]),
            '&' => Expr::And(Box::new(lhs), Box::new(rhs)),
            '|' => Expr::Or(Box::new(lhs), Box::new(rhs)),
            '^' => Expr::Xor(Box::new(lhs), Box::new(rhs)),
            _ => return None,
        })
    }
}

impl std::fmt::Display for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Var(c) => write!(f, "{}", c),
            Self::Add(v) => {
                write!(f, "{}", v[0])?;
                for e in &v[1..] {
                    write!(f, " + {}", e)?;
                }
                Ok(())
            },
            Self::Mul(i, e)     =>
                if let Self::Var(c) = e.as_ref() {
                    write!(f, "{}*{}", i, c)
                } else {
                    write!(f, "{}*({})", i, e)
                },
            Self::And(e1, e2)   => Self::write_safe(e1, e2, '&', f),
            Self::Or(e1, e2)    => Self::write_safe(e1, e2, '|', f),
            Self::Xor(e1, e2)   => Self::write_safe(e1, e2, '^', f),
            Self::Not(e)        =>
                if let Self::Var(c) = e.as_ref() {
                    write!(f, "~{}", c)
                } else {
                    write!(f, "~({})", e)
                },
        }
    }
}

/// Rewrite an expression using a set of operations modulo n.
fn rewrite(expr: &Expr, ops: &[Expr], n: &Integer) -> Option<Expr> {
    // Find all variables we have access to.
    // This includes variables in the expression as well as potentially the ops.

    let mut v = Vec::new();
    expr.vars(&mut v);

    ops.iter().for_each(|e| e.vars(&mut v));

    // Remove duplicates and sort.
    v.sort();
    v.dedup();

    assert!(v.len() <= 63, "More than 63 variables are currently not supported
            (You wouldn't be able to run this anyways).");

    let mut val = Valuation::zero(&v);

    let rows = (1 as usize) << v.len();
    let cols = ops.len();

    let mut a = Matrix::zero(rows, cols);
    let mut b = Vector::zero(rows);

    // Build up the matrix.
    for i in 0..rows {
        let row = a.row_mut(i);

        // Initialize the valuation.
        for (j, c) in v.iter().enumerate() {
            val[*c] = ((i >> j) & 1).into();
        }

        // println!("{:?}", val);

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = e.eval(&val);
        }

        // Write the desired result into the vector.
        b[i] = expr.eval(&val);
    }

    println!("a: {:?}", a);
    println!("b: {:?}", b);

    // Solve the system.
    let l = diophantine::solve_modular(&a, &b, n);

    println!("l: {:?}", l);

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = l.sample_point() % n;
    println!("{:?}", solution);

    let v = solution.iter()
        .zip(ops.iter())
        .map(|(m, e)| Expr::Mul(m.clone(), Box::new(e.clone())))
        .collect();

    Some(Expr::Add(v))
}

fn main() {
    let n = ((1 as usize) << 8).into();

    let expr = Expr::from_string("x+y").unwrap();

    let ops = [
        Expr::from_string("x&y").unwrap(),
        Expr::from_string("x^y").unwrap(),
        Expr::from_string("~(x^y)").unwrap(),
        Expr::from_string("~x").unwrap(),
        Expr::from_string("~y").unwrap(),
        Expr::from_string("y").unwrap(),
    ];

    if let Some(e) = rewrite(&expr, &ops, &n) {
        println!("{}", e);
    } else {
        println!("Rewriting the expression with the given operations isn't possible.");
    }


    // let a = Matrix::from_array([
    //     [0, 0, 1, 1, 1, 0],
    //     [1, 0, 0, 0, 1, 0],
    //     [1, 0, 0, 1, 0, 1],
    //     [0, 1, 1, 0, 0, 1],
    // ]);

    // let b = Vector::from_array(
    //     [0, 1, 1, 2]
    // );

    // println!("A: {:?}", a);

    // println!("b: {:?}", b);

    // let l = diophantine::solve_modular(&a, &b, &n);

    // if l.is_empty() {
    //     println!("System does not have a solution");
    //     return;
    // }

    // println!("x: {:?}", l);

    // for _ in 0..10 {
    //     let mut s = l.offset.clone();

    //     for b in &l.basis {
    //         let f = Integer::from(random::<u8>());
    //         s += b * &f;
    //     }

    //     s %= &n;

    //     println!("s = {:?}", s);

    //     let r = (&a * &s) % &n;
    //     println!("r = {:?}", r);

    //     assert!(r == b, "Invalid solution found");
    // }
}
