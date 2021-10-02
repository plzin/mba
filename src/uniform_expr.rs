#![allow(dead_code)]

use std::ops::{Index, IndexMut};

use rug::Integer;

use crate::int_from_it;

/// LUExpr is short for "Linear combination of Uniform Expressions"
/// These are the expressions for which rewrite rules can be efficiently
/// generated.
#[derive(Clone, Debug)]
pub struct LUExpr(pub Vec<(Integer, UniformExpr)>);

impl LUExpr {
    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self, v: &mut Vec<char>) {
        for (_, e) in &self.0 {
            e.vars(v);
        }
    }

    /// Evaluate an expression with a valuation for the occuring variables.
    pub fn eval(&self, v: &Valuation) -> Integer {
        self.0.iter()
            .map(|(i, e)| i * e.eval(v))
            .sum()
    }

    /// Parse a string to an expression.
    /// Note that this function is extremly limited
    /// and expects very specific syntax.
    /// It is used for convenience when testing things and
    /// not really meant to be used by something outside this crate.
    pub(crate) fn from_string<T: ToString>(s: T) -> Option<Self> {
        let mut s = s.to_string();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        // This stores the current linear combination.
        let mut v = Vec::new();

        // Loop over the string/the summands.
        loop {
            // Are there still characters left?
            // If not then we're done.
            let c = match it.peek() {
                None => return Some(Self(v)),
                Some(c) => *c,
            };

            // If this is a digit then we expect num*UExpr.
            if c.is_ascii_digit() {
                // Parse the number.
                let num = int_from_it(&mut it)?;

                // Is it the expected '*'?
                let lookahead = it.peek();
                match lookahead {
                    Some('*') => {
                        it.next();

                        // Parse the UniformExpr.
                        let e = UniformExpr::parse(&mut it, 0)?;

                        // Push it.
                        v.push((num, e));
                    },

                    // If this is a different character then we push num*1.
                    _ => v.push((num, UniformExpr::One)),
                }
            } else {
                // We don't have a factor so just parse the UniformExpr.
                let e = UniformExpr::parse(&mut it, 0)?;

                // Push 1*e.
                v.push((1.into(), e));
            }

            // If the next character is not a plus then we are done.
            match it.peek() {
                Some('+') => it.next(), // Skip the +.
                _ => return Some(Self(v)),
            };
        }
    }
}

impl std::fmt::Display for LUExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.0.iter().filter(|(i, _)| *i != 0);
        let (i, e) = match iter.next() {
            Some(t) => t,
            None => return write!(f, "0"),
        };

        if *i == 1 {
            write!(f, "{}", e)?;
        } else {
            write!(f, "{}*({})", i, e)?;
        }

        for (i, e) in iter {
            write!(f, " + ")?;
            if *i == 1 {
                write!(f, "{}", e)?;
            } else {
                write!(f, "{}*({})", i, e)?;
            }
        }

        Ok(())
    }
}

/// Represents an expression that is uniform on all bits.
#[derive(Clone, Debug)]
pub enum UniformExpr {
    One,
    Var(char),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Xor(Box<Self>, Box<Self>),
    Not(Box<Self>),
}

impl UniformExpr {
    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self, v: &mut Vec<char>) {
        use UniformExpr::*;
        match self {
            One      => {},
            Var(c)          => v.push(*c),
            And(e1, e2)     => { e1.vars(v); e2.vars(v) },
            Or(e1, e2)      => { e1.vars(v); e2.vars(v) },
            Xor(e1, e2)     => { e1.vars(v); e2.vars(v) },
            Not(e)          => e.vars(v),
        }
    }

    /// Evaluate an expression with a valuation for the occuring variables.
    pub fn eval(&self, v: &Valuation) -> Integer {
        use UniformExpr::*;
        match self {
            One             => 1.into(),
            Var(c)          => v[*c].clone(),
            And(e1, e2)     => (e1.eval(v) & e2.eval(v)) & 1,
            Or(e1, e2)      => (e1.eval(v) | e2.eval(v)) & 1,
            Xor(e1, e2)     => (e1.eval(v) ^ e2.eval(v)) & 1,
            Not(e)          => (!e.eval(v)) & 1,
        }
    }

    fn write_safe(
        e1: &Self, e2: &Self, op: char, f: &mut std::fmt::Formatter<'_>
    ) -> std::fmt::Result {
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

    fn precedence(op: char) -> usize {
        match op {
            '|' => 1,
            '^' => 2,
            '&' => 3,
            _ => 0,
        }
    }

    /// Parse a string to an expression.
    pub(crate) fn from_string<T: ToString>(s: T) -> Option<Self> {
        let mut s = s.to_string();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    pub(self) fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize
    ) -> Option<Self> {
        use UniformExpr::*;

        let c = *it.peek()?;

        let mut e = if c == '(' {
            it.next();
            let e = Self::parse(it, 0)?;
            match it.next() {
                Some(')') => e,
                _ => return None,
            }
        } else if c == '~' {
            it.next();
            let e = Self::parse(it, 15)?;
            Not(Box::new(e))
        } else if c.is_alphabetic() {
            it.next();
            Var(c)
        } else if c == '1' {
            One
        } else {
            return None;
        };

        loop {
            let c = match it.peek() {
                None => return Some(e),
                Some(c) => *c,
            };

            let op_pre = Self::precedence(c);
            if op_pre <= pre {
                return Some(e);
            }

            // If the current operators precedence is higher than
            // the one whos subexpression we are currently parsing
            // then we need to finish this operator first.
            it.next();
            let rhs = Box::new(Self::parse(it, op_pre)?);
            let lhs = Box::new(e);
            e = match c {
                '&' => And(lhs, rhs),
                '|' => Or(lhs, rhs),
                '^' => Xor(lhs, rhs),
                _ => return None,
            };
        }
    }
}

impl std::fmt::Display for UniformExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use UniformExpr::*;
        match self {
            One => write!(f, "1"),
            Var(c) => write!(f, "{}", c),
            And(e1, e2)   => Self::write_safe(e1, e2, '&', f),
            Or(e1, e2)    => Self::write_safe(e1, e2, '|', f),
            Xor(e1, e2)   => Self::write_safe(e1, e2, '^', f),
            Not(e) =>
                if let Self::Var(c) = e.as_ref() {
                    write!(f, "~{}", c)
                } else {
                    write!(f, "~({})", e)
                },
        }
    }
}

/// Stores values that should be substituted into variables.
#[derive(Debug)]
pub struct Valuation {
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
