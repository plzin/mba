#![allow(unused)]

use rug::{Integer, Complete};
use rand::{
    Rng, distributions::Uniform,
    prelude::Distribution, rngs::ThreadRng
};

use std::rc::Rc;

mod matrix;
use matrix::Matrix;

mod vector;
use vector::Vector;

mod expr;
use expr::*;

mod uniform_expr;
use uniform_expr::*;

mod diophantine;
mod poly;
mod perm_poly;

/// Rewrite a linear combination of uniform expression using a set of uniform expressions modulo `n`.
/// `n` has to be a power of two.
fn rewrite(expr: &LUExpr, ops: &[LUExpr], n: &Integer) -> Option<LUExpr> {
    // Find all variables we have access to.
    // This includes variables in the expression as well as potentially the ops.

    let mut v = Vec::new();
    expr.vars(&mut v);

    ops.iter().for_each(|e| e.vars(&mut v));

    // Remove duplicates and sort.
    v.sort();
    v.dedup();

    assert!(v.len() <= 63, "More than 63 variables are currently not supported \
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
            val[*c] = -Integer::from((i >> j) & 1);
        }

        // println!("{:?}", val);

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = e.eval(&val);
        }

        // Write the desired result into the vector.
        b[i] = expr.eval(&val);
    }

    // Solve the system.
    let l = diophantine::solve_modular(&a, &b, n);

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = l.sample_point() % n;

    // Put it in a LUExpr.
    let mut v = Vec::new();
    for (c, o) in solution.iter().zip(ops.iter()) {
        for (d, e) in &o.0 {
            // Is the UExpr already in the linear combination?
            match v.iter_mut().find(|(_, f)| f == e) {
                Some((f, _)) => *f += c * d,
                None => v.push(((c * d).complete(), e.clone())),
            }
        }
    }

    Some(LUExpr(v))
}

// This generates an 8-bit permutation polynomial of degree 3 and its inverse.
//fn main() {
//    let r = perm_poly::QuotientRing::init(8);
//    let (p, q) = perm_poly::perm_pair(&r, 3);
//    println!("p(x) = {}", p);
//    println!("q(x) = {}", q);
//}

// Rewrite an expression using linear MBA.
//fn main() {
//    let expr = LUExpr::from_string("x+y").unwrap();
//    let ops = &[
//        LUExpr::from_string("x&y").unwrap(),
//        LUExpr::from_string("x^y").unwrap(),
//    ];
//
//    match rewrite(&expr, ops, &256.into()) {
//        None => println!("Can't rewrite."),
//        Some(e) => println!("{}", e),
//    }
//}

// Rewrite a constant using non-linear MBA.
fn main() {
    // Initialize stuff.
    let bits = 8;
    let n = Integer::from(1) << bits;
    let zi = perm_poly::ZeroIdeal::init(bits);

    // This is the constant that will be obfuscated.
    let input = Integer::from(42);

    // Obfuscate a different constant and apply a
    // polynomial at the end to get the real one.
    //let p = poly::Poly::from_vec(vec![0, 1, 2]);
    //let q = perm_poly::compute_inverse(&p, &zi);
    //let input = q.eval_bits(&input, bits);

    // The expression to obfuscate.
    let expr = LUExpr::constant(input);

    // The operations for rewriting.
    let ops = &[
        LUExpr::from_string("x^y").unwrap(),
        LUExpr::from_string("x&y").unwrap(),
        LUExpr::from_string("~(x|y)").unwrap(),
        LUExpr::from_string("~x").unwrap(),
        LUExpr::from_string("~y").unwrap(),
        LUExpr::from_string("z").unwrap(),
        LUExpr::from_string("y&z").unwrap(),
        LUExpr::from_string("x|z").unwrap(),
        LUExpr::from_string("~x&z").unwrap(),
        LUExpr::from_string("y|~z").unwrap(),
    ];

    // Rewrite the constant.
    let obf = rewrite(&expr, ops, &n).unwrap();

    // Contains the linear combination of the multiplication.
    let mut v: Vec<(Integer, UExpr, UExpr)> = Vec::new();

    // Obfuscate each coefficient.
    for (c, e) in &obf.0 {

        // Rewrite the coefficient.
        let coeff = LUExpr::constant(c.clone());
        let mut ops = ops.to_vec();
        ops.retain(|f| &f.0[0].1 != e);
        let coeff = rewrite(&coeff, &ops, &n).unwrap();

        // Add all coefficients of the obfuscated coefficient
        // to the linear combination.
        for (d, f) in coeff.0 {
            let e = e.clone();
            let (fst, snd) = if e <= f { (e, f) } else { (f, e) };

            match v.iter_mut().find(|(_, a, b)| &fst == a && &snd == b) {
                Some((c, _, _)) => *c += d,
                None => v.push((d, fst, snd)),
            }
        }
    }

    for (c, _, _) in &mut v {
        c.keep_bits_mut(bits);
    }

    // Convert a summand to an expression.
    let term_to_expr = |c: &Integer, e: &UExpr, f: &UExpr| {
        Expr::Mul(
            Rc::new(Expr::Const(c.clone())),
            Rc::new(Expr::Mul(Rc::new(e.to_expr()), Rc::new(f.to_expr())))
        )
    };

    // Convert the linear combination to an `Expr`.
    let mut iter = v.iter().filter(|(c, _, _)| c != &Integer::ZERO);
    let mut expr = match iter.next() {
        Some((c, e, f)) => term_to_expr(c, e, f),
        None => Expr::Const(Integer::new()),
    };

    for (c, e, f) in iter {
        let a = term_to_expr(c, e, f);
        expr = Expr::Add(Rc::new(expr), Rc::new(a));
    }

    //let mut re = Rc::new(expr);
    //expr = p.to_expr().substitute('x', &mut re);

    expr.print_simple();
}

// Solve a system of linear congruences.
//fn main() {
//    let solution = diophantine::solve_modular(
//        &Matrix::from_array([[3, 5], [4, 2]]),
//        &Vector::from_slice(&[2.into(), 0.into()]),
//        &256.into()
//    );
//
//    if (solution.is_empty()) {
//        println!("No solution");
//    } else {
//        println!("Off: {:?}\nBasis:", solution.offset);
//        for b in &solution.basis {
//            println!("{:?}", b);
//        }
//    }
//}

//fn main() {
//    // The example will obfuscate x+y where x and y are 32-bit integers.
//    // It will include a 32-bit integer z whose value doesn't matter.
//    // The variables in the string should be lower case letters.
//    // Note that the example obfuscation is quite limited and only meant as
//    // a starting point and an example of how to combine the primitives.
//    let e = Expr::from_string("x+y").unwrap();
//    let e = ExampleObfuscation::obfuscate(e, 32);
//    e.print_simple();
//}

/// Provides an algorithm for obfuscating "any" kind of expression
/// using the two primitives of rewriting and permutation polynomials.
/// The members of this struct are used
/// to cache things for the recursive algorithm.
struct ExampleObfuscation {
    /// 2^n.
    pow: Integer,

    /// The operations used for rewriting things.
    ops: Vec<LUExpr>,

    /// The variables occurring in the expression.
    vars: Vec<char>,

    /// The quotient ring of permutation polynomials.
    qr: perm_poly::ZeroIdeal,

    /// Random number generator.
    rng: ThreadRng,
}

impl ExampleObfuscation {
    pub fn obfuscate(e: Expr, n: u32) -> Expr {
        let mut o = ExampleObfuscation::init(&e, n);

        // Yes we need to create an unnecessary Rc once.
        // But the code is much simpler this way.
        let mut e = Rc::new(e);

        o.obfuscate_impl(&mut e);

        let mut x = Rc::new(Expr::Var(*o.vars.get(0).unwrap_or(&'x')));
        let mut y = Rc::new(Expr::Var(*o.vars.get(1).unwrap_or(&'y')));
        let mut z = Rc::new(Expr::Var(*o.vars.get(2).unwrap_or(&'z')));

        e.as_ref().clone().substitute('X', &mut x)
            .substitute('Y', &mut y)
            .substitute('Z', &mut z)
    }

    fn obfuscate_impl(&mut self, e: &mut Rc<Expr>) {
        if self.rng.gen::<f32>() < 0.9 {
            match Rc::make_mut(e) {
                Expr::Const(i) => {
                    let u = LUExpr::constant(i.clone());
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr().into();
                    }
                },
                Expr::Var(v) => {
                    let v = *v;
                    let u = LUExpr::var('X');
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', &mut Rc::new(Expr::Var(v)))
                            .into();
                    }
                },
                Expr::Add(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X+Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', l)
                            .substitute('Y', r)
                            .into();
                    }
                },
                Expr::Sub(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X-Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', l)
                            .substitute('Y', r)
                            .into();
                    }
                },
                Expr::Mul(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let f = Expr::from_string("(X&Y)*(X|Y)+(X&~Y)*(~X&Y)")
                        .unwrap();

                    *e = f.substitute('X', l).substitute('Y', r).into();
                },
                Expr::Div(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                },
                Expr::Neg(i) => {
                    self.obfuscate_impl(i);
                    let u = LUExpr::from_string("-X").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', i)
                            .into();
                    }
                },
                Expr::And(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X&Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', l)
                            .substitute('Y', r)
                            .into();
                    }
                },
                Expr::Or(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X|Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', l)
                            .substitute('Y', r)
                            .into();
                    }
                },
                Expr::Xor(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X^Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', l)
                            .substitute('Y', r)
                            .into();
                    }
                },
                Expr::Not(i) => {
                    self.obfuscate_impl(i);
                    let u = LUExpr::from_string("~X").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, &self.pow) {
                        *e = f.to_expr()
                            .substitute('X', i)
                            .into();
                    }
                },
            };
        }

        if self.rng.gen::<f32>() < 0.3 {
            let degree = Uniform::from(2..4).sample(&mut self.rng);
            let (p, q) = perm_poly::perm_pair(&self.qr, degree);
            let p = p.to_expr().substitute('x', e);
            *e = Rc::new(q.to_expr().substitute('x', &mut Rc::new(p)));
        }
    }

    fn init(e: &Expr, n: u32) -> Self {
        // Get the variables in the expression.
        let mut vars = Vec::new();
        e.vars(&mut vars);
        vars.sort_unstable();
        vars.dedup();

        let ops = vec![
            LUExpr::from_string("X&Y").unwrap(),
            LUExpr::from_string("X^Y").unwrap(),
            LUExpr::from_string("~(X^Y^Z)").unwrap(),
            LUExpr::from_string("~X").unwrap(),
            LUExpr::from_string("~Y").unwrap(),
            LUExpr::from_string("Y").unwrap(),
            LUExpr::from_string("Z").unwrap(),
            LUExpr::from_string("Y&Z").unwrap(),
            LUExpr::from_string("X|Z").unwrap(),
            LUExpr::from_string("~X&Z").unwrap(),
            LUExpr::from_string("Y|~Z").unwrap(),
        ];

        let pow = Integer::from(1) << n;
        let qr = perm_poly::ZeroIdeal::init(n);
        let rng = rand::thread_rng();

        Self { pow, ops, vars, qr, rng }
    }
}
