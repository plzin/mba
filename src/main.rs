use rug::Integer;
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

/// Rewrite a linear combination of uniform expression using a set of uniform expressions modulo n.
fn rewrite(expr: &LUExpr, ops: &[UExpr], n: &Integer) -> Option<LUExpr> {
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

    // Solve the system.
    let l = diophantine::solve_modular(&a, &b, n);

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = l.sample_point() % n;

    // Put it in a LUExpr.
    let v = solution.iter()
        .zip(ops.iter())
        .map(|(m, e)| (m.clone(), e.clone()))
        .collect();

    Some(LUExpr(v))
}

/// Provides an algorithm for obfuscating "any" kind of expression
/// using the two primitives of rewriting and permutation polynomials.
/// The members of this struct are used
/// to cache things for the recursive algorithm.
struct ExampleObfuscation {
    /// The integer n, everything is mod 2^n.
    n: u32,

    /// 2^n.
    pow: Integer,

    /// The operations used for rewriting things.
    ops: Vec<UExpr>,

    /// The variables occuring in the expression.
    vars: Vec<char>,

    /// The quotient ring of permutation polynomials.
    qr: perm_poly::QuotientRing,

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
            let (p, q) = perm_poly::perm_pair(&self.qr,
                Uniform::from(1..3).sample(&mut self.rng));

            let p = p.to_expr().substitute('x', e);
            *e = Rc::new(q.to_expr().substitute('x', &mut Rc::new(p)));
        }
    }

    fn init(e: &Expr, n: u32) -> Self {
        // Get the variables in the expression.
        let mut vars = Vec::new();
        e.vars(&mut vars);
        vars.sort();
        vars.dedup();

        let ops = vec![
            UExpr::from_string("X&Y").unwrap(),
            UExpr::from_string("X^Y").unwrap(),
            UExpr::from_string("~(X^Y^Z)").unwrap(),
            UExpr::from_string("~X").unwrap(),
            UExpr::from_string("~Y").unwrap(),
            UExpr::from_string("Y").unwrap(),
            UExpr::from_string("Z").unwrap(),
            UExpr::from_string("Y&Z").unwrap(),
            UExpr::from_string("X|Z").unwrap(),
            UExpr::from_string("~X&Z").unwrap(),
            UExpr::from_string("Y|~Z").unwrap(),
        ];

        let pow = Integer::from(1) << n;
        let qr = perm_poly::QuotientRing::init(n);
        let rng = rand::thread_rng();

        Self { n, pow, ops, vars, qr, rng }
    }
}

// fn main() {
//     let r = perm_poly::QuotientRing::init(8);
//     let (p, q) = perm_poly::perm_pair(&r, 3);
//     println!("p(x) = {}", p);
//     println!("q(x) = {}", q);
//     //perm_poly::check_inverse_64();
// }

fn main() {
    let e = Expr::from_string("2*x+y+2").unwrap();
    let e = ExampleObfuscation::obfuscate(e, 32);
    //println!("{:?}", e);
    e.print_simple();
    // let sub = std::rc::Rc::new(Expr::from_string("x+y").unwrap());
    // let mut e = Expr::Mul(sub.clone(), sub);
    // println!("{:?}", e);
    // e.print_simple();

    // e.substitute('x', Expr::from_string("x+~z|y").unwrap());
    // e.print_simple();

    //let e = Expr::from_string("~x&y+z");
    //println!("{:?}", e);
    // let e = LUExpr::from_string("~x&y|z+123");
    // println!("{:?}", e);

    // let n = 32;

    // let expr = LUExpr::from_string("2*x+y+2").unwrap();
    // // let expr = LUExpr::from_string("128").unwrap();

    // let ops = [
    //     UExpr::from_string("x&y").unwrap(),
    //     UExpr::from_string("x^y").unwrap(),
    //     UExpr::from_string("~(x^y^z)").unwrap(),
    //     UExpr::from_string("~x").unwrap(),
    //     UExpr::from_string("~y").unwrap(),
    //     UExpr::from_string("y").unwrap(),
    //     UExpr::from_string("z").unwrap(),
    //     UExpr::from_string("y&z").unwrap(),
    //     UExpr::from_string("x|z").unwrap(),
    //     UExpr::from_string("~x&z").unwrap(),
    //     UExpr::from_string("y|(~z)").unwrap(),
    // ];

    // let e = rewrite(&expr, &ops, &Integer::new().set_bit(n, true));
    // if let Some(e) = e {
    //     println!("{}", e);
    // } else {
    //     println!(
    //         "Rewriting the expression with the \
    //         given operations isn't possible.");
    // }


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
