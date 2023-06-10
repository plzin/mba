#![allow(unused)]
#![feature(ptr_metadata)]

// You can only compile this for 64-bits because I use some hacky stuff
// with vector views. This can be removed once rust allows supports
// custom dynamically sized types.
#[cfg(not(target_pointer_width = "64"))]
compile_error!("This crate only works on 64-bit systems.");


use rug::{Integer, Complete};
use rand::{
    Rng, distributions::Uniform,
    prelude::Distribution, rngs::ThreadRng
};
#[cfg(feature = "z3")]
use z3::{self, ast::Ast};

use std::rc::Rc;

mod matrix;
use matrix::*;

mod vector;
use vector::*;

mod expr;
use expr::*;

mod uniform_expr;
use uniform_expr::*;

mod diophantine;
mod lattice;
mod poly;
mod perm_poly;

/// Rewrite a linear combination of uniform expression
/// using a set of uniform expressions modulo `2^bits`.
fn rewrite(
    expr: &LUExpr,
    ops: &[LUExpr],
    bits: u32,
    randomize: bool
) -> Option<LUExpr> {
    // Find all variables we have access to.
    // This includes variables in the expression as well as potentially the ops.
    let mut v = std::collections::BTreeSet::new();
    expr.vars_impl(&mut v);
    ops.iter().for_each(|e| e.vars_impl(&mut v));
    let v: Vec<_> = v.into_iter().collect();

    assert!(v.len() <= 63, "More than 63 variables are currently not supported \
            (You wouldn't be able to run this anyways).");

    let mut val = Valuation::zero(v.clone());

    let rows = 1 << v.len();
    let cols = ops.len();

    let mut a = Matrix::zero(rows, cols);
    let mut b = OwnedVector::zero(rows);

    // Build up the matrix.
    for i in 0..rows {
        let row = a.row_mut(i);

        // Initialize the valuation.
        for (j, c) in v.iter().enumerate() {
            val[c] = -Integer::from((i >> j) & 1);
        }

        // println!("{:?}", val);

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = e.eval(&val).keep_bits(bits);
        }

        // Write the desired result into the vector.
        b[i] = expr.eval(&val).keep_bits(bits);
    }

    // Solve the system.
    let l = diophantine::solve_modular(&a, &b, &(Integer::from(1) << bits));

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = if randomize {
        l.sample_point(bits)
    } else {
        l.offset
    };

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

    for (i, _) in &mut v {
        i.keep_signed_bits_mut(bits);
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
//fn main() {
//    // Initialize stuff.
//    let bits = 8;
//    let n = Integer::from(1) << bits;
//    let zi = perm_poly::ZeroIdeal::init(bits);
//
//    // This is the constant that will be obfuscated.
//    let input = Integer::from(42);
//
//    // Obfuscate a different constant and apply a
//    // polynomial at the end to get the real one.
//    //let p = poly::Poly::from_vec(vec![0, 1, 2]);
//    //let q = perm_poly::compute_inverse(&p, &zi);
//    //let input = q.eval_bits(&input, bits);
//
//    // The expression to obfuscate.
//    let expr = LUExpr::constant(input.clone());
//
//    // The operations for rewriting.
//    let ops = &[
//        LUExpr::from_string("x^y").unwrap(),
//        LUExpr::from_string("x&y").unwrap(),
//        LUExpr::from_string("~(x|y)").unwrap(),
//        LUExpr::from_string("~x").unwrap(),
//        LUExpr::from_string("~y").unwrap(),
//        LUExpr::from_string("z").unwrap(),
//        LUExpr::from_string("y&z").unwrap(),
//        LUExpr::from_string("x|z").unwrap(),
//        LUExpr::from_string("~x&z").unwrap(),
//        LUExpr::from_string("y|~z").unwrap(),
//    ];
//
//    // Rewrite the constant.
//    let obf = rewrite(&expr, ops, bits, true).unwrap();
//
//    // Contains the linear combination of the multiplication.
//    let mut v: Vec<(Integer, UExpr, UExpr)> = Vec::new();
//
//    // Obfuscate each coefficient.
//    for (c, e) in &obf.0 {
//
//        // Rewrite the coefficient.
//        let coeff = LUExpr::constant(c.clone());
//        let mut ops = ops.to_vec();
//        ops.retain(|f| &f.0[0].1 != e);
//        let coeff = rewrite(&coeff, &ops, bits, true).unwrap();
//
//        // Add all coefficients of the obfuscated coefficient
//        // to the linear combination.
//        for (d, f) in coeff.0 {
//            let e = e.clone();
//            let (fst, snd) = if e <= f { (e, f) } else { (f, e) };
//
//            match v.iter_mut().find(|(_, a, b)| &fst == a && &snd == b) {
//                Some((c, _, _)) => *c += d,
//                None => v.push((d, fst, snd)),
//            }
//        }
//    }
//
//    for (c, _, _) in &mut v {
//        c.keep_bits_mut(bits);
//    }
//
//    // Convert a summand to an expression.
//    let term_to_expr = |c: &Integer, e: &UExpr, f: &UExpr| {
//        Expr::Mul(
//            Rc::new(Expr::Const(c.clone())),
//            Rc::new(Expr::Mul(Rc::new(e.to_expr()), Rc::new(f.to_expr())))
//        )
//    };
//
//    // Convert the linear combination to an `Expr`.
//    let mut iter = v.iter().filter(|(c, _, _)| c != &Integer::ZERO);
//    let mut expr = match iter.next() {
//        Some((c, e, f)) => term_to_expr(c, e, f),
//        None => Expr::Const(Integer::new()),
//    };
//
//    for (c, e, f) in iter {
//        let a = term_to_expr(c, e, f);
//        expr = Expr::Add(Rc::new(expr), Rc::new(a));
//    }
//
//    //let mut re = Rc::new(expr);
//    //expr = p.to_expr().substitute("x", &mut re);
//
//
//    //if cfg!(feature = "z3") {
//    //    let cfg = z3::Config::new();
//    //    let ctx = z3::Context::new(&cfg);
//    //    let bv = expr.to_z3_bv(&ctx, bits);
//
//    //    let solver = z3::Solver::new(&ctx);
//    //    let assertion = bv._eq(&int_to_bv(&ctx, bits, &input)).not();
//    //    println!("assertion: {}", assertion);
//    //    solver.assert(&assertion);
//    //    println!("Using z3 to check that the rewritten expression is equivalent.");
//    //    let now = std::time::Instant::now();
//    //    let result = solver.check();
//    //    let dur = now.elapsed().as_secs_f32();
//    //    match result {
//    //        z3::SatResult::Sat => println!("Found counterexample:\n{}",
//    //            solver.get_model().unwrap()),
//    //        z3::SatResult::Unsat => println!("Verification successful."),
//    //        z3::SatResult::Unknown => println!("z3 aborted."),
//    //    }
//    //    println!("z3 verification took {dur:.2} seconds.");
//    //}
//
//    expr.print_simple();
//}

fn main() {
    let s = "4071272 * w + 3590309086 * (w | z & z & (y |
        w)) + 3690425492 * (w & x ^ ~z ^ ~y) + 3735539420 *
        (y ^ (w & z | y)) + 3176111544 * ((x & y | x ^ z) ^
        ~y) + 90227856 * (y & x & x & (x & w | ~w)) +
        2231609705 * (~z ^ w & x | z ^ x | x | w) +
        263864489 * (w ^ z | x | w | z) + 17029904 *
        (~w ^ w ^ (x | z) ^ ~x) + 2987805917 * (z & x) +
        1280785607 * z + 2092481170 * (y & (w & y | ~z)) +
        637019375 * (~w | w & z & ~x) + 3221225472 * ((x |
        w) ^ x ^ x ^ y) + 3985988880 * x + 263864489 * (~~w &
        x) + 469200020 * ((z ^ w & w) & ~(x ^ y)) +
        1774328762 * ((w | x) & (x ^ z) & z) + 3645311564 *
        (~(z | w) | w) + 3194849700 * ~((y | y) ^ y ^ z)
        + 1678283628 * ~(~y & ~w) + 1083630375 * y";
    //let e = LUExpr::from_string("(x ^ y) + 2 * (x & y)").unwrap();
    let e = LUExpr::from_string(s).unwrap();
    let d = deobfuscate_linear(e, 32, false);
    println!("{}", d);
}

fn deobfuscate_linear(e: LUExpr, bits: u32, fast: bool) -> LUExpr {
    // Get all the variables.
    let vars = e.vars();

    // Insert some default operations.
    let mut ops = Vec::new();
    ops.push(UExpr::Ones);
    for v in &vars {
        ops.push(UExpr::Var(v.clone()));
    }
    for v in &vars {
        for w in &vars {
            if std::ptr::eq(v, w) {
                break;
            }

            ops.push(UExpr::and(UExpr::Var(w.clone()), UExpr::Var(v.clone())));
            ops.push(UExpr::or(UExpr::Var(w.clone()), UExpr::Var(v.clone())));
            ops.push(UExpr::xor(UExpr::Var(w.clone()), UExpr::Var(v.clone())));
        }
    }

    // Insert the original operations so we always have a solution.
    for u in &e.0 {
        if !ops.contains(&u.1) {
            ops.push(u.1.clone());
        }
    }

    //for o in &ops {
    //    println!("{o}");
    //}

    if fast {
        let ops: Vec<_> = ops.iter().map(|e| LUExpr::from(e.clone())).collect();
        return rewrite(&e, &ops, bits, false).unwrap_or(e);
    }

    assert!(vars.len() <= 63, "More than 63 variables are currently not supported \
            (You wouldn't be able to run this anyways).");

    let mut val = Valuation::zero(vars.clone());

    let rows = 1 << vars.len();
    let cols = ops.len();

    let mut a = Matrix::zero(rows, cols);
    let mut b = Vector::zero(rows);

    // Build up the matrix.
    for i in 0..rows {
        let row = a.row_mut(i);

        // Initialize the valuation.
        for (j, c) in vars.iter().enumerate() {
            val[c] = -Integer::from((i >> j) & 1);
        }

        // println!("{:?}", val);

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = e.eval(&val).keep_bits(bits);
        }

        // Write the desired result into the vector.
        b[i] = e.eval(&val).keep_bits(bits);
    }

    // Solve the system.
    let mut l = diophantine::solve_modular(&a, &b, &(Integer::from(1) << bits));

    // If I did everything correctly, this system should always have a solution.
    assert!(!l.is_empty());
    assert!(l.offset.dim() == ops.len());

    let complexity: Vec<_> = ops.iter()
        .map(|e| e.complexity())
        .collect();

    for v in l.lattice.basis.rows_mut() {
        for (e, c) in v.iter_mut().zip(&complexity) {
            *e *= c;
        }
    }

    l.lattice.lll(&rug::Rational::from((99, 100)));

    for (e, c) in l.offset.iter_mut().zip(&complexity) {
        *e *= c;
    }

    let mut solution = l.offset.clone();
    let norm = solution.norm_sqr().to_f64().sqrt();
    solution -= &l.lattice.cvp_planes_f64(l.offset.view(), Some(norm))
        .unwrap();

    for (e, c) in solution.iter_mut().zip(&complexity) {
        debug_assert!(e.is_divisible_u(*c));
        e.div_exact_u_mut(*c);
    }

    // Put it in a LUExpr.
    let v: Vec<_> = solution.iter()
        .cloned()
        .map(|i| i.keep_signed_bits(bits))
        .zip(ops)
        .collect();

    LUExpr(v)
}

#[cfg(feature = "z3")]
fn deobfuscate_nonlinear(e: Rc<Expr>) {
    let mut v = Vec::new();
    let mut q = std::collections::VecDeque::new();
    q.push_back(e);
    while let Some(e) = q.pop_front() {
        match e.as_ref() {
            Expr::Add(l, r) => {
                q.push_back(l.clone());
                q.push_back(r.clone());
            },
            Expr::Mul(l, r) => {
                if let Expr::Const(i) = l.as_ref() {
                    v.push((i.clone(), r.clone()));
                } else if let Expr::Const(i) = r.as_ref() {
                    v.push((i.clone(), l.clone()));
                } else {
                    println!("Could not deobfuscate expression, \
                        because a term in the linear combination \
                        is not in a form we recognize.");
                }
            },
            _ => {
                println!("Could not deobfuscate the expression, \
                    because it is not in a form we recognize.");
                return;
            }
        }
    }
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
    /// Integer width.
    bits: u32,

    /// The operations used for rewriting things.
    ops: Vec<LUExpr>,

    /// The variables occurring in the expression.
    vars: Vec<String>,

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

        let mut x = Rc::new(Expr::Var(
            o.vars.get(0).map_or_else(|| "x".to_owned(), |v| v.clone())));
        let mut y = Rc::new(Expr::Var(
            o.vars.get(0).map_or_else(|| "y".to_owned(), |v| v.clone())));
        let mut z = Rc::new(Expr::Var(
            o.vars.get(0).map_or_else(|| "z".to_owned(), |v| v.clone())));

        e.as_ref().clone().substitute("X", &mut x)
            .substitute("Y", &mut y)
            .substitute("Z", &mut z)
    }

    fn obfuscate_impl(&mut self, e: &mut Rc<Expr>) {
        if self.rng.gen::<f32>() < 0.9 {
            match Rc::make_mut(e) {
                Expr::Const(i) => {
                    let u = LUExpr::constant(i.clone());
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr().into();
                    }
                },
                Expr::Var(v) => {
                    let u = LUExpr::var("X".to_owned());
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", &mut Rc::new(Expr::Var(v.clone())))
                            .into();
                    }
                },
                Expr::Add(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X+Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", l)
                            .substitute("Y", r)
                            .into();
                    }
                },
                Expr::Sub(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X-Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", l)
                            .substitute("Y", r)
                            .into();
                    }
                },
                Expr::Mul(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let f = Expr::from_string("(X&Y)*(X|Y)+(X&~Y)*(~X&Y)")
                        .unwrap();

                    *e = f.substitute("X", l).substitute("Y", r).into();
                },
                Expr::Div(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                },
                Expr::Neg(i) => {
                    self.obfuscate_impl(i);
                    let u = LUExpr::from_string("-X").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", i)
                            .into();
                    }
                },
                Expr::And(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X&Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", l)
                            .substitute("Y", r)
                            .into();
                    }
                },
                Expr::Or(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X|Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", l)
                            .substitute("Y", r)
                            .into();
                    }
                },
                Expr::Xor(l, r) => {
                    self.obfuscate_impl(l);
                    self.obfuscate_impl(r);
                    let u = LUExpr::from_string("X^Y").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", l)
                            .substitute("Y", r)
                            .into();
                    }
                },
                Expr::Not(i) => {
                    self.obfuscate_impl(i);
                    let u = LUExpr::from_string("~X").unwrap();
                    if let Some(f) = rewrite(&u, &self.ops, self.bits, true) {
                        *e = f.to_expr()
                            .substitute("X", i)
                            .into();
                    }
                },
            };
        }

        if self.rng.gen::<f32>() < 0.3 {
            let degree = Uniform::from(2..4).sample(&mut self.rng);
            let (p, q) = perm_poly::perm_pair(&self.qr, degree);
            let p = p.to_expr().substitute("x", e);
            *e = Rc::new(q.to_expr().substitute("x", &mut Rc::new(p)));
        }
    }

    fn init(e: &Expr, bits: u32) -> Self {
        // Get the variables in the expression.
        let vars = e.vars();

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
        let qr = perm_poly::ZeroIdeal::init(bits);
        let rng = rand::thread_rng();

        Self { bits, ops, vars, qr, rng }
    }
}

/// Hacky function to convert a rug::Integer to a z3 bitvector.
#[cfg(feature = "z3")]
pub(crate) fn int_to_bv<'ctx>(
    ctx: &'ctx z3::Context, width: u32, i: &Integer
) -> z3::ast::BV<'ctx> {
    use z3::ast::Ast;
    let mut bits = i.to_digits::<bool>(rug::integer::Order::Lsf);
    bits.resize(width as usize, false);
    let ast = unsafe {
        z3_sys::Z3_mk_bv_numeral(
            *(ctx as *const _ as *const z3_sys::Z3_context),
            width,
            bits.as_ptr()
        )
    };
    z3::ast::BV::new(ctx, ast)
}

/// Very hacky macro that acts as a match statement for types.
/// ```
/// use crate::select;
/// let i = select!(f32,
///     f32 => { 1 },
///     f64 => { 2 },
///     default => { 3 },
/// );
/// assert_eq!(i, 1);
/// ```
/// This is used to generate different code for different types.
/// Ideally we would use traits, but I felt I would run into too many
/// problems. We want to support `Copy` types, like `f32`, and types
/// like `rug::Float` which own memory.
/// That means we want to use Copy for f64 where possible
/// instead of passing around references, but for rug::Float
/// we want to avoid copying at all costs.
/// In addition operations like &rug::Float + &rug::Float
/// return "Incomplete computation-values",
/// which you have to `.complete` or assign in some way.
/// I felt like doing all this with traits would be
/// challenging, but obviously the superior solution.
/// Instead we use this very hacky macro.
/// A proc macro would be preferable, but this does
/// the job for now.
macro_rules! select {
    (f32, f32 => { $($b:tt)* }, $($o:tt)*) => { $($b)* };
    (f64, f64 => { $($b:tt)* }, $($o:tt)*) => { $($b)* };
    (Float, Float => { $($b:tt)* }, $($o:tt)*) => { $($b)* };
    (Integer, Integer => { $($b:tt)* }, $($o:tt)*) => { $($b)* };
    (Rational, Rational => { $($b:tt)* }, $($o:tt)*) => { $($b)* };
    ($t:tt, default => {$($b:tt)* }, $($o:tt)*) => {
        $($b)*
    };
    ($t:tt, $n:ty => { $($nb:tt)* }, $($o:tt)*) => {
        select!($t, $($o)*)
    };
}

pub(crate) use select;