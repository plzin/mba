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
use itertools::Itertools;

#[cfg(feature = "z3")]
use z3::{self, ast::Ast};

// It would be nicer to import the symbol_table crate ourselves,
// and re-export GlobalSymbol, but we don't want to run into
// version conflicts.
pub use egg::Symbol;
use std::rc::Rc;

mod matrix;
use matrix::*;

mod vector;
use vector::*;

mod valuation;

mod expr;
use expr::*;

mod uniform_expr;
use uniform_expr::*;

use crate::linear_mba::{DeobfuscationConfig, SolutionAlgorithm};

mod diophantine;
mod lattice;
mod poly;
mod perm_poly;
mod linear_mba;
mod simplify_boolean;

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
    let Some((_, bits, num)) = std::env::args().next_tuple() else {
        println!("\
            The current `main` is just one example of how to use the crate.\n\
            Eventually, this crate should be a library.\n\
            The current `main` obfuscates an n bit constant using non-linear MBA.\n\
            Usage: mba <bits> <constant>\
        ");
        return;
    };

    let Ok(bits) = bits.parse() else {
        println!("Invalid number of bits.");
        return;
    };

    // This is the constant that will be obfuscated.
    let Ok(input) = Integer::from_str_radix(&num, 10) else {
        println!("Invalid constant.");
        return;
    };

    // Initialize stuff.
    let n = Integer::from(1) << bits;
    let zi = perm_poly::ZeroIdeal::init(bits);

    // Obfuscate a different constant and apply a
    // polynomial at the end to get the real one.
    //let p = poly::Poly::from_vec(vec![0, 1, 2]);
    //let q = perm_poly::compute_inverse(&p, &zi);
    //let input = q.eval_bits(&input, bits);

    // The expression to obfuscate.
    let expr = LUExpr::constant(input.clone());

    // The operations for rewriting.
    let ops = &[
        LUExpr::from(UExpr::Ones),
        LUExpr::from_string("x^y").unwrap(),
        LUExpr::from_string("x&y").unwrap(),
        LUExpr::from_string("~(x|y)").unwrap(),
        LUExpr::from_string("~x").unwrap(),
        LUExpr::from_string("~y").unwrap(),
        LUExpr::from_string("z").unwrap(),
        LUExpr::from_string("x").unwrap(),
        LUExpr::from_string("y").unwrap(),
        LUExpr::from_string("y&z").unwrap(),
        LUExpr::from_string("x|z").unwrap(),
        LUExpr::from_string("~x&z").unwrap(),
        LUExpr::from_string("y|~z").unwrap(),
    ];

    // Rewrite the constant.
    let obf = linear_mba::rewrite(&expr, ops, bits, true).unwrap();
    println!("{obf}");

    // Contains the linear combination of the multiplication.
    let mut v: Vec<(Integer, UExpr, UExpr)> = Vec::new();

    // Obfuscate each coefficient.
    for (c, e) in &obf.0 {
        // Rewrite the coefficient.
        let coeff = LUExpr::constant(c.clone());
        let mut ops = ops.to_vec();
        ops.retain(|f| &f.0[0].1 != e);
        let coeff = linear_mba::rewrite(&coeff, &ops, bits, true).unwrap();
        println!("{e}: {coeff}");

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
        Expr::new(ExprOp::Mul(
            Expr::new(ExprOp::Const(c.clone())),
            Expr::new(ExprOp::Mul(Expr::new(e.to_expr()), Expr::new(f.to_expr())))
        ))
    };

    // Convert the linear combination to an `Expr`.
    let mut iter = v.iter().filter(|(c, _, _)| c != &Integer::ZERO);
    let mut expr = match iter.next() {
        Some((c, e, f)) => term_to_expr(c, e, f),
        None => Expr::new(ExprOp::Const(Integer::new())),
    };

    for (c, e, f) in iter {
        let a = term_to_expr(c, e, f);
        expr = Expr::new(ExprOp::Add(expr, a));
    }

    //let mut re = Rc::new(expr);
    //expr = p.to_expr().substitute("x", &mut re);


    //if cfg!(feature = "z3") {
    //    let cfg = z3::Config::new();
    //    let ctx = z3::Context::new(&cfg);
    //    let bv = expr.to_z3_bv(&ctx, bits);

    //    let solver = z3::Solver::new(&ctx);
    //    let assertion = bv._eq(&int_to_bv(&ctx, bits, &input)).not();
    //    println!("assertion: {}", assertion);
    //    solver.assert(&assertion);
    //    println!("Using z3 to check that the rewritten expression is equivalent.");
    //    let now = std::time::Instant::now();
    //    let result = solver.check();
    //    let dur = now.elapsed().as_secs_f32();
    //    match result {
    //        z3::SatResult::Sat => println!("Found counterexample:\n{}",
    //            solver.get_model().unwrap()),
    //        z3::SatResult::Unsat => println!("Verification successful."),
    //        z3::SatResult::Unknown => println!("z3 aborted."),
    //    }
    //    println!("z3 verification took {dur:.2} seconds.");
    //}

    println!("{expr}");
}

// fn main() {
//     env_logger::init();
//
//     let s = "4071272 * w + 3590309086 * (w | z & z & (y |
//         w)) + 3690425492 * (w & x ^ ~z ^ ~y) + 3735539420 *
//         (y ^ (w & z | y)) + 3176111544 * ((x & y | x ^ z) ^
//         ~y) + 90227856 * (y & x & x & (x & w | ~w)) +
//         2231609705 * (~z ^ w & x | z ^ x | x | w) +
//         263864489 * (w ^ z | x | w | z) + 17029904 *
//         (~w ^ w ^ (x | z) ^ ~x) + 2987805917 * (z & x) +
//         1280785607 * z + 2092481170 * (y & (w & y | ~z)) +
//         637019375 * (~w | w & z & ~x) + 3221225472 * ((x |
//         w) ^ x ^ x ^ y) + 3985988880 * x + 263864489 * (~~w &
//         x) + 469200020 * ((z ^ w & w) & ~(x ^ y)) +
//         1774328762 * ((w | x) & (x ^ z) & z) + 3645311564 *
//         (~(z | w) | w) + 3194849700 * ~((y | y) ^ y ^ z)
//         + 1678283628 * ~(~y & ~w) + 1083630375 * y";
//     //let e = LUExpr::from_string("(x ^ y) + 2 * (x & y)").unwrap();
//     let cfg = DeobfuscationConfig {
//         alg: SolutionAlgorithm::LeastComplexTerms,
//         boolean: true,
//     };
//     let e = LUExpr::from_string(s).unwrap();
//     let d = linear_mba::deobfuscate_luexpr(
//         e, 32, &cfg
//     );
//     println!("{}", d);
// }

// Solve a system of linear congruences.
// fn main() {
//     let solution = diophantine::solve_modular(
//         Matrix::from_array([[3, 5], [4, 2]]).view(),
//         Vector::from_slice(&[2.into(), 0.into()]),
//         &256.into()
//     );
//
//     if (solution.is_empty()) {
//         println!("No solution");
//     } else {
//         println!("Off: {:?}\nBasis:", solution.offset);
//         for b in solution.lattice.basis.rows() {
//             println!("{:?}", b);
//         }
//     }
// }

//fn main() {
//    // The example will obfuscate x+y where x and y are 32-bit integers.
//    // This should be about the same as the WASM obfuscation.
//    // Note that the example obfuscation is quite limited and only meant as
//    // a starting point and an example of how to combine the primitives.
//    let mut e = Expr::from_string("x+y").unwrap();
//    let cfg = linear_mba::ObfuscationConfig::default();
//    linear_mba::obfuscate(&mut e, 32, &cfg);
//    println!("{e}");
//}

/// Hacky function to convert a rug::Integer to a z3 bitvector.
#[cfg(feature = "z3")]
pub(crate) fn int_to_bv<'ctx>(
    ctx: &'ctx z3::Context, width: u32, i: &Integer
) -> z3::ast::BV<'ctx> {
    use z3::ast::Ast;
    let mut bits = i.to_digits::<bool>(rug::integer::Order::Lsf);
    bits.resize(width as usize, false);
    unsafe {
        let ast = z3_sys::Z3_mk_bv_numeral(
            *(ctx as *const _ as *const z3_sys::Z3_context),
            width,
            bits.as_ptr()
        );
        z3::ast::BV::wrap(ctx, ast)
    }
}

/// Used internally to convert an integer from and iterator over characters.
pub(crate) fn int_from_it(
    it: &mut std::iter::Peekable<std::str::Chars>
) -> Option<Integer> {
    let c = it.next()?;
    let mut num: Integer = c.to_digit(10)?.into();

    let mut lookahead = it.peek();
    while let Some(c) = lookahead {
        if !c.is_ascii_digit() {
            break;
        }

        num *= 10;
        num += c.to_digit(10).unwrap();
        it.next();
        lookahead = it.peek();
    }

    Some(num)
}