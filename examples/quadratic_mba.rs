use itertools::Itertools;
use rug::Integer;
use mba::uniform_expr::{UExpr, LUExpr};
use mba::linear_mba;
use mba::expr::{Expr, ExprOp};


// Rewrite a linear MBA expression using quadratic MBA.
fn main() {
    let Some((_, bits, expr)) = std::env::args().next_tuple() else {
        println!("Arguments: <bits> <expr>");
        return;
    };

    let Ok(bits) = bits.parse() else {
        println!("Invalid number of bits.");
        return;
    };

    let Some(expr) = LUExpr::from_string(expr) else {
        println!("Invalid expression.");
        return;
    };

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

    // Contains the linear combination of the multiplication.
    let mut v: Vec<(Integer, UExpr, UExpr)> = Vec::new();

    // Obfuscate each coefficient.
    for (c, e) in &obf.0 {
        // Rewrite the coefficient.
        let coeff = LUExpr::constant(c.clone());
        let mut ops = ops.to_vec();
        ops.retain(|f| &f.0[0].1 != e);
        let coeff = linear_mba::rewrite(&coeff, &ops, bits, true).unwrap();

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

    expr.simplify();

    println!("{expr}");
}
