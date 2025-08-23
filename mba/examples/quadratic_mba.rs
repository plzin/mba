use itertools::Itertools;
use mba::{
    bitwise_expr::{BExpr, LBExpr},
    expr::{Expr, ExprOp},
    formatter::Formatter,
    linear_mba::rewrite,
    rings::{
        BinaryBigInt, BinaryRing, IntDivRing, OrderedRing, RingElement as _,
        U8, U16, U32, U64, U128,
    },
};

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

    match bits {
        8 => quadratic_mba(expr, &U8),
        16 => quadratic_mba(expr, &U16),
        32 => quadratic_mba(expr, &U32),
        64 => quadratic_mba(expr, &U64),
        128 => quadratic_mba(expr, &U128),
        _ => quadratic_mba(expr, &BinaryBigInt::new(bits)),
    }
}

fn quadratic_mba<R: BinaryRing + OrderedRing + IntDivRing>(
    expr: String,
    r: &R,
) {
    let expr = match LBExpr::from_string(expr, r) {
        Ok(expr) => expr,
        Err(e) => {
            println!("Invalid expression: {e}");
            return;
        },
    };

    let make_op = |s: &str| LBExpr::from_string(s.to_owned(), r).unwrap();

    // The operations for rewriting.
    let ops = &[
        LBExpr::from(BExpr::Ones),
        make_op("x^y"),
        make_op("x&y"),
        make_op("~(x|y)"),
        make_op("~x"),
        make_op("~y"),
        make_op("z"),
        make_op("x"),
        make_op("y"),
        make_op("y&z"),
        make_op("x|z"),
        make_op("~x&z"),
        make_op("y|~z"),
    ];

    // Rewrite the constant.
    let rng = &mut rand::rng();
    let obf = rewrite(&expr, ops, Some(rng), r).unwrap();

    // Contains the linear combination of the multiplication.
    let mut v: Vec<(R::Element, BExpr, BExpr)> = Vec::new();

    // Obfuscate each coefficient.
    for (c, e) in &obf.0 {
        // Rewrite the coefficient.
        let coeff = LBExpr::constant(c.clone(), r);
        let mut ops = ops.to_vec();
        ops.retain(|f| &f.0[0].1 != e);
        let coeff = rewrite(&coeff, &ops, Some(rng), r).unwrap();

        // Add all coefficients of the obfuscated coefficient
        // to the linear combination.
        for (d, f) in coeff.0 {
            let e = e.clone();
            let (fst, snd) = if e <= f {
                (e, f)
            } else {
                (f, e)
            };

            match v.iter_mut().find(|(_, a, b)| &fst == a && &snd == b) {
                Some((c, _, _)) => r.add_assign(c, &d),
                None => v.push((d, fst, snd)),
            }
        }
    }

    // Convert a summand to an expression.
    let term_to_expr = |c: &R::Element, e: &BExpr, f: &BExpr| {
        Expr::new(ExprOp::Mul(
            Expr::new(ExprOp::Const(c.clone())),
            Expr::new(ExprOp::Mul(
                Expr::new(e.to_expr(r)),
                Expr::new(f.to_expr(r)),
            )),
        ))
    };

    // Convert the linear combination to an `Expr`.
    let mut iter = v.iter().filter(|(c, _, _)| !c.is_zero());
    let mut expr = match iter.next() {
        Some((c, e, f)) => term_to_expr(c, e, f),
        None => Expr::new(ExprOp::Const(R::zero())),
    };

    for (c, e, f) in iter {
        let a = term_to_expr(c, e, f);
        expr = Expr::new(ExprOp::Add(expr, a));
    }

    expr.simplify();

    println!("{}", expr.display(Formatter::C, 0, "", r));
}
