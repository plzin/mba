use mba::uniform_expr::LUExpr;
use mba::linear_mba::rewrite;

// Rewrite an expression using linear MBA.
fn main() {
    let expr = LUExpr::from_string("x+y").unwrap();
    let ops = &[
        LUExpr::from_string("x&y").unwrap(),
        LUExpr::from_string("x^y").unwrap(),
    ];

    match rewrite(&expr, ops, 8, true) {
        None => println!("Can't rewrite."),
        Some(e) => println!("{e}"),
    }
}