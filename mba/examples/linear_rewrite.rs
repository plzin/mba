use mba::bitwise_expr::LBExpr;
use mba::linear_mba::rewrite;
use mba::rings::U8;

// Rewrite an expression using linear MBA.
fn main() {
    let mut rng = rand::rng();
    let r = &U8;
    let expr = LBExpr::from_string("x+y".to_owned(), r).unwrap();
    let ops = &[
        LBExpr::from_string("x&y".to_owned(), r).unwrap(),
        LBExpr::from_string("x^y".to_owned(), r).unwrap(),
    ];

    match rewrite(&expr, ops, Some(&mut rng), r) {
        None => println!("Can't rewrite."),
        Some(e) => println!("{e}"),
    }
}
