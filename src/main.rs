use rug::Integer;

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

/// Rewrite an expression using a set of operations modulo n.
fn rewrite(expr: &LUExpr, ops: &[UniformExpr], n: &Integer) -> Option<LUExpr> {
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
        .map(|(m, e)| (m.clone(), e.clone()))
        .collect();

    Some(LUExpr(v))
}

// fn main() {
//     let r = perm_poly::QuotientRing::init(8);
//     let (p, q) = perm_poly::perm_pair(&r, 3);
//     println!("p(x) = {}", p);
//     println!("q(x) = {}", q);
//     //perm_poly::check_inverse_64();
// }

fn main() {
    //let e = Expr::from_string("~x&y+z");
    //println!("{:?}", e);
    // let e = LUExpr::from_string("~x&y|z+123");
    // println!("{:?}", e);

    let n = 32;

    let expr = LUExpr::from_string("x+y").unwrap();
    // let expr = Expr::from_string("128").unwrap();

    let ops = [
        UniformExpr::from_string("x&y").unwrap(),
        UniformExpr::from_string("x^y").unwrap(),
        UniformExpr::from_string("~(x^y^z)").unwrap(),
        UniformExpr::from_string("~x").unwrap(),
        UniformExpr::from_string("~y").unwrap(),
        UniformExpr::from_string("y").unwrap(),
        UniformExpr::from_string("z").unwrap(),
        UniformExpr::from_string("y&z").unwrap(),
        UniformExpr::from_string("x|z").unwrap(),
        UniformExpr::from_string("(~x)&z").unwrap(),
        UniformExpr::from_string("y|(~z)").unwrap(),
    ];

    let e = rewrite(&expr, &ops, &Integer::new().set_bit(n, true));
    if let Some(e) = e {
        println!("{}", e);
    } else {
        println!(
            "Rewriting the expression with the \
            given operations isn't possible.");
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
