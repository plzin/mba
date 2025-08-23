//! Functions related to rewriting an [`LBExpr`] using a set of [`LBExpr`]s.

use rand::Rng;

use crate::{
    Symbol,
    bitwise_expr::LBExpr,
    lattice::AffineLattice,
    matrix::OwnedMatrix,
    rings::BinaryRing,
    solver::solve_via_modular_diagonalize,
    valuation::Valuation,
    vector::{OwnedVector, VectorView},
};

/// Rewrite an [`LBExpr`] using a set of [`LBExpr`]s.
///
/// If `rng` is not `None`, a random solution among all possible solutions will
/// be sampled. Otherwise, the solution is the deterministic solution returned
/// by the solver, but it is not necessarily the smallest solution.
pub fn rewrite<R: BinaryRing, Rand: Rng>(
    expr: &LBExpr<R>,
    ops: &[LBExpr<R>],
    rng: Option<&mut Rand>,
    ring: &R,
) -> Option<LBExpr<R>> {
    // Find all variables we have access to.
    // This includes variables in the expression as well as potentially the ops.
    let mut v = std::collections::BTreeSet::new();
    expr.vars_impl(&mut v);
    ops.iter().for_each(|e| e.vars_impl(&mut v));
    let v: Vec<_> = v.into_iter().collect();

    // Solve the system.
    let l = solve_linear_system(expr, ops, &v, ring);

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = match rng {
        Some(rng) => l.sample_point(rng, ring),
        None => l.offset,
    };

    // Collect the solution into an LBExpr.
    Some(collect_solution(solution.view(), ops, ring))
}

/// Solves the linear system for rewriting `expr` using the rewrite operations
/// `ops`.
///
/// `vars` is all the variables that occur in `expr` or `ops`.
pub fn solve_linear_system<R: BinaryRing>(
    expr: &LBExpr<R>,
    ops: &[LBExpr<R>],
    vars: &[Symbol],
    ring: &R,
) -> AffineLattice<R> {
    assert!(
        vars.len() < usize::BITS as usize,
        "More than {} variables are currently not supported on your system.",
        usize::BITS - 1
    );

    // Allocate a valuation.
    let mut val = Valuation::vars_zero_or_panic(vars);

    let rows = 1 << vars.len();
    let cols = ops.len();

    let mut a = OwnedMatrix::<R>::zero(rows, cols);
    let mut b = OwnedVector::zero(rows);

    // Build up the matrix.
    for i in 0..rows {
        let row = a.row_mut(i);

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = e.eval(&mut val, ring);
        }

        // Write the desired result into the vector.
        b[i] = expr.eval(&mut val, ring);

        val.inc_bin(ring);
    }

    solve_via_modular_diagonalize(a, b, ring)
}

/// Converts a solution to the linear system and the rewrite operations to a
/// single [`LBExpr`] with no duplicate [`crate::bitwise_expr::BExpr`]s.
pub fn collect_solution<R: BinaryRing>(
    solution: &VectorView<R>,
    ops: &[LBExpr<R>],
    r: &R,
) -> LBExpr<R> {
    // Put it in a LBExpr.
    let mut l = LBExpr::zero();
    for (c, o) in solution.iter().zip(ops.iter()) {
        for (d, e) in &o.0 {
            // Is the BExpr already in the linear combination?
            match l.0.iter_mut().find(|(_, f)| f == e) {
                Some((f, _)) => r.mul_add_assign(f, c, d),
                None => l.0.push((r.mul(c.clone(), d), e.clone())),
            }
        }
    }

    l.remove_zero_terms();
    l
}

#[test]
fn rewrite_test() {
    use rand::{SeedableRng as _, rngs::StdRng};

    use crate::rings::U8;

    let mut rng = StdRng::seed_from_u64(0);
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
