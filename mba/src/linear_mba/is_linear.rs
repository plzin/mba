//! Functions for checking if an [`Expr`] implements a linear MBA function.

use rand::Rng;

use crate::{Symbol, expr::Expr, rings::BinaryRing, valuation::Valuation};

/// Checks if an [`Expr`] implements a linear MBA function by checking if the
/// fundamental theorem of linear MBA holds for all inputs.
///
/// <div class="warning">
/// This checks *ALL* inputs, which for almost all cases is way too expensive.
/// This function should only be used for testing.
/// </div>
pub fn is_linear_mba<R: BinaryRing>(expr: &Expr<R>, ring: &R) -> bool {
    let (vars, b, mut val) = build_table(expr, ring);

    // Loop over all possible values of the variables.
    while {
        // Does this follow the fundamental theorem of linear MBA?
        let actual = expr.eval(&mut val, ring);
        let expected = compute_expected_result(&vars, &b, &mut val, ring);
        if actual != expected {
            return false;
        }

        !val.inc_full(ring)
    } {}

    true
}

/// Checks if an [`Expr`] implements a linear MBA function by randomly sampling
/// `num_inputs` inputs and checking if the fundamental theorem of linear MBA
/// holds for those inputs.
pub fn is_probably_linear_mba<R: BinaryRing, Rand: Rng>(
    expr: &Expr<R>,
    num_inputs: usize,
    rng: &mut Rand,
    ring: &R,
) -> bool {
    let (vars, b, mut val) = build_table(expr, ring);

    for _ in 0..num_inputs {
        // Update the valuation to random values.
        val.update_random(rng, ring);

        // Does this follow the fundamental theorem of linear MBA?
        let actual = expr.eval(&mut val, ring);
        let expected = compute_expected_result(&vars, &b, &mut val, ring);
        if actual != expected {
            return false;
        }
    }

    true
}

/// Computes the table of values of an [`Expr`] for all 0/-1 inputs.
fn build_table<R: BinaryRing>(
    expr: &Expr<R>,
    ring: &R,
) -> (Vec<Symbol>, Vec<R::Element>, Valuation<R>) {
    let vars = expr.vars();
    assert!(
        vars.len() < usize::BITS as usize,
        "More than {} variables are currently not supported on your system.",
        usize::BITS - 1
    );

    // Allocate a valuation.
    let mut val = Valuation::vars_zero_or_panic(&vars);

    let rows = 1 << vars.len();

    let mut b = vec![R::zero(); rows];

    // Evaluate the expression for all -1/0 inputs.
    // Store the results of the expression for all inputs.
    for b in &mut b {
        // Write the desired result into the vector.
        *b = expr.eval(&mut val, ring);

        // Increment the valuation.
        val.inc_bin(ring);
    }

    (vars, b, val)
}

/// Computes the expected result based on the fundamental theorem of linear MBA.
fn compute_expected_result<R: BinaryRing>(
    vars: &[Symbol],
    b: &[R::Element],
    val: &mut Valuation<R>,
    r: &R,
) -> R::Element {
    (0..r.bits()).fold(R::zero(), |acc, j| {
        r.sub(acc, &r.shl(b[compute_index(vars, val, j, r)].clone(), j))
    })
}

/// This function is only meaningful in the context of the fundamental theorem
/// and [`build_table`]. In the fundamental theorem, for each bit index `j`,
/// we need to evaluate the expression at the inputs 0/-1 depending on whether
/// but `j` is set or not in each input variable. [`build_table`] evaluates the
/// expression at all such inputs and stores the results in a table. This
/// function computes the index into that table in this scenario.
fn compute_index<R: BinaryRing>(
    vars: &[Symbol],
    val: &mut Valuation<R>,
    j: u32,
    r: &R,
) -> usize {
    vars.iter().enumerate().fold(0, |acc, (i, &v)| {
        acc | (R::bit(val.value(v, r), j) as usize) << i
    })
}
