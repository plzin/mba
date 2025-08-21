use crate::Formatter;
use mba::{bitwise_expr::LBExpr, linear_mba, rings::BinaryRing};
use wasm_bindgen::prelude::*;

/// Obfuscates a linear MBA expression `expr` using the given rewrite operations
/// `rewrite_ops` for the given integer width `bits`.
///
/// If `randomize` is true, the resulting linear combination will be chosed
/// randomly from the lattice of solutions. Otherwise, the particular solution
/// found by the solver will be used.
///
/// The result will be formatted using the given `formatter`.
///
/// This function can fail for the following reasons:
/// - `expr` cannot be parsed into an [`LBExpr`].
/// - Any of the `rewrite_ops` cannot be parsed into an [`LBExpr`].
/// - There is no solution to the linear system.
/// - `bits` is 0.
#[wasm_bindgen(js_name = "obfuscateLinear")]
pub fn obfuscate_linear(
    expr: String,
    #[wasm_bindgen(js_name = "rewriteOps")] rewrite_ops: Vec<String>,
    bits: u32,
    randomize: bool,
    formatter: Formatter,
) -> Result<String, String> {
    mba::choose_binary_ring!(
        obfuscate_linear_impl(expr, rewrite_ops, randomize, formatter, &r),
        r = bits
    )
}

fn obfuscate_linear_impl<R: BinaryRing>(
    expr: String,
    rewrite_ops: Vec<String>,
    randomize: bool,
    formatter: Formatter,
    r: &R,
) -> Result<String, String> {
    // Parse the expression.
    let expr =
        LBExpr::from_string(expr, r).map_err(|e| format!("Failed to parse expression: {e}"))?;

    // Parse the rewrite operations.
    let rewrite_ops: Vec<_> = rewrite_ops
        .into_iter()
        .map(|op| {
            // A bit wasteful to clone the `op` here,
            // but we potentially need it on the error message.
            LBExpr::from_string(op.clone(), r)
                .map_err(|e| format!("Failed to parse rewrite operation {op}: {e}"))
        })
        .try_collect()?;

    // Try to rewrite the expression.
    let e = linear_mba::rewrite(&expr, &rewrite_ops, randomize.then(rand::rng).as_mut(), r)
        .ok_or("There is no solution")?;

    // Format the result.
    Ok(e.display_function(formatter.to_rust(), "f", r).to_string())
}
