use mba::{expr::Expr, linear_mba, rings::BinaryRing};
use wasm_bindgen::prelude::*;

use crate::Formatter;

/// Obfuscates a general expression.
#[wasm_bindgen(js_name = "obfuscate")]
pub fn obfuscate(
    expr: String,
    bits: u32,
    #[wasm_bindgen(js_name = "auxiliaryVars")] auxiliary_vars: usize,
    #[wasm_bindgen(js_name = "rewriteExprDepth")] rewrite_expr_depth: usize,
    #[wasm_bindgen(js_name = "rewriteExprCount")] rewrite_expr_count: usize,
    formatter: Formatter,
) -> Result<String, String> {
    mba::choose_binary_ring!(
        obfuscate_impl(
            expr,
            auxiliary_vars,
            rewrite_expr_depth,
            rewrite_expr_count,
            formatter,
            &r
        ),
        r = bits
    )
}

pub fn obfuscate_impl<R: BinaryRing>(
    expr: String,
    auxiliary_vars: usize,
    rewrite_expr_depth: usize,
    rewrite_expr_count: usize,
    formatter: Formatter,
    r: &R,
) -> Result<String, String> {
    let mut e = Expr::from_string(expr, r)
        .map_err(|e| format!("Failed to parse expression: {e}"))?;

    let cfg = linear_mba::ObfuscationConfig {
        auxiliary_vars,
        rewrite_expr_depth,
        rewrite_expr_count,
        ..Default::default()
    };

    linear_mba::obfuscate(&mut e, &cfg, &mut rand::rng(), r);
    Ok(e.display_function(formatter.to_rust(), "f", r).to_string())
}
