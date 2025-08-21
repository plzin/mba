use mba::{
    bitwise_expr::LBExpr,
    choose_binary_ring,
    linear_mba::{self, DeobfuscationConfig},
    rings::BinaryRing,
};
use wasm_bindgen::prelude::*;

use crate::{Formatter, SolutionAlgorithm};

/// TODO: Add a `allow_non_lbexpr` parameter so you can pass in non-LBExprs that
/// are essentially asserted to be equivalent to an LBExpr.
#[wasm_bindgen(js_name = "deobfuscateLinear")]
pub fn deobfuscate_linear(
    expr: String,
    bits: u32,
    alg: SolutionAlgorithm,
    #[wasm_bindgen(js_name = "detectBoolean")] detect_boolean: bool,
    formatter: Formatter,
) -> Result<String, String> {
    let cfg = DeobfuscationConfig {
        alg: alg.to_rust(),
        boolean: detect_boolean,
    };
    choose_binary_ring!(deobfuscate_linear_impl(expr, &cfg, formatter, &r), r = bits)
}

fn deobfuscate_linear_impl<R: BinaryRing>(
    expr: String,
    cfg: &DeobfuscationConfig,
    formatter: Formatter,
    ring: &R,
) -> Result<String, String> {
    // Parse the expression.
    let expr =
        LBExpr::from_string(expr, ring).map_err(|e| format!("Failed to parse expression: {e}"))?;

    Ok(linear_mba::deobfuscate_lbexpr(expr, cfg, ring)
        .display_function(formatter.to_rust(), "f", ring)
        .to_string())
}
