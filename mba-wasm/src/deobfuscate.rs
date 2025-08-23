use mba::{
    choose_binary_ring,
    expr::Expr,
    linear_mba::{self, DeobfuscationConfig},
    rings::BinaryRing,
};
use wasm_bindgen::prelude::*;

use crate::{Formatter, SolutionAlgorithm};

#[wasm_bindgen(js_name = "deobfuscate")]
pub fn deobfuscate(
    expr: String,
    bits: u32,
    alg: SolutionAlgorithm,
    #[wasm_bindgen(js_name = "detectBoolean")] detect_boolean: bool,
    formatter: Formatter,
) -> Result<String, String> {
    let cfg =
        DeobfuscationConfig { alg: alg.to_rust(), boolean: detect_boolean };
    choose_binary_ring!(deobfuscate_impl(expr, &cfg, formatter, &r), r = bits)
}

fn deobfuscate_impl<R: BinaryRing>(
    expr: String,
    cfg: &DeobfuscationConfig,
    formatter: Formatter,
    ring: &R,
) -> Result<String, String> {
    // Parse the expression.
    let mut expr = Expr::from_string(expr, ring)
        .map_err(|e| format!("Failed to parse expression: {e}"))?;

    linear_mba::deobfuscate(&mut expr, cfg, ring);
    Ok(expr.display_function(formatter.to_rust(), "f", ring).to_string())
}
