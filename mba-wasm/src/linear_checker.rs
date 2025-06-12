use mba::{expr::Expr, linear_mba, rings::BinaryRing};
use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = "probabilisticLinearMBACheck")]
pub fn probabilistic_linear_mba_check(
    expr: String,
    #[wasm_bindgen(js_name = "numInputs")]
    num_inputs: usize,
    bits: u32,
) -> Result<bool, String> {
    mba::choose_binary_ring!(probabilistic_linear_mba_check_impl(
        expr, num_inputs, &r
    ), r = bits)
}

fn probabilistic_linear_mba_check_impl<R: BinaryRing>(
    expr: String,
    num_inputs: usize,
    r: &R,
) -> Result<bool, String> {
    let expr = Expr::from_string(expr, r).map_err(|e| format!("Failed to parse expression: {e}"))?;

    let mut rng = rand::rng();

    Ok(linear_mba::is_probably_linear_mba(&expr, num_inputs, &mut rng, r))
}

#[wasm_bindgen(js_name = "fullLinearMBACheck")]
pub fn full_linear_mba_check(
    expr: String,
    bits: u32,
) -> Result<bool, String> {
    mba::choose_binary_ring!(full_linear_mba_check_impl(expr, &r), r = bits)
}

fn full_linear_mba_check_impl<R: BinaryRing>(
    expr: String,
    r: &R,
) -> Result<bool, String> {
    let expr = Expr::from_string(expr, r).map_err(|e| format!("Failed to parse expression: {e}"))?;

    Ok(linear_mba::is_linear_mba(&expr, r))
}