use mba::{
    choose_binary_ring,
    perm_poly::{ZeroIdeal, compute_inverse, is_perm_poly, random_perm_poly},
    poly::Poly,
    rings::BinaryRing,
};
use wasm_bindgen::prelude::*;

#[wasm_bindgen(js_name = "invertPermutationPolynomial")]
pub fn invert_permutation_polynomial(poly: String, bits: u32) -> Result<String, String> {
    choose_binary_ring!(invert_permutation_polynomial_impl(poly, &r), r = bits)
}

fn invert_permutation_polynomial_impl<R: BinaryRing>(
    poly: String,
    r: &R,
) -> Result<String, String> {
    let zi = ZeroIdeal::init(r);
    let poly = Poly::parse(poly, r)?;

    if !is_perm_poly(&poly) {
        return Err("Polynomial is not a permutation polynomial".into());
    }

    let inv = compute_inverse(&poly, &zi, r);
    Ok(inv.to_tex().to_string())
}

#[wasm_bindgen(js_name = "randomPermutationPolynomial")]
pub fn random_permutation_polynomial(bits: u32) -> String {
    choose_binary_ring!(random_permutation_polynomial_impl(&r), r = bits)
}

fn random_permutation_polynomial_impl<R: BinaryRing>(r: &R) -> String {
    let rng = &mut rand::rng();
    let zi = ZeroIdeal::init(r);
    // This is the smallest degree possible that can represent any permutation
    // that has a polynomial representation.
    let degree = zi.generators().last().unwrap().len() - 1;

    let poly = random_perm_poly(rng, degree, r).simplified(&zi, r);
    poly.display("X").to_string()
}
