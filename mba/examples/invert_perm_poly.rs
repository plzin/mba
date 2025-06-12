use itertools::Itertools;
use mba::rings::BinaryRing;
use mba::poly::Poly;
use mba::perm_poly::{
    ZeroIdeal, compute_inverse, compute_inverse_generator,
    compute_inverse_interpolation,
};

enum InverseAlgorithm {
    Newton,
    Fermat,
    Lagrange,
}

fn main() {
    let mut args = std::env::args();
    let Some((_, bits, poly)) = args.next_tuple() else {
        println!("Arguments: <bits> <expr>");
        return;
    };

    let alg = args.next();
    let alg = match alg.as_deref() {
        None | Some("newton") => InverseAlgorithm::Newton,
        Some("fermat") => InverseAlgorithm::Fermat,
        Some("lagrange") => InverseAlgorithm::Lagrange,
        Some(s) => return println!("Unknown algorithm: {s}."),
    };

    let Ok(bits) = bits.parse() else {
        println!("Invalid number of bits: {bits}.");
        return;
    };

    mba::choose_binary_ring!(do_compute_inverse(poly, alg, &r), r = bits);
}

fn do_compute_inverse<R: BinaryRing>(
    poly: String,
    alg: InverseAlgorithm,
    r: &R,
) {
    let p = match Poly::parse(poly, r) {
        Ok(p) => p,
        Err(e) => return println!("Invalid polynomial: {e}"),
    };

    let zi = ZeroIdeal::init(r);

    let q = match alg {
        InverseAlgorithm::Newton => compute_inverse(&p, &zi, r),
        InverseAlgorithm::Fermat => compute_inverse_generator(&p, &zi, r),
        InverseAlgorithm::Lagrange => compute_inverse_interpolation(&p, &zi, r),
    };

    println!("{q}")
}