use itertools::Itertools;
use mba::poly::Poly;
use mba::perm_poly::{
    ZeroIdeal, compute_inverse, compute_inverse_generator,
    compute_inverse_interpolation,
};

fn main() {
    let mut args = std::env::args();
    let Some((_, bits, poly)) = args.next_tuple() else {
        println!("Arguments: <bits> <expr>");
        return;
    };

    let alg = args.next();
    let alg = match alg.as_ref().map(String::as_str) {
        None | Some("newton") => compute_inverse,
        Some("fermat") => compute_inverse_generator,
        Some("lagrange") => compute_inverse_interpolation,
        _ => return println!("Unknown algorithm."),
    };

    let Ok(bits) = bits.parse() else {
        println!("Invalid number of bits.");
        return;
    };

    let p = match Poly::parse(poly) {
        Ok(p) => p,
        Err(e) => return println!("Invalid polynomial: {e}"),
    };

    let zi = ZeroIdeal::init(bits);
    let q = alg(&p, &zi);
    println!("{q}");
}