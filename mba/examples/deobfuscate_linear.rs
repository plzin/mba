use mba::bitwise_expr::LBExpr;
use mba::linear_mba::{deobfuscate_lbexpr, DeobfuscationConfig, SolutionAlgorithm};
use mba::rings::U8;

fn main() {
    //let e = LBExpr::from_string("(x ^ y) + 2 * (x & y)".to_owned(), &U8).unwrap();
    let e = LBExpr::from_string("\
        145 + 81 * (aux1 & y) + 198 * ~~(aux1 ^ x) + 105 * (~(y & y) | ~aux0 ^ \
        (aux1 | y)) + 242 * (~aux1 | y | aux1 | ~x) + 189 * aux0 + 104 * (aux1 \
        & aux0) + 156 * y + 71 * ~aux0 + 27 * (aux1 & aux1 & aux1 | aux1 & \
        ~aux0) + 154 * (aux0 ^ ~y ^ ~x) + 249 * x + 37 * ~~(aux0 | aux1) + 67 \
        * aux1 + 116 * (aux1 | x) + 186 * ~((aux1 | aux1) & y & y) + 172 * \
        (aux0 ^ (aux1 | aux1 | aux0 | aux0)) + 151 * ((aux0 | aux1 ^ y) ^ y ^ \
        aux0 ^ aux1) + 152 * (~y ^ aux0 & x & (y ^ x)) + 105 * (~(aux1 | y) ^ \
        aux0 ^ aux1 ^ y ^ aux1) + 104 * ~~(aux0 | x) + 126 * (aux0 & y & (aux0 \
        ^ y) & ~aux1) + 52 * ~(aux0 ^ x | x ^ y)".to_owned(), &U8).unwrap();

    let cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::ShortVector,
        boolean: false,
    };
    let d = deobfuscate_lbexpr(
        e, &cfg, &U8
    );
    println!("{d}");
}
