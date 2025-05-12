use mba::uniform_expr::LUExpr;
use mba::linear_mba::{self, *};

fn main() {
    let e = LUExpr::from_string("(x ^ y) + 2 * (x & y)").unwrap();
    let cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::LeastComplexTerms,
        boolean: false,
    };
    let d = linear_mba::deobfuscate_luexpr(
        e, 8, &cfg
    );
    println!("{d}");
}
