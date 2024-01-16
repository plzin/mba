use mba::linear_mba;
use mba::expr::Expr;

fn main() {
    // The example will obfuscate x+y where x and y are 8 bit integers.
    // This should be about the same as the WASM obfuscation.
    // Note that the example obfuscation is quite limited and only meant as
    // a starting point and an example of how to combine the primitives.
    let mut e = Expr::from_string("x+y").unwrap();
    let cfg = linear_mba::ObfuscationConfig::default();
    linear_mba::obfuscate(&mut e, 8, &cfg);
    println!("{e}");
}