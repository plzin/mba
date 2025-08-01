use mba::formatter::Formatter;
use mba::linear_mba::{obfuscate, ObfuscationConfig};
use mba::expr::Expr;
use mba::rings::U8;
use rand::{SeedableRng as _, rngs::StdRng};

fn main() {
    let mut rng = StdRng::seed_from_u64(0);
    // The example will obfuscate x+y where x and y are 8 bit integers.
    // This should be about the same as the WASM obfuscation.
    // Note that the example obfuscation is quite limited and only meant as
    // a starting point and an example of how to combine the primitives.
    let r = &U8;
    let mut e = Expr::from_string("x+y".to_owned(), r).unwrap();
    let cfg = ObfuscationConfig::default();
    obfuscate(&mut e, &cfg, &mut rng, r);
    println!("{}", e.display_function(Formatter::C, "f", r));
}