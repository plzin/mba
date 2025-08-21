#![feature(iterator_try_collect)]

mod deobfuscate;
mod deobfuscate_linear;
mod linear_checker;
mod linear_system;
mod obfuscate;
mod obfuscate_linear;
mod perm_poly;

use mba::bitwise_expr::LBExpr;
use mba::rings::Z;
use wasm_bindgen::prelude::*;

/// Sets the panic hook to display useful error messages.
#[wasm_bindgen(js_name = "setPanicHook")]
pub fn set_panic_hook() {
    std::panic::set_hook(Box::new(console_error_panic_hook::hook));
}

/// Parses a linear MBA expression and prints it.
///
/// This is used by the linear obfuscation and deobfuscation page to sanitize
/// and normalize the rewrite operations.
///
/// We parse the expressions into integer ([`mba::rings::Z`]) expressions so
/// that arbitrarily large coefficients are supported. They will be parsed by
/// the particular ring into their correct elements, so if a constant 1025
/// appears, it will be 1 in [`mba::rings::U8`].
#[wasm_bindgen(js_name = "parseAndPrintLBExpr")]
pub fn parse_and_print_lbexpr(expr: String) -> Result<String, String> {
    LBExpr::from_string(expr, &Z).map(|expr| expr.to_string())
}

#[wasm_bindgen]
pub enum Formatter {
    C,
    Rust,
    Tex,
    LLVM,
}

impl Formatter {
    pub fn to_rust(&self) -> mba::formatter::Formatter {
        match self {
            Formatter::C => mba::formatter::Formatter::C,
            Formatter::Rust => mba::formatter::Formatter::Rust,
            Formatter::Tex => mba::formatter::Formatter::Tex,
            Formatter::LLVM => mba::formatter::Formatter::LLVM,
        }
    }
}

#[wasm_bindgen]
pub enum SolutionAlgorithm {
    Fast,
    LeastComplexTerms,
    ShortVector,
}

impl SolutionAlgorithm {
    fn to_rust(&self) -> mba::linear_mba::SolutionAlgorithm {
        match self {
            SolutionAlgorithm::Fast => mba::linear_mba::SolutionAlgorithm::Fast,
            SolutionAlgorithm::LeastComplexTerms => {
                mba::linear_mba::SolutionAlgorithm::LeastComplexTerms
            }
            SolutionAlgorithm::ShortVector => mba::linear_mba::SolutionAlgorithm::ShortVector,
        }
    }
}
