//! Mixed Boolean-Arithmetic expressions and operations on them.

#![feature(ptr_metadata)]

// You can only compile this for 64-bits because I use some hacky stuff
// with vector views. This can be removed once rust allows supports
// custom dynamically sized types.
#[cfg(not(target_pointer_width = "64"))]
compile_error!("This crate only works on 64-bit systems.");

use rug::Integer;

#[cfg(feature = "z3")]
use z3::{self, ast::Ast};

// It would be nicer to import the symbol_table crate ourselves,
// and re-export GlobalSymbol, but we don't want to run into
// version conflicts.
pub use egg::Symbol;

pub mod matrix;

pub mod vector;

pub mod valuation;

pub mod expr;
use expr::*;

pub mod uniform_expr;


pub mod diophantine;
pub mod lattice;
pub mod poly;
pub mod perm_poly;
pub mod linear_mba;
pub mod simplify_boolean;

/// Hacky function to convert a rug::Integer to a z3 bitvector.
#[cfg(feature = "z3")]
pub(crate) fn int_to_bv<'ctx>(
    ctx: &'ctx z3::Context, width: u32, i: &Integer
) -> z3::ast::BV<'ctx> {
    let mut bits = i.to_digits::<bool>(rug::integer::Order::Lsf);
    bits.resize(width as usize, false);
    unsafe {
        let ast = z3_sys::Z3_mk_bv_numeral(
            *(ctx as *const _ as *const z3_sys::Z3_context),
            width,
            bits.as_ptr()
        );
        z3::ast::BV::wrap(ctx, ast)
    }
}

/// Used internally to convert an integer from and iterator over characters.
pub(crate) fn int_from_it(
    it: &mut std::iter::Peekable<std::str::Chars>
) -> Option<Integer> {
    let c = it.next()?;
    let mut num: Integer = c.to_digit(10)?.into();

    let mut lookahead = it.peek();
    while let Some(c) = lookahead {
        if !c.is_ascii_digit() {
            break;
        }

        num *= 10;
        num += c.to_digit(10).unwrap();
        it.next();
        lookahead = it.peek();
    }

    Some(num)
}