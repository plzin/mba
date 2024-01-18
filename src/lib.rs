//! Mixed Boolean-Arithmetic expressions and operations on them.

#![feature(ptr_metadata)]

// You can only compile this for 64-bits because I use some hacky stuff
// with vector views. This can be removed once rust allows supports
// custom dynamically sized types.
#[cfg(not(target_pointer_width = "64"))]
compile_error!("This crate only works on 64-bit systems.");

// It would be nicer to import the symbol_table crate ourselves,
// and re-export GlobalSymbol, but we don't want to run into
// version conflicts.
pub use egg::Symbol;

pub mod matrix;

pub mod vector;

pub mod valuation;

pub mod expr;

pub mod uniform_expr;


pub mod diophantine;
pub mod lattice;
pub mod poly;
pub mod perm_poly;
pub mod linear_mba;
pub mod simplify_boolean;

use std::ops::Neg;
use num_bigint::BigInt;
use num_traits::{One, Euclid};
use expr::*;

// The following functions were supported by rug but not num_bigint.
// They are way more inefficient than the rug versions, so they should
// ideally be implemented in num_bigint directly.
// If you see code that uses them and it could easily be improved,
// then that is why.
pub(crate) fn keep_bits(i: &BigInt, n: u32) -> BigInt {
    let m = BigInt::one() << n;
    i.rem_euclid(&m)
}

pub(crate) fn keep_bits_mut(i: &mut BigInt, n: u32) {
    *i = keep_bits(i, n);
}

pub(crate) fn keep_signed_bits(i: &BigInt, n: u32) -> BigInt {
    let m = BigInt::one() << n;
    let mut r = i.rem_euclid(&m);
    if i.bit(n as u64 - 1) {
        r -= m;
    }

    r
}

pub(crate) fn keep_signed_bits_mut(i: &mut BigInt, n: u32) {
    *i = keep_signed_bits(i, n);
}

/// The `rug` crate had a `NegAssign` trait for this.
/// This is my workaround for now.
/// We could add a `NegAssign` trait here, but it would
/// make way more sense in the [num_traits] crate.
pub(crate) fn neg_mut<T: Neg<Output = T>>(i: &mut T) {
    unsafe {
        let temp = std::ptr::read(i);
        let temp = -temp;
        std::ptr::write(i, temp);
    }
}

/// Hacky function to convert a rug::Integer to a z3 bitvector.
#[cfg(feature = "z3")]
pub(crate) fn int_to_bv<'ctx>(
    ctx: &'ctx z3::Context, width: u32, i: &BigInt
) -> z3::ast::BV<'ctx> {
    use z3::ast::Ast;
    let bits: Vec<_> = (0..width as u64).map(|j| i.bit(j)).collect();
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
) -> Option<BigInt> {
    let c = it.next()?;
    let mut num: BigInt = c.to_digit(10)?.into();

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