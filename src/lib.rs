//! Mixed Boolean-Arithmetic expressions and operations on them.

// See [CustomMetadataSlice] for an explanation of why I use this.
#![feature(ptr_metadata)]

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

/// I do some hacky stuff with vector and matrix views.
/// The other option would be having two structs for immutable and mutable
/// views, but that is just annoying and rust also doesn't have two structs for
/// immutable and mutable slices. Whether you can mutate the data is determined
/// by whether you have a mutable reference to it.
/// The way this is implemented in rust is that a pointer (and thus reference)
/// has metadata attached to it. In the case of a slice, this is the length.
/// But in the case of more complicated view types, this data could be much
/// more complicated, like a stride and a size, or in the case of a matrix
/// just the number of rows and columns. Unfortunately, rust doesn't let you
/// choose the metadata type and in our case it is a usize
/// (see https://doc.rust-lang.org/stable/std/ptr/trait.Pointee.html),
/// which is not enough space to store the data we need.
/// So I am abusing the usize to store a pointer to the data.
struct CustomMetadataSlice<T, M> {
    phantom: std::marker::PhantomData<M>,
    data: [T]
}

trait CustomMetadata {
    /// The size of the slice.
    fn size(&self) -> usize;
}

impl<T, M: CustomMetadata> CustomMetadataSlice<T, M> {
    pub fn new<'a>(data: *const T, metadata: M) -> &'a Self {
        let metadata = Box::into_raw(Box::new(metadata));
        unsafe {
            &*std::ptr::from_raw_parts(
                data as _,
                metadata as usize
            )
        }
    }

    pub fn new_mut<'a>(data: *mut T, metadata: M) -> &'a mut Self {
        let metadata = Box::into_raw(Box::new(metadata));
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(
                data as _,
                metadata as usize
            )
        }
    }

    pub fn metadata(&self) -> &M {
        unsafe {
            &*(std::ptr::metadata(self) as *const M)
        }
    }

    pub fn as_ptr(&self) -> *const T {
        self.data.as_ptr()
    }

    pub fn as_mut_ptr(&mut self) -> *mut T {
        self.data.as_mut_ptr()
    }

    pub fn slice(&self) -> &[T] {
        unsafe {
            std::slice::from_raw_parts(
                self.data.as_ptr() as _,
                self.metadata().size()
            )
        }
    }

    pub fn slice_mut(&mut self) -> &mut [T] {
        unsafe {
            std::slice::from_raw_parts_mut(
                self.data.as_mut_ptr() as _,
                self.metadata().size()
            )
        }
    }
}

impl<T, M> Drop for CustomMetadataSlice<T, M> {
    fn drop(&mut self) {
        unsafe {
            let _ = Box::from_raw(
                self.data.as_ptr() as *mut M
            );
        }
    }
}

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