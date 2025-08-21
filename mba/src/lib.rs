//! Mixed Boolean-Arithmetic expressions and operations on them.

// See `CustomMetadataSlice` for an explanation of why I use this.
#![allow(internal_features)]
#![feature(ptr_metadata)]
#![feature(core_intrinsics)]

// It would be nicer to import the symbol_table crate ourselves,
// and re-export GlobalSymbol, but we don't want to run into
// version conflicts.
pub use egg::Symbol;

pub mod bitwise_expr;
pub mod expr;
pub mod formatter;
pub mod lattice;
pub mod linear_mba;
pub mod matrix;
pub mod perm_poly;
pub mod poly;
pub mod rings;
pub mod simplify_boolean;
pub mod solver;
pub mod tex;
pub mod valuation;
pub mod vector;

use expr::*;

// Also used for the hacky vector/matrix views.
// See below.
#[cfg(target_pointer_width = "64")]
type Half = u32;
#[cfg(target_pointer_width = "32")]
type Half = u16;

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
/// So I am abusing the usize to store the data.
struct CustomMetadataSlice<T, M> {
    phantom: std::marker::PhantomData<M>,
    data: [T],
}

trait CustomMetadata {
    /// The size of the slice.
    fn size(&self) -> usize;
}

impl<T, M: CustomMetadata> CustomMetadataSlice<T, M> {
    pub unsafe fn new<'a>(data: *const T, metadata: M) -> &'a Self {
        assert_eq!(
            std::mem::size_of::<M>(),
            std::mem::size_of::<<[T] as std::ptr::Pointee>::Metadata>()
        );
        unsafe {
            &*std::ptr::from_raw_parts(data as _, core::intrinsics::transmute_unchecked(metadata))
        }
    }

    pub unsafe fn new_mut<'a>(data: *mut T, metadata: M) -> &'a mut Self {
        assert_eq!(
            std::mem::size_of::<M>(),
            std::mem::size_of::<<[T] as std::ptr::Pointee>::Metadata>()
        );
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(
                data as _,
                core::intrinsics::transmute_unchecked(metadata),
            )
        }
    }

    pub fn metadata(&self) -> M {
        unsafe { core::intrinsics::transmute_unchecked(std::ptr::metadata(self)) }
    }

    pub fn as_ptr(&self) -> *const T {
        self.data.as_ptr()
    }

    pub fn as_mut_ptr(&mut self) -> *mut T {
        self.data.as_mut_ptr()
    }

    pub fn slice(&self) -> &[T] {
        unsafe { std::slice::from_raw_parts(self.data.as_ptr() as _, self.metadata().size()) }
    }

    pub fn slice_mut(&mut self) -> &mut [T] {
        unsafe {
            std::slice::from_raw_parts_mut(self.data.as_mut_ptr() as _, self.metadata().size())
        }
    }
}

/// Hacky function to convert a rug::Integer to a z3 bitvector.
#[cfg(feature = "z3")]
pub(crate) fn int_to_bv<'ctx, R: rings::BinaryRing>(
    ctx: &'ctx z3::Context,
    e: &R::Element,
    r: &R,
) -> z3::ast::BV<'ctx> {
    use z3::ast::Ast;
    let bits: Vec<_> = (0..r.bits()).map(|j| R::bit(e, j)).collect();
    unsafe {
        let ast = z3_sys::Z3_mk_bv_numeral(
            *(ctx as *const _ as *const z3_sys::Z3_context),
            r.bits(),
            bits.as_ptr(),
        );
        z3::ast::BV::wrap(ctx, ast)
    }
}
