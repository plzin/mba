//! Owned and non-owned vectors constant runtime dimension.

use std::borrow::{Borrow, BorrowMut};
use std::marker::PhantomData;
use std::ops::{
    Deref, DerefMut, Neg, Add, Sub, Mul, Div,
    AddAssign, SubAssign, MulAssign, DivAssign, Index, Range, IndexMut
};
use std::fmt::Debug;
use num_traits::Zero;
use rug::{Integer, Complete, Float, Rational, ops::NegAssign};
use crate::matrix::{RowVector, ColumnVector};
use crate::select;

/// How are the entries of a vector stored?
pub trait VectorStorage<T> {
    type Iter<'a>: Iterator<Item = &'a T> where Self: 'a, T: 'a;
    type IterMut<'a>: Iterator<Item = &'a mut T> where Self: 'a, T: 'a;

    /// The dimension of the vector.
    fn dim(&self) -> usize;

    /// Reference to an entry at an index.
    fn entry(&self, idx: usize) -> &T;

    /// Mutable reference to an entry at an index.
    fn entry_mut(&mut self, idx: usize) -> &mut T;

    /// Returns an iterator over the elements.
    fn iter(&self) -> Self::Iter<'_>;

    /// Returns an iterator over the mutable elements.
    fn iter_mut(&mut self) -> Self::IterMut<'_>;
}

pub trait ContiguousVectorStorage<T>: VectorStorage<T> {
    /// Returns a slice of the vector.
    fn as_slice(&self) -> &[T];

    /// Returns a mutable slice of the vector.
    fn as_slice_mut(&mut self) -> &mut [T];
}

pub struct Vector<T, S: VectorStorage<T> + ?Sized> {
    phantom: PhantomData<T>,
    storage: S,
}

impl<T, S: VectorStorage<T> + ?Sized> Vector<T, S> {
    /// Returns the dimension of the vector.
    pub fn dim(&self) -> usize {
        self.storage.dim()
    }

    /// Returns a reference to an entry at an index.
    pub fn entry(&self, idx: usize) -> &T {
        self.storage.entry(idx)
    }

    /// Returns a mutable reference to an entry at an index.
    pub fn entry_mut(&mut self, idx: usize) -> &mut T {
        self.storage.entry_mut(idx)
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> S::Iter<'_> {
        self.storage.iter()
    }

    /// Returns an iterator over the mutable elements.
    pub fn iter_mut(&mut self) -> S::IterMut<'_> {
        self.storage.iter_mut()
    }

    /// Is the vector empty, i.e. dimension zero?
    pub fn is_empty(&self) -> bool {
        self.dim() == 0
    }

    /// Apply a function to each entry.
    pub fn transform<U, F: FnMut(&T) -> U>(&self, f: F) -> OwnedVector<U> {
        Vector::from_iter(self.dim(), self.iter().map(f))
    }

    /// Apply a function to each entry.
    pub fn map_mut(&mut self, mut f: impl FnMut(&mut T)) {
        for e in self.iter_mut() {
            f(e);
        }
    }

    /// Swap the two entries at index i and j.
    pub fn swap(&mut self, i : usize, j: usize) {
        unsafe {
            std::ptr::swap(
                self.entry_mut(i),
                self.entry_mut(j)
            );
        }
    }

    /// Convert the entries of this vector to a different type.
    /// This function is a bit restrictive because
    /// (e.g.) u32 does not implement From<&u32>.
    /// So what we do is clone the entries in the vector
    /// and then call into() on them.
    /// So because this clone is almost always unnecessary,
    /// we only implement this for Copy types.
    pub fn to<U>(&self) -> OwnedVector<U> where T: Copy + Into<U> {
        #[allow(clippy::clone_on_copy)]
        self.transform(|e| e.clone().into())
    }

    /// Converts the entries of this vector to rug::Floats.
    pub fn to_float(&self, prec: u32) -> OwnedVector<Float>
    where
        for<'a> Float: rug::Assign<&'a T>
    {
        self.transform(|e| Float::with_val(prec, e))
    }
}

impl<T, S: ContiguousVectorStorage<T> + ?Sized> Vector<T, S> {
    /// Returns a slice of the vector.
    pub fn as_slice(&self) -> &[T] {
        self.storage.as_slice()
    }

    /// Returns a mutable slice of the vector.
    pub fn as_slice_mut(&mut self) -> &mut [T] {
        self.storage.as_slice_mut()
    }

    /// Converts the vector into a row vector view, i.e.
    /// a matrix with one row.
    pub fn row_vector(&self) -> &RowVector<T> {
        unsafe { &*(self.as_slice() as *const _ as *const _) }
    }

    /// Converts the vector into a mutable row vector view, i.e.
    /// a matrix with one row.
    pub fn row_vector_mut(&mut self) -> &mut RowVector<T> {
        unsafe { &mut *(self.as_slice_mut() as *mut _ as *mut _) }
    }

    /// Converts the vector into a column vector view, i.e.
    /// a matrix with one column.
    pub fn column_vector(&self) -> &ColumnVector<T> {
        unsafe { &*(self.as_slice() as *const _ as *const _) }
    }

    /// Converts the vector into a mutable column vector view, i.e.
    /// a matrix with one column.
    pub fn column_vector_mut(&mut self) -> &mut ColumnVector<T> {
        unsafe { &mut *(self.as_slice_mut() as *mut _ as *mut _) }
    }
}

impl<T, S: ContiguousVectorStorage<T>> AsRef<[T]> for Vector<T, S> {
    fn as_ref(&self) -> &[T] {
        self.as_slice()
    }
}

impl<T, S: ContiguousVectorStorage<T>> AsMut<[T]> for Vector<T, S> {
    fn as_mut(&mut self) -> &mut [T] {
        self.as_slice_mut()
    }
}

impl<T, S: VectorStorage<T>> Vector<T, S> {
    /// Construct the vector from its storage.
    pub fn from_storage(storage: S) -> Self {
        Self {
            phantom: std::marker::PhantomData,
            storage,
        }
    }
}

impl<S: VectorStorage<Float> + ?Sized> Vector<Float, S> {
    pub fn precision(&self) -> u32 {
        assert!(self.dim() > 0,
            "Cannot get precision of an empty vector");
        if cfg!(debug_assert) {
            self.assert_precision();
        }
        self.entry(0).prec()
    }

    pub fn assert_precision(&self) {
        let mut iter = self.iter();
        let Some(prec) = iter.next().map(|f| f.prec()) else {
            return
        };
        assert!(iter.all(|f| f.prec() == prec));
    }
}

impl<T, S: VectorStorage<T> + ?Sized> Index<usize> for Vector<T, S> {
    type Output = T;

    fn index(&self, idx: usize) -> &Self::Output {
        self.entry(idx)
    }
}

impl<T, S: VectorStorage<T> + ?Sized> IndexMut<usize> for Vector<T, S> {
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        self.entry_mut(idx)
    }
}

impl<T, S: ContiguousVectorStorage<T> + ?Sized> Index<Range<usize>> for Vector<T, S> {
    type Output = VectorView<T>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        self.as_slice().index(index).into()
    }
}

impl<T, S: ContiguousVectorStorage<T> + ?Sized> IndexMut<Range<usize>> for Vector<T, S> {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        self.as_slice_mut().index_mut(index).into()
    }
}

impl<'a, T, S: VectorStorage<T> + ?Sized> IntoIterator for &'a Vector<T, S> {
    type Item = &'a T;
    type IntoIter = S::Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, T, S: VectorStorage<T> + ?Sized> IntoIterator for &'a mut Vector<T, S> {
    type Item = &'a mut T;
    type IntoIter = S::IterMut<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

impl<T, S, I: ?Sized, U> PartialEq<I> for Vector<T, S> where
    T: PartialEq<U>,
    S: VectorStorage<T> + ?Sized,
    for<'a> &'a I: IntoIterator<Item = &'a U>,
{
    fn eq(&self, other: &I) -> bool {
        self.iter().eq(other)
    }
}

impl<T: std::fmt::Debug, S: VectorStorage<T> + ?Sized> std::fmt::Debug for Vector<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<T: Eq, S: VectorStorage<T> + ?Sized> Eq for Vector<T, S> {}

/// Technically this is not a vector view, because it owns
/// the elements. However, you can only get references
/// to this type, which are the views.
pub type VectorView<T> = Vector<T, SliceVectorStorage<T>>;
pub type FVectorView = VectorView<Float>;
pub type IVectorView = VectorView<Integer>;
pub type RVectorView = VectorView<Rational>;

impl<T> VectorView<T> {
    /// Get a vector view from a slice.
    pub fn from_slice(s: &[T]) -> &Self {
        unsafe {
            &*std::ptr::from_raw_parts(s.as_ptr() as _, s.len())
        }
    }

    /// Get a mutable vector view from a slice.
    pub fn from_slice_mut(s: &mut [T]) -> &mut Self {
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(s.as_mut_ptr() as _, s.len())
        }
    }
}

impl<'a, T> From<&'a [T]> for &'a VectorView<T> {
    fn from(s: &'a [T]) -> Self {
        VectorView::from_slice(s)
    }
}

impl<'a, T> From<&'a mut [T]> for &'a mut VectorView<T> {
    fn from(s: &'a mut [T]) -> Self {
        VectorView::from_slice_mut(s)
    }
}

impl<T: Clone> ToOwned for VectorView<T> {
    type Owned = OwnedVector<T>;

    fn to_owned(&self) -> Self::Owned {
        OwnedVector::from_iter(self.dim(), self.iter().cloned())
    }
}

#[derive(Debug)]
#[repr(transparent)]
pub struct SliceVectorStorage<T>([T]);

impl<T> VectorStorage<T> for SliceVectorStorage<T> {
    type Iter<'a> = std::slice::Iter<'a, T> where T: 'a;
    type IterMut<'a> = std::slice::IterMut<'a, T> where T: 'a;

    fn dim(&self) -> usize {
        self.0.len()
    }

    fn entry(&self, idx: usize) -> &T {
        &self.0[idx]
    }

    fn entry_mut(&mut self, idx: usize) -> &mut T {
        &mut self.0[idx]
    }

    fn iter(&self) -> Self::Iter<'_> {
        self.0.iter()
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        self.0.iter_mut()
    }
}

impl<T> ContiguousVectorStorage<T> for SliceVectorStorage<T> {
    fn as_slice(&self) -> &[T] {
        &self.0
    }

    fn as_slice_mut(&mut self) -> &mut [T] {
        &mut self.0
    }
}

pub type OwnedVector<T> = Vector<T, OwnedVectorStorage<T>>;
pub type FOwnedVector = OwnedVector<Float>;
pub type IOwnedVector = OwnedVector<Integer>;
pub type ROwnedVector = OwnedVector<Rational>;

impl<T> OwnedVector<T> {
    /// Returns a view of the owned vector.
    /// Currently, this has the same layout as the owned vector,
    /// so this is pretty much a no-op.
    pub fn view(&self) -> &VectorView<T> {
        VectorView::from_slice(self.as_slice())
    }

    /// Returns a mutable view of the owned vector.
    /// Currently, this has the same layout as the owned vector,
    /// so this is pretty much a no-op.
    pub fn view_mut(&mut self) -> &mut VectorView<T> {
        VectorView::from_slice_mut(self.as_slice_mut())
    }

    /// Returns an empty vector.
    pub fn empty() -> Self {
        Self::from_raw_parts(core::ptr::null_mut(), 0)
    }

    /// Creates a vector from a slice.
    pub fn from_entries<U, V>(a: U) -> Self
    where
        U: AsRef<[V]>,
        V: Into<T> + Clone,
    {
        let a = a.as_ref();
        Self::from_iter(a.len(), a.iter().cloned().map(Into::into))
    }

    /// Creates a vector from an iterator.
    pub fn from_iter<U: Iterator<Item = T>>(
        dim: usize, mut iter: U
    ) -> Self {
        let v = Self::uninit(dim);
        for i in 0..dim {
            unsafe {
                let e = iter.next()
                    .expect("Iter needs to return `dim` elements.");
                v.storage.entries.add(i).write(e);
            }
        }

        v
    }

    /// Creates a vector from an array.
    pub fn from_array<U: Into<T>, const D: usize>(a: [U; D]) -> Self {
        Self::from_iter(D, a.into_iter().map(U::into))
    }

    /// Applies a function to each entry.
    pub fn map<F: FnMut(T) -> T>(self, mut f: F) -> OwnedVector<T> {
        for i in 0..self.dim() {
            unsafe {
                let e = self.storage.entries.add(i).read();
                self.storage.entries.add(i).write(f(e));
            }
        }

        self
    }

    /// Appends an element to the end.
    pub fn append(&mut self, e: T) {
        self.storage.append(e)
    }

    /// Returns an owned vector from a pointer and dimension.
    pub(self) fn from_raw_parts(entries: *mut T, dim: usize) -> Self {
        Self::from_storage(OwnedVectorStorage { entries, dim })
    }

    /// Returns an uninitialized vector.
    pub(self) fn uninit(dim: usize) -> Self {
        if dim == 0 {
            return Self::empty();
        }

        let layout = std::alloc::Layout::from_size_align(
            dim * core::mem::size_of::<T>(),
            core::mem::align_of::<T>()
        ).unwrap();

        let entries = unsafe {
            std::alloc::alloc(layout) as *mut T
        };

        Self::from_raw_parts(entries, dim)
    }
}

impl OwnedVector<Float> {
    pub fn zero_prec(dim: usize, prec: u32) -> Self {
        Self::from_iter(dim, core::iter::repeat(Float::with_val(prec, 0)))
    }
}

impl<T: Zero> OwnedVector<T> {
    /// Returns a zero vector.
    pub fn zero(dim: usize) -> Self {
        Self::from_iter(dim, core::iter::repeat_with(|| T::zero()))
    }
}

impl<T> Borrow<VectorView<T>> for OwnedVector<T> {
    fn borrow(&self) -> &VectorView<T> {
        self.view()
    }
}

impl<T> BorrowMut<VectorView<T>> for OwnedVector<T> {
    fn borrow_mut(&mut self) -> &mut VectorView<T> {
        self.view_mut()
    }
}

impl<T: Clone> Clone for OwnedVector<T> {
    fn clone(&self) -> Self {
        Self::from_iter(self.dim(), self.iter().cloned())
    }
}

pub struct OwnedVectorStorage<T> {
    /// Memory that holds the entries.
    pub entries: *mut T,

    /// The dimension (number of entries) of the vector.
    pub dim: usize,
}

impl<T> OwnedVectorStorage<T> {
    /// Appends an element to the end.
    pub fn append(&mut self, e: T) {
        let layout = std::alloc::Layout::from_size_align(
            self.dim * core::mem::size_of::<T>(),
            core::mem::align_of::<T>()
        ).unwrap();
        self.dim += 1;
        unsafe {
            // Reallocate the entries.
            self.entries = std::alloc::realloc(
                self.entries as _,
                layout,
                core::mem::size_of::<T>() * self.dim
            ) as _;
            // Write the new element into the last entry.
            self.entries.add(self.dim - 1).write(e);
        };
    }
}

impl<T> VectorStorage<T> for OwnedVectorStorage<T> {
    type Iter<'a> = std::slice::Iter<'a, T> where T: 'a;
    type IterMut<'a> = std::slice::IterMut<'a, T> where T: 'a;

    fn dim(&self) -> usize {
        self.dim
    }

    fn entry(&self, idx: usize) -> &T {
        assert!(idx < self.dim);
        unsafe {
            &*self.entries.add(idx)
        }
    }

    fn entry_mut(&mut self, idx: usize) -> &mut T {
        assert!(idx < self.dim);
        unsafe {
            &mut *self.entries.add(idx)
        }
    }

    fn iter(&self) -> Self::Iter<'_> {
        self.as_slice().iter()
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        self.as_slice_mut().iter_mut()
    }
}

impl<T> ContiguousVectorStorage<T> for OwnedVectorStorage<T> {
    fn as_slice(&self) -> &[T] {
        unsafe {
            std::slice::from_raw_parts(self.entries, self.dim)
        }
    }

    fn as_slice_mut(&mut self) -> &mut [T] {
        unsafe {
            std::slice::from_raw_parts_mut(self.entries, self.dim)
        }
    }
}

impl<T> Drop for OwnedVectorStorage<T> {
    fn drop(&mut self) {
        if self.entries.is_null() {
            return;
        }

        unsafe {
            // Call the destructor on each element.
            for i in 0..self.dim {
                core::ptr::drop_in_place(self.entries.add(i));
            }

            let layout = std::alloc::Layout::from_size_align(
                self.dim * core::mem::size_of::<T>(),
                core::mem::align_of::<T>()
            ).unwrap();

            // Free the memory.
            std::alloc::dealloc(self.entries as _, layout);
        }
    }
}

/// A vector view where elements are skipped.
/// The stride is the amount of elements to skip
/// not the number of bytes to skip.
pub type StrideVectorView<T> = Vector<T, StrideStorage<T>>;
pub type FStrideVectorView = StrideVectorView<Float>;
pub type IStrideVectorView = StrideVectorView<Integer>;
pub type RStrideVectorView = StrideVectorView<Rational>;

impl<T> StrideVectorView<T> {
    /// Returns a vector view with stride from a raw pointer and dimension.
    pub fn from_raw_parts<'a>(ptr: *const T, dim: usize, stride: usize) -> &'a Self {
        unsafe {
            &*std::ptr::from_raw_parts(
                ptr as _,
                Self::metadata(dim, stride)
            )
        }
    }

    /// Returns a mutable vector view with stride from a raw pointer and dimension.
    pub fn from_raw_parts_mut<'a>(ptr: *mut T, dim: usize, stride: usize) -> &'a mut Self {
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(
                ptr as _,
                Self::metadata(dim, stride)
            )
        }
    }

    /// Returns a vector view with stride from a slice.
    /// The vector view will have the maximum possible dimension
    /// that fits into the slice.
    pub fn from_stride_slice(s: &[T], stride: usize) -> &Self {
        assert!(stride > 0);
        let dim = (s.len() + (stride - 1)) / stride;
        unsafe {
            &*std::ptr::from_raw_parts(
                s.as_ptr() as _,
                Self::metadata(dim, stride)
            )
        }
    }

    /// Returns a vector view with stride from a slice.
    /// The vector view will have the maximum possible dimension
    /// that fits into the slice.
    pub fn from_stride_slice_mut(s: &mut [T], stride: usize) -> &mut Self {
        assert!(stride > 0);
        let dim = s.len() + (stride - 1) / stride;
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(
                s.as_mut_ptr() as _,
                Self::metadata(dim, stride)
            )
        }
    }

    /// Returns the pointer metadata for the given dimension and stride.
    fn metadata(dim: usize, stride: usize) -> usize {
        assert!(dim <= u32::MAX as usize);
        assert!(stride <= u32::MAX as usize);
        dim | stride << 32
    }

}

/// This is extremely hacky, because rust currently does not support
/// custom dynamically sized types (DST). We make this a DST by having
/// a slice as the member, just as for `ContiguousStorage`.
/// We abuse the metadata of the slice to store the dimension and stride!
/// (See https://doc.rust-lang.org/std/ptr/trait.Pointee.html, for what
/// metadata is). Since we cannot implement `Pointee` ourselves, we
/// have to use the usize metadata of the slice.
/// The first 32 bits are the dimension, the second 32 bits are the stride.
/// As I said, very hacky, but in practice the both the dimension and
/// stride should fit into 32 bits.
#[repr(transparent)]
pub struct StrideStorage<T>([T]);

impl<T> StrideStorage<T> {
    /// The stride of the vector.
    pub fn stride(&self) -> usize {
        unsafe {
            self.0.len() >> 32
        }
    }
}

impl<T> VectorStorage<T> for StrideStorage<T> {
    type Iter<'a> = StrideStorageIter<'a, T> where T: 'a;
    type IterMut<'a> = StrideStorageIterMut<'a, T> where T: 'a;

    fn dim(&self) -> usize {
        self.0.len() & 0xffff_ffff
    }

    fn entry(&self, idx: usize) -> &T {
        assert!(idx < self.dim());
        unsafe {
            self.0.get_unchecked(idx * self.stride())
        }
    }

    fn entry_mut(&mut self, idx: usize) -> &mut T {
        assert!(idx < self.dim());
        unsafe {
            self.0.get_unchecked_mut(idx * self.stride())
        }
    }

    fn iter(&self) -> Self::Iter<'_> {
        StrideStorageIter {
            ptr: self.0.as_ptr(),
            end: unsafe { self.0.as_ptr().add(self.dim() * self.stride()) },
            stride: self.stride(),
            marker: PhantomData,
        }
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        StrideStorageIterMut {
            ptr: self.0.as_mut_ptr(),
            end: unsafe { self.0.as_mut_ptr().add(self.dim() * self.stride()) },
            stride: self.stride(),
            marker: PhantomData,
        }
    }
}

pub struct StrideStorageIter<'a, T> {
    pub(self) ptr: *const T,
    pub(self) end: *const T,
    pub(self) stride: usize,
    pub(self) marker: PhantomData<&'a T>,
}

impl<'a, T> Iterator for StrideStorageIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr >= self.end {
            None
        } else {
            let ptr = self.ptr;
            self.ptr = unsafe { ptr.add(self.stride) };
            Some(unsafe { &*ptr })
        }
    }
}

pub struct StrideStorageIterMut<'a, T> {
    pub(self) ptr: *mut T,
    pub(self) end: *mut T,
    pub(self) stride: usize,
    pub(self) marker: PhantomData<&'a mut T>,
}

impl<'a, T> Iterator for StrideStorageIterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr >= self.end {
            None
        } else {
            let ptr = self.ptr;
            self.ptr = unsafe { ptr.add(self.stride) };
            Some(unsafe { &mut *ptr })
        }
    }
}

pub trait InnerProduct: Sized {
    /// Computes the inner product of two vectors.
    fn dot<R, S>(v: &Vector<Self, R>, w: &Vector<Self, S>) -> Self
    where
        R: VectorStorage<Self> + ?Sized,
        S: VectorStorage<Self> + ?Sized;

    /// Computes the inner product of a vector with itself.
    fn norm_sqr<S>(v: &Vector<Self, S>) -> Self
        where S: VectorStorage<Self> + ?Sized,
    {
        Self::dot(v, v)
    }
}

pub trait VectorNorm: Sized {
    type Scalar;

    /// Computes the norm of the vector.
    fn norm<S>(v: &Vector<Self, S>) -> Self::Scalar
        where S: VectorStorage<Self> + ?Sized;
}

impl<S, T> Vector<T, S>
where
    S: VectorStorage<T> + ?Sized,
    T: InnerProduct,
{
    /// Computes the inner product of two vectors.
    pub fn dot<R>(&self, v: &Vector<T, R>) -> T
        where R: VectorStorage<T> + ?Sized,
    {
        T::dot(self, v)
    }

    /// Computes the inner product of a vector with itself.
    pub fn norm_sqr(&self) -> T {
        T::norm_sqr(self)
    }
}

impl<S, T> Vector<T, S>
where
    S: VectorStorage<T> + ?Sized,
    T: VectorNorm,
{
    /// Returns the norm of the vector.
    pub fn norm(&self) -> T::Scalar {
        T::norm(self)
    }
}

impl<S, T> Vector<T, S>
where
    S: VectorStorage<T> + ?Sized,
    T: VectorNorm<Scalar = T>,
    for<'a> Vector<T, S>: DivAssign<&'a T>,
{
    pub fn normalize(&mut self) {
        *self /= &T::norm(self);
    }
}

impl InnerProduct for Integer {
    fn dot<R, S>(v: &Vector<Self, R>, w: &Vector<Self, S>) -> Self
    where
        R: VectorStorage<Self> + ?Sized,
        S: VectorStorage<Self> + ?Sized
    {
        assert!(v.dim() == w.dim());
        v.iter().zip(w.iter())
            .map(|(c, d)| c * d)
            .sum()
    }

    fn norm_sqr<S>(v: &Vector<Self, S>) -> Self
        where S: VectorStorage<Self> + ?Sized,
    {
        v.iter().map(|x| x.square_ref()).sum()
    }
}

impl VectorNorm for Integer {
    type Scalar = Float;

    fn norm<S>(v: &Vector<Self, S>) -> Self::Scalar
        where S: VectorStorage<Self> + ?Sized,
    {
        let ns = v.norm_sqr();
        Float::with_val(ns.signed_bits(), ns).sqrt()
    }
}

impl InnerProduct for Rational {
    fn dot<R, S>(v: &Vector<Self, R>, w: &Vector<Self, S>) -> Self
    where
        R: VectorStorage<Self> + ?Sized,
        S: VectorStorage<Self> + ?Sized
    {
        assert!(v.dim() == w.dim());
        v.iter().zip(w.iter())
            .map(|(c, d)| (c * d).complete())
            .sum()
    }

    fn norm_sqr<S>(v: &Vector<Self, S>) -> Self
        where S: VectorStorage<Self> + ?Sized,
    {
        v.iter().map(|x| x.square_ref().complete()).sum()
    }
}

impl VectorNorm for Rational {
    type Scalar = Float;
    fn norm<S>(v: &Vector<Self, S>) -> Self::Scalar
        where S: VectorStorage<Self> + ?Sized
    {
        let ns = v.norm_sqr();
        Float::with_val(ns.numer().signed_bits(), ns).sqrt()
    }
}

impl InnerProduct for Float {
    fn dot<R, S>(v: &Vector<Self, R>, w: &Vector<Self, S>) -> Self
    where
        R: VectorStorage<Self> + ?Sized,
        S: VectorStorage<Self> + ?Sized,
    {
        assert!(v.dim() == w.dim());
        let prec = v.precision();
        assert_eq!(prec, w.precision(), "Cannot compute the \
            inner product of two vectors of different precision.");
        v.iter().zip(w.iter())
            .map(|(a, b)| a * b)
            .fold(Float::with_val(prec, 0), |acc, x| acc + x)
    }

    fn norm_sqr<S>(v: &Vector<Self, S>) -> Self
        where S: VectorStorage<Self> + ?Sized,
    {
        let prec = v.precision();
        v.iter()
            .map(|x| x * x)
            .fold(Float::with_val(prec, 0), |acc, x| acc + x)
    }
}

impl VectorNorm for Float {
    type Scalar = Float;

    fn norm<S>(v: &Vector<Self, S>) -> Self::Scalar
        where S: VectorStorage<Self> + ?Sized,
    {
        v.norm_sqr().sqrt()
    }
}

macro_rules! impl_innerproduct_norm {
    ($t:ty) => {
        impl InnerProduct for $t {
            fn dot<R, S>(v: &Vector<Self, R>, w: &Vector<Self, S>) -> Self
            where
                R: VectorStorage<Self> + ?Sized,
                S: VectorStorage<Self> + ?Sized
            {
                assert!(v.dim() == w.dim());
                v.iter().zip(w.iter())
                    .map(|(c, d)| c * d)
                    .sum()
            }

            fn norm_sqr<S>(v: &Vector<Self, S>) -> Self
                where S: VectorStorage<Self> + ?Sized,
            {
                v.iter().map(|x| x * x).sum()
            }
        }

        impl VectorNorm for $t {
            type Scalar = Self;

            fn norm<S>(v: &Vector<Self, S>) -> Self::Scalar
                where S: VectorStorage<Self> + ?Sized,
            {
                v.norm_sqr().sqrt()
            }
        }
    }
}

impl_innerproduct_norm!(f32);
impl_innerproduct_norm!(f64);

macro_rules! impl_addsub {
    ($t:tt) => {
        impl_addsub!(impl, $t, add, Add);
        impl_addsub!(impl, $t, sub, Sub);
    };
    (impl, $t:tt, $op:tt, $class:tt) => {
        impl<R, S> $class<&Vector<$t, R>> for &Vector<$t, S>
        where
            R: VectorStorage<$t> + ?Sized,
            S: VectorStorage<$t> + ?Sized,
        {
            type Output = OwnedVector<$t>;
            fn $op(self, rhs: &Vector<$t, R>) -> Self::Output {
                assert!(self.dim() == rhs.dim(), "Can not perform operation \
                    for vectors of incompatible sizes");
                select!($t,
                    Float => {
                        let prec = self.precision();
                        assert!(prec == rhs.precision(),
                            "Can not add vectors of different precision.");
                    },
                    default => {},
                );
                OwnedVector::from_iter(self.dim(),
                    self.iter().zip(rhs.iter()).map(|(a, b)| select!($t,
                        Float => { Float::with_val(prec, <&$t>::$op(a, b)) },
                        Integer => { <&$t>::$op(a, b).complete() },
                        Rational => { <&$t>::$op(a, b).complete() },
                        default => { <&$t>::$op(a, b) },
                    ))
                )
            }
        }
    };
    ($t:tt, $($o:tt),+) => {
        impl_addsub!($t);
        impl_addsub!($($o),+);
    };
}

impl_addsub!(Integer, Rational, Float, f32, f64);

macro_rules! impl_addsub_reuse {
    ($t:ty) => {
        impl<S: VectorStorage<$t> + ?Sized> Add<&Vector<$t, S>> for OwnedVector<$t> {
            type Output = OwnedVector<$t>;
            fn add(mut self, rhs: &Vector<$t, S>) -> Self::Output {
                self += rhs;
                self
            }
        }

        impl<S: VectorStorage<$t> + ?Sized> Add<OwnedVector<$t>> for &Vector<$t, S> {
            type Output = OwnedVector<$t>;
            fn add(self, mut rhs: OwnedVector<$t>) -> Self::Output {
                rhs += self;
                rhs
            }
        }

        impl<S: VectorStorage<$t>  + ?Sized> Sub<&Vector<$t, S>> for OwnedVector<$t> {
            type Output = OwnedVector<$t>;
            fn sub(mut self, rhs: &Vector<$t, S>) -> Self::Output {
                self -= rhs;
                self
            }
        }

        impl<S: VectorStorage<$t> + ?Sized> Sub<OwnedVector<$t>> for &Vector<$t, S> {
            type Output = OwnedVector<$t>;
            fn sub(self, mut rhs: OwnedVector<$t>) -> Self::Output {
                for i in 0..self.dim() {
                    let ptr = &mut rhs[i] as *mut $t;
                    unsafe {
                        let r = ptr.read();
                        let v = &self[i] - r;
                        ptr.write(v);
                    }
                }
                rhs
            }
        }
    };
    ($t:ty, $($o:ty),+) => {
        impl_addsub_reuse!($t);
        impl_addsub_reuse!($($o),+);
    };
}

impl_addsub_reuse!(Integer, Rational, Float, f32, f64);

macro_rules! impl_muldiv {
    ($t:tt) => {
        impl_muldiv!(impl, $t, mul, Mul);
        impl_muldiv!(impl, $t, div, Div);
        impl_muldiv!(invert, $t);
    };
    (invert, $t:ty) => {
        impl<S: VectorStorage<$t> + ?Sized> Mul<&Vector<$t, S>> for &$t {
            type Output = OwnedVector<$t>;
            fn mul(self, rhs: &Vector<$t, S>) -> Self::Output {
                rhs * self
            }
        }
    };
    (impl, $t:tt, $op:tt, $class:tt) => {
        impl<S: VectorStorage<$t> + ?Sized> $class<&$t> for &Vector<$t, S> {
            type Output = OwnedVector<$t>;
            fn $op(self, rhs: &$t) -> Self::Output {
                select!($t,
                    Float => {
                        let prec = self.precision();
                            assert!(prec == rhs.prec(), "Can not \
                                multiply/divide vectors with float of \
                                different precision of different precision.");
                    },
                    default => {},
                );
                OwnedVector::from_iter(self.dim(),
                    self.iter().map(|a| select!($t,
                        Float => { Float::with_val(prec, <&$t>::$op(a, rhs)) },
                        Integer => { <&$t>::$op(a, rhs).complete() },
                        Rational => { <&$t>::$op(a, rhs).complete() },
                        default => { <&$t>::$op(a, rhs) },
                    ))
                )
            }
        }
    };
    ($t:tt, $($o:tt),+) => {
        impl_muldiv!($t);
        impl_muldiv!($($o),+);
    };
}

impl_muldiv!(Integer, Rational, Float, f32, f64);

macro_rules! check_prec {
    (Float, $l:expr, $r:expr) => {
        assert!($l.precision() == $r.precision(),
            "Can't add subtract vectors of different precision.");
    };
    ($_:tt, $l:expr, $r:expr) => {};
}

macro_rules! impl_assign_addsub {
    (impl, $t:tt, $op:tt, $class:tt) => {
        impl<R, S> $class<&Vector<$t, R>> for Vector<$t, S>
        where
            R: VectorStorage<$t> + ?Sized,
            S: VectorStorage<$t> + ?Sized,
        {
            fn $op(&mut self, rhs: &Vector<$t, R>) {
                assert!(self.dim() == rhs.dim(),
                    "Can not add/subtract vectors of different dimensions.");
                check_prec!($t, self, rhs);
                for i in 0..self.dim() {
                    macro_rules! op_impl {
                        (add_assign) => { self[i] += &rhs[i]; };
                        (sub_assign) => { self[i] -= &rhs[i]; };
                    }

                    op_impl!($op);
                }
            }
        }
    };
    ($t:tt, $($o:tt),+) => {
        impl_assign_addsub!($t);
        impl_assign_addsub!($($o),+);
    };
    ($t:tt) => {
        impl_assign_addsub!(impl, $t, add_assign, AddAssign);
        impl_assign_addsub!(impl, $t, sub_assign, SubAssign);
    };
}

impl_assign_addsub!(Integer, Rational, Float, f32, f64);

macro_rules! impl_assign_muldiv {
    ($t:ty) => {
        impl<S: VectorStorage<$t> + ?Sized> MulAssign<&$t> for Vector<$t, S> {
            fn mul_assign(&mut self, rhs: &$t) {
                check_prec!($t, self, rhs);
                self.map_mut(|i| *i *= rhs);
            }
        }

        impl<S: VectorStorage<$t> + ?Sized> DivAssign<&$t> for Vector<$t, S> {
            fn div_assign(&mut self, rhs: &$t) {
                check_prec!($t, self, rhs);
                self.map_mut(|i| *i /= rhs);
            }
        }
    };
    ($t:ty, $($o:ty),+) => {
        impl_assign_muldiv!($t);
        impl_assign_muldiv!($($o),+);
    };
}

impl_assign_muldiv!(Integer, Rational, Float, f32, f64);

macro_rules! impl_neg {
    ($t:ty) => {
        impl Neg for OwnedVector<$t> {
            type Output = OwnedVector<$t>;
            fn neg(mut self) -> Self::Output {
                for e in self.iter_mut() {
                    e.neg_assign();
                }
                self
            }
        }
    };
    ($t:ty, $($o:ty),+) => {
        impl_neg!($t);
        impl_neg!($($o),+);
    }
}

impl_neg!(Integer, Rational, Float, f32, f64);