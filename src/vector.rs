//! Owned and non-owned vectors constant runtime dimension.

use std::{ops::{Deref, DerefMut}, fmt::Debug};
use num_traits::Zero;
use rug::{Integer, Complete, Float};

/// Vector view.
#[derive(Debug)]
#[repr(transparent)]
pub struct VV<T>([T]);

pub type IVV = VV<Integer>;
pub type FVV = VV<Float>;

impl<T> VV<T> {
    /// Get a vector view from a slice.
    pub fn from_slice<'a>(s: &'a [T]) -> &'a Self {
        unsafe {
            &*std::ptr::from_raw_parts(s.as_ptr() as _, s.len())
        }
    }

    /// Get a mutable vector view from a slice.
    pub fn from_slice_mut<'a>(s: &'a mut [T]) -> &'a mut Self {
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(s.as_mut_ptr() as _, s.len())
        }
    }

    /// Pointer to the entries.
    pub fn as_ptr(&self) -> *const T {
        self.0.as_ptr()
    }

    /// Mutable pointer to the entries.
    pub fn as_mut_ptr(&mut self) -> *mut T {
        self.0.as_mut_ptr()
    }

    /// The dimension of the vector.
    pub fn dim(&self) -> usize {
        self.0.len()
    }

    /// Is the vector empty, i.e. dimension zero?
    pub fn is_empty(&self) -> bool {
        self.dim() == 0
    }

    /// Interpret the vector as a slice.
    pub fn as_slice(&self) -> &[T] {
        &self.0
    }

    /// Interpret the vector as a mutable slice.
    pub fn as_slice_mut(&mut self) -> &mut [T] {
        &mut self.0
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.as_slice().iter()
    }

    /// Returns an iterator over the mutable elements.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, T> {
        self.as_slice_mut().iter_mut()
    }

    /// Reference to an entry at an index.
    pub fn entry(&self, idx: usize) -> &T {
        assert!(idx < self.dim());
        unsafe {
            &*self.as_ptr().add(idx)
        }
    }

    /// Mutable reference to an entry at an index.
    pub fn entry_mut(&mut self, idx: usize) -> &mut T {
        assert!(idx < self.dim());
        unsafe {
            &mut *self.as_mut_ptr().add(idx)
        }
    }

    /// Apply a function to each entry.
    pub fn map_mut<F: FnMut(&mut T)>(&mut self, mut f: F) {
        for e in self.iter_mut() {
            f(e);
        }
    }

    /// Swap the two entries i and j.
    pub fn swap(&mut self, i : usize, j: usize) {
        unsafe {
            std::ptr::swap(
                self.as_mut_ptr().add(i),
                self.as_mut_ptr().add(j)
            );
        }
    }
}

impl<T: Zero> VV<T> {
    /// Is this the zero vector.
    pub fn is_zero(&self) -> bool {
        self.iter().all(|i| i.is_zero())
    }
}

impl VV<Integer> {
    /// Computes the dot product of two vectors.
    pub fn dot(&self, other: &Self) -> Integer {
        assert!(self.dim() == other.dim());
        self.iter().zip(other.iter())
            .map(|(c, d)| (c * d).complete())
            .sum()
    }

    /// Computes the square of the l2 norm of the vector.
    pub fn norm_sqr(&self) -> Integer {
        self.iter().map(|i| i.square_ref()).sum()
    }

    /// Computes the l2 norm of the vector.
    pub fn norm(&self) -> Float {
        let ns = self.norm_sqr();
        Float::with_val(ns.signed_bits(), ns).sqrt()
    }
}

impl<'a, T> IntoIterator for &'a VV<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a, T> IntoIterator for &'a mut VV<T> {
    type Item = &'a mut T;
    type IntoIter = std::slice::IterMut<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl<T> std::ops::Index<usize> for VV<T> {
    type Output = T;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T> std::ops::IndexMut<usize> for VV<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T> AsRef<[T]> for VV<T> {
    fn as_ref(&self) -> &[T] {
        self.as_slice()
    }
}

impl<T> AsMut<[T]> for VV<T> {
    fn as_mut(&mut self) -> &mut [T] {
        self.as_slice_mut()
    }
}

/// An owned vector whose entries are integers.
pub struct Vector<T> {
    /// Memory that holds the entries.
    pub(self) entries: *mut T,

    /// The dimension (number of entries) of the vector.
    pub dim: usize,

}

pub type IVector = Vector<Integer>;
pub type FVector = Vector<Float>;

impl<T> Vector<T> {
    /// Returns an empty vector.
    pub fn empty() -> Self {
        Self {
            dim: 0,
            entries: core::ptr::null_mut(),
        }
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

        Self {
            dim,
            entries,
        }
    }

    /// Creates a vector from a slice.
    pub fn from_entries<U, V>(a: U) -> Self
    where
        U: AsRef<[V]>,
        V: Into<T> + Clone,
    {
        let a = a.as_ref();
        let v = Self::uninit(a.len());
        for i in 0..v.dim {
            unsafe {
                v.entries.add(i).write(a[i].clone().into());
            }
        }

        v
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
                v.entries.add(i).write(e);
            }
        }

        v
    }

    /// Creates a vector from an array.
    pub fn from_array<U: Into<T>, const D: usize>(a: [U; D]) -> Self {
        let v = Self::uninit(D);

        let mut ptr = v.entries;

        for e in a {
            unsafe {
                ptr.write(e.into());
                ptr = ptr.add(1);
            }
        }

        v
    }

    /// Free the entries.
    pub(self) unsafe fn free(&mut self) {
        let layout = std::alloc::Layout::from_size_align(
            self.dim * core::mem::size_of::<Integer>(),
            core::mem::align_of::<Integer>()
        ).unwrap();

        // Free the memory.
        std::alloc::dealloc(self.entries as _, layout);
    }

    pub fn map<Fn: FnMut(&mut T)>(mut self, f: Fn) -> Self {
        self.map_mut(f);
        self
    }
}

impl Vector<Integer> {
        /// Returns a zero vector.
    pub fn zero(dim: usize) -> Self {
        let v = Self::uninit(dim);

        for i in 0..dim {
            unsafe {
                v.entries.add(i).write(Integer::new());
            }
        }

        v
    }

    /// Returns a vector with index `i` set to `v`
    /// and every other entry to zero.
    pub fn ith(dim: usize, i: usize, v: Integer) -> Self {
        assert!(i < dim);
        let mut r = Self::zero(dim);
        r[i] = v;
        r
    }
}

impl<T> Deref for Vector<T> {
    type Target = VV<T>;

    fn deref(&self) -> &Self::Target {
        let slice = unsafe {
            std::slice::from_raw_parts(self.entries, self.dim)
        };
        VV::from_slice(slice)
    }
}

impl<T> DerefMut for Vector<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        let slice = unsafe {
            std::slice::from_raw_parts_mut(self.entries, self.dim)
        };
        VV::from_slice_mut(slice)
    }
}

impl<T> std::ops::Index<usize> for Vector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        self.entry(index)
    }
}

impl<T> std::ops::IndexMut<usize> for Vector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.entry_mut(index)
    }
}

impl<T: PartialEq> PartialEq for Vector<T> {
    fn eq(&self, other: &Self) -> bool {
        assert!(self.dim == other.dim, "Can not compare vectors of different sizes");
        self.iter()
            .zip(other.iter())
            .all(|(l, r)| l == r)
    }
}

impl<T: Eq> Eq for Vector<T> {}

impl<T: Clone> Clone for Vector<T> {
    fn clone(&self) -> Self {
        let v = Self::uninit(self.dim);

        for i in 0..self.dim {
            unsafe {
                v.entries.add(i).write(self[i].clone());
            }
        }

        return v;
    }
}

impl<T: Debug> Debug for Vector<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries(self.as_slice())
            .finish()
    }
}

impl<T> AsRef<[T]> for Vector<T> {
    fn as_ref(&self) -> &[T] {
        self.as_slice()
    }
}

impl<T> AsMut<[T]> for Vector<T> {
    fn as_mut(&mut self) -> &mut [T] {
        self.as_slice_mut()
    }
}

impl<T> Drop for Vector<T> {
    fn drop(&mut self) {
        if self.entries.is_null() {
            return;
        }

        unsafe {
            // Call the destructor on each element.
            for i in 0..self.dim {
                core::ptr::drop_in_place(self.entries.add(i));
            }

            self.free();
        }
    }
}


macro_rules! impl_addsub {
    ($t:ty, $u:ty) => {
        impl std::ops::Add<$u> for $t {
            type Output = Vector<Integer>;
            fn add(self, rhs: $u) -> Self::Output {
                assert!(self.dim() == rhs.dim(),
                    "Can not add vectors of incompatible sizes");
                Vector::from_iter(self.dim(),
                    self.iter().zip(rhs.iter()).map(|(a, b)| (a + b).complete())
                )
            }
        }

        impl std::ops::Sub<$u> for $t {
            type Output = Vector<Integer>;
            fn sub(self, rhs: $u) -> Self::Output {
                assert!(self.dim() == rhs.dim(),
                    "Can not subtract vectors of incompatible sizes");
                Vector::from_iter(self.dim(),
                    self.iter().zip(rhs.iter()).map(|(a, b)| (a - b).complete())
                )
            }
        }
    }
}

impl_addsub!(&Vector<Integer>, &Vector<Integer>);
impl_addsub!(&Vector<Integer>, &VV<Integer>);
impl_addsub!(&VV<Integer>, &Vector<Integer>);
impl_addsub!(&VV<Integer>, &VV<Integer>);

macro_rules! impl_addsub_reuse {
    ($t:ty) => {
        impl std::ops::Add<$t> for Vector<Integer> {
            type Output = Vector<Integer>;
            fn add(mut self, rhs: $t) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                self += rhs;
                self
            }
        }

        impl std::ops::Add<Vector<Integer>> for $t {
            type Output = Vector<Integer>;
            fn add(self, mut rhs: Vector<Integer>) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                rhs += self;
                rhs
            }
        }

        impl std::ops::Sub<$t> for Vector<Integer> {
            type Output = Vector<Integer>;
            fn sub(mut self, rhs: $t) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                self -= rhs;
                self
            }
        }

        impl std::ops::Sub<Vector<Integer>> for $t {
            type Output = Vector<Integer>;
            fn sub(self, mut rhs: Vector<Integer>) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                // We can't get away without any allocation because
                // sub is not commutative.
                for i in 0..self.dim() {
                    let s = &self[i] - &rhs[i];
                    let s = s.complete();
                    rhs[i] = s;
                }
                rhs
            }
        }
    }
}

impl_addsub_reuse!(&Vector<Integer>);
impl_addsub_reuse!(&VV<Integer>);

macro_rules! impl_muldiv {
    ($t:ty) => {
        impl std::ops::Mul<&Integer> for $t {
            type Output = Vector<Integer>;
            fn mul(self, rhs: &Integer) -> Self::Output {
                Vector::from_iter(self.dim(),
                    self.iter().map(|i| (i * rhs).complete())
                )
            }
        }

        impl std::ops::Mul<$t> for &Integer {
            type Output = Vector<Integer>;
            fn mul(self, rhs: $t) -> Self::Output {
                Vector::from_iter(rhs.dim(),
                    rhs.iter().map(|i| (i * self).complete())
                )
            }
        }

        impl std::ops::Div<&Integer> for $t {
            type Output = Vector<Integer>;
            fn div(self, rhs: &Integer) -> Self::Output {
                Vector::from_iter(self.dim(),
                    self.iter().map(|i| (i / rhs).complete())
                )
            }
        }
    }
}

impl_muldiv!(&Vector<Integer>);
impl_muldiv!(&VV<Integer>);

macro_rules! impl_assign_addsub {
    ($t:ty, $u:ty) => {
        impl std::ops::AddAssign<$u> for $t {
            fn add_assign(&mut self, rhs: $u) {
                assert!(self.dim() == rhs.dim(), "Can not add vectors of different dimensions.");
                for i in 0..self.dim() {
                    self[i] += &rhs[i];
                }
            }
        }

        impl std::ops::SubAssign<$u> for $t {
            fn sub_assign(&mut self, rhs: $u) {
                assert!(self.dim() == rhs.dim(), "Can not subtract vectors of different dimensions.");
                for i in 0..self.dim() {
                    self[i] -= &rhs[i];
                }
            }
        }
    }
}

impl_assign_addsub!(Vector<Integer>, &Vector<Integer>);
impl_assign_addsub!(Vector<Integer>, &VV<Integer>);
impl_assign_addsub!(VV<Integer>, &Vector<Integer>);
impl_assign_addsub!(VV<Integer>, &VV<Integer>);

macro_rules! impl_assign_muldiv {
    ($t:ty) => {
        impl std::ops::MulAssign<&Integer> for $t {
            fn mul_assign(&mut self, rhs: &Integer) {
                self.map_mut(|i| *i *= rhs);
            }
        }

        impl std::ops::DivAssign<&Integer> for $t {
            fn div_assign(&mut self, rhs: &Integer) {
                self.map_mut(|i| *i /= rhs);
            }
        }
    }
}

impl_assign_muldiv!(Vector<Integer>);
impl_assign_muldiv!(&mut VV<Integer>);

impl std::ops::Neg for Vector<Integer> {
    type Output = Vector<Integer>;
    fn neg(mut self) -> Self::Output {
        for e in self.iter_mut() {
            *e *= -1;
        }

        self
    }
}