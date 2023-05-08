//! Owned vector with integer entries and
//! trait to interpret any slice of integers as one.

use num_traits::Zero;
use rug::{Integer, Complete, Float};

/// Operations implemented for an integer vector.
pub trait VectorOps {
    /// The dimension (number of entries) of the vector.
    fn dim(&self) -> usize;

    /// Pointer to the entries.
    fn as_ptr(&self) -> *const Integer;

    /// Pointer to mutable entries.
    fn as_mut_ptr(&mut self) -> *mut Integer;

    /// Is this vector of dimension zero?
    fn is_empty(&self) -> bool {
        self.dim() == 0
    }

    /// Interpret the vector as a slice.
    fn as_slice(&self) -> &[Integer] {
        unsafe {
            std::slice::from_raw_parts(self.as_ptr(), self.dim())
        }
    }

    /// Interpret the vector as a mutable slice.
    fn as_slice_mut(&mut self) -> &mut [Integer] {
        unsafe {
            std::slice::from_raw_parts_mut(self.as_mut_ptr(), self.dim())
        }
    }

    /// Returns an iterator over the elements.
    fn iter(&self) -> std::slice::Iter<'_, Integer> {
        self.as_slice().iter()
    }

    /// Returns an iterator over the mutable elements.
    fn iter_mut(&mut self) -> std::slice::IterMut<'_, Integer> {
        self.as_slice_mut().iter_mut()
    }

    /// Is this the zero vector.
    fn is_zero(&self) -> bool {
        self.iter().all(|i| i.is_zero())
    }

    /// Reference to an entry at an index.
    fn entry(&self, idx: usize) -> &Integer {
        assert!(idx < self.dim());
        unsafe {
            &*self.as_ptr().add(idx)
        }
    }

    /// Mutable reference to an entry at an index.
    fn entry_mut(&mut self, idx: usize) -> &mut Integer {
        assert!(idx < self.dim());
        unsafe {
            &mut *self.as_mut_ptr().add(idx)
        }
    }

    /// Computes the dot product of two vectors.
    fn dot(&self, other: &Self) -> Integer {
        assert!(self.dim() == other.dim());
        self.iter().zip(other.iter())
            .map(|(c, d)| (c * d).complete())
            .sum()
    }

    /// Computes the square of the l2 norm of the vector.
    fn norm_sqr(&self) -> Integer {
        self.iter().map(|i| i.square_ref()).sum()
    }

    /// Computes the l2 norm of the vector.
    fn norm(&self) -> Float {
        let ns = self.norm_sqr();
        Float::with_val(ns.signed_bits(), ns).sqrt()
    }

    /// Apply a function to each entry.
    fn map_mut<F: FnMut(&mut Integer)>(&mut self, mut f: F) {
        for e in self.iter_mut() {
            f(e);
        }
    }

    /// Swap the two entries i and j.
    fn swap(&mut self, i : usize, j: usize) {
        unsafe {
            std::ptr::swap(
                self.as_mut_ptr().add(i),
                self.as_mut_ptr().add(j)
            );
        }
    }
}

#[derive(Debug)]
#[repr(transparent)]
pub struct VV([Integer]);

impl VV {
    pub fn from_slice<'a>(s: &'a [Integer]) -> &'a Self {
        unsafe {
            &*std::ptr::from_raw_parts(s.as_ptr() as _, s.len())
        }
    }

    pub fn from_slice_mut<'a>(s: &'a mut [Integer]) -> &'a mut Self {
        unsafe {
            &mut *std::ptr::from_raw_parts_mut(s.as_mut_ptr() as _, s.len())
        }
    }
}

impl VectorOps for VV {
    fn as_ptr(&self) -> *const Integer {
        self.0.as_ptr()
    }

    fn as_mut_ptr(&mut self) -> *mut Integer {
        self.0.as_mut_ptr()
    }

    fn dim(&self) -> usize {
        self.0.len()
    }
}

impl<'a> IntoIterator for &'a VV {
    type Item = &'a Integer;
    type IntoIter = std::slice::Iter<'a, Integer>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

impl<'a> IntoIterator for &'a mut VV {
    type Item = &'a mut Integer;
    type IntoIter = std::slice::IterMut<'a, Integer>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter_mut()
    }
}

impl std::ops::Index<usize> for VV {
    type Output = Integer;
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for VV {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl AsRef<[Integer]> for VV {
    fn as_ref(&self) -> &[Integer] {
        self.as_slice()
    }
}

impl AsMut<[Integer]> for VV {
    fn as_mut(&mut self) -> &mut [Integer] {
        self.as_slice_mut()
    }
}

/// An owned vector whose entries are integers.
pub struct Vector {
    /// The dimension (number of entries) of the vector.
    pub dim: usize,

    /// Memory that holds the entries.
    pub(self) entries: *mut Integer,
}

impl Vector {
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
            dim * core::mem::size_of::<Integer>(),
            core::mem::align_of::<Integer>()
        ).unwrap();

        let entries = unsafe {
            std::alloc::alloc(layout) as *mut Integer
        };

        Self {
            dim,
            entries,
        }
    }

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
    pub fn ith(dim: usize, i: usize, v: Integer) -> Vector {
        assert!(i < dim);
        let mut r = Self::zero(dim);
        r[i] = v;
        r
    }

    /// Creates a vector from a slice.
    pub fn from_entries<T, U>(a: T) -> Self
    where
        T: AsRef<[U]>,
        U: Into<Integer> + Clone,
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
    pub fn from_iter<T: Iterator<Item = Integer>>(
        dim: usize, mut iter: T
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
    pub fn from_array<T: Into<Integer>, const D: usize>(a: [T; D]) -> Self {
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

    pub fn map<Fn: FnMut(&mut Integer)>(mut self, f: Fn) -> Self {
        self.map_mut(f);
        self
    }
}

impl VectorOps for Vector {
    fn dim(&self) -> usize {
        self.dim
    }

    fn as_ptr(&self) -> *const Integer {
        self.entries
    }

    fn as_mut_ptr(&mut self) -> *mut Integer {
        self.entries
    }
}

impl std::ops::Index<usize> for Vector {
    type Output = Integer;

    fn index(&self, index: usize) -> &Self::Output {
        self.entry(index)
    }
}

impl std::ops::IndexMut<usize> for Vector {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.entry_mut(index)
    }
}

impl PartialEq for Vector {
    fn eq(&self, other: &Self) -> bool {
        assert!(self.dim == other.dim, "Can not compare vectors of different sizes");
        self.iter()
            .zip(other.iter())
            .all(|(l, r)| l == r)
    }
}

impl Eq for Vector {}

impl Clone for Vector {
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

impl std::fmt::Debug for Vector {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries(self.as_slice())
            .finish()
    }
}

impl AsRef<[Integer]> for Vector {
    fn as_ref(&self) -> &[Integer] {
        self.as_slice()
    }
}

impl AsMut<[Integer]> for Vector {
    fn as_mut(&mut self) -> &mut [Integer] {
        self.as_slice_mut()
    }
}

impl Drop for Vector {
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
            type Output = Vector;
            fn add(self, rhs: $u) -> Self::Output {
                assert!(self.dim() == rhs.dim(),
                    "Can not add vectors of incompatible sizes");
                Vector::from_iter(self.dim(),
                    self.iter().zip(rhs.iter()).map(|(a, b)| (a + b).complete())
                )
            }
        }

        impl std::ops::Sub<$u> for $t {
            type Output = Vector;
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

impl_addsub!(&Vector, &Vector);
impl_addsub!(&Vector, &VV);
impl_addsub!(&VV, &Vector);
impl_addsub!(&VV, &VV);

macro_rules! impl_addsub_reuse {
    ($t:ty) => {
        impl std::ops::Add<$t> for Vector {
            type Output = Vector;
            fn add(mut self, rhs: $t) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                self += rhs;
                self
            }
        }

        impl std::ops::Add<Vector> for $t {
            type Output = Vector;
            fn add(self, mut rhs: Vector) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                rhs += self;
                rhs
            }
        }

        impl std::ops::Sub<$t> for Vector {
            type Output = Vector;
            fn sub(mut self, rhs: $t) -> Self::Output {
                assert!(self.dim() == rhs.dim());
                self -= rhs;
                self
            }
        }

        impl std::ops::Sub<Vector> for $t {
            type Output = Vector;
            fn sub(self, mut rhs: Vector) -> Self::Output {
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

impl_addsub_reuse!(&Vector);
impl_addsub_reuse!(&VV);

macro_rules! impl_muldiv {
    ($t:ty) => {
        impl std::ops::Mul<&Integer> for $t {
            type Output = Vector;
            fn mul(self, rhs: &Integer) -> Self::Output {
                Vector::from_iter(self.dim(),
                    self.iter().map(|i| (i * rhs).complete())
                )
            }
        }

        impl std::ops::Mul<$t> for &Integer {
            type Output = Vector;
            fn mul(self, rhs: $t) -> Self::Output {
                Vector::from_iter(rhs.dim(),
                    rhs.iter().map(|i| (i * self).complete())
                )
            }
        }

        impl std::ops::Div<&Integer> for $t {
            type Output = Vector;
            fn div(self, rhs: &Integer) -> Self::Output {
                Vector::from_iter(self.dim(),
                    self.iter().map(|i| (i / rhs).complete())
                )
            }
        }
    }
}

impl_muldiv!(&Vector);
impl_muldiv!(&VV);

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

impl_assign_addsub!(Vector, &Vector);
impl_assign_addsub!(Vector, &VV);
impl_assign_addsub!(VV, &Vector);
impl_assign_addsub!(VV, &VV);

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

impl_assign_muldiv!(Vector);
impl_assign_muldiv!(&mut VV);

impl std::ops::Neg for Vector {
    type Output = Vector;
    fn neg(mut self) -> Self::Output {
        for e in self.iter_mut() {
            *e *= -1;
        }

        self
    }
}