//! Owned and non-owned vectors constant runtime dimension.

use std::borrow::{Borrow, BorrowMut};
use std::ops::{Deref, DerefMut, Neg, Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign, Index, Range, IndexMut};
use std::fmt::Debug;
use num_traits::Zero;
use rug::{Integer, Complete, Float, Rational, ops::NegAssign};

/// Vector view.
#[derive(Debug)]
#[repr(transparent)]
pub struct VV<T>([T]);

#[allow(clippy::upper_case_acronyms)]
pub type IVV = VV<Integer>;
#[allow(clippy::upper_case_acronyms)]
pub type QVV = VV<Rational>;
#[allow(clippy::upper_case_acronyms)]
pub type FVV = VV<Float>;

impl<T> VV<T> {
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
    pub fn map<U, F: FnMut(&T) -> U>(&self, f: F) -> Vector<U> {
        Vector::from_iter(self.dim(), self.iter().map(f))
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

impl<T> VV<T>
    where for<'a> Float: rug::Assign<&'a T>
{
    pub fn to_float(&self, prec: u32) -> FVector {
        self.map(|e| Float::with_val(prec, e))
    }
}

impl IVV {
    /// Computes the dot product of two vectors.
    pub fn dot(&self, other: &Self) -> Integer {
        assert!(self.dim() == other.dim());
        self.iter().zip(other.iter())
            .map(|(c, d)| c * d)
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

impl FVV {
    pub fn precision(&self) -> u32 {
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

    /// Computes the dot product of two vectors.
    pub fn dot(&self, other: &Self) -> Float {
        assert!(self.dim() == other.dim());
        self.iter().zip(other.iter())
            .map(|(c, d)| c * d)
            .fold(Float::with_val(self.precision(), 0), |acc, f| acc + f)
    }

    /// Computes the square of the l2 norm of the vector.
    pub fn norm_sqr(&self) -> Float {
        self.iter().map(|i| i * i)
            .fold(Float::with_val(self.precision(), 0), |acc, f| acc + f)
    }

    /// Computes the l2 norm of the vector.
    pub fn norm(&self) -> Float {
        self.norm_sqr().sqrt()
    }
}

impl<T> Index<Range<usize>> for VV<T> {
    type Output = VV<T>;
    fn index(&self, index: Range<usize>) -> &Self::Output {
        VV::from_slice(&self.as_slice()[index])
    }
}

impl<T> IndexMut<Range<usize>> for VV<T> {
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        VV::from_slice_mut(&mut self.as_slice_mut()[index])
    }
}

impl<T: Clone> ToOwned for VV<T> {
    type Owned = Vector<T>;
    fn to_owned(&self) -> Self::Owned {
        Vector::from_iter(self.dim(), self.iter().cloned())
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

impl<T, U: ?Sized> PartialEq<U> for VV<T>
    where
        [T]: PartialEq<U>,
{
    fn eq(&self, other: &U) -> bool {
        self.as_slice().eq(other)
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
pub type QVector = Vector<Rational>;
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

    pub fn view(&self) -> &VV<T> {
        let slice = unsafe {
            std::slice::from_raw_parts(self.entries, self.dim)
        };
        VV::from_slice(slice)
    }

    pub fn view_mut(&mut self) -> &mut VV<T> {
        let slice = unsafe {
            std::slice::from_raw_parts_mut(self.entries, self.dim)
        };
        VV::from_slice_mut(slice)
    }

    /// Free the entries.
    pub(self) unsafe fn free(&mut self) {
        let layout = std::alloc::Layout::from_size_align(
            self.dim * core::mem::size_of::<T>(),
            core::mem::align_of::<T>()
        ).unwrap();

        // Free the memory.
        std::alloc::dealloc(self.entries as _, layout);
    }

    pub fn map<Fn: FnMut(&mut T)>(mut self, f: Fn) -> Self {
        self.map_mut(f);
        self
    }

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

impl<T: Zero> Vector<T> {
    /// Returns a zero vector.
    pub fn zero(dim: usize) -> Self {
        let v = Self::uninit(dim);

        for i in 0..dim {
            unsafe {
                v.entries.add(i).write(T::zero());
            }
        }

        v
    }

    /// Returns a vector with index `i` set to `v`
    /// and every other entry to zero.
    pub fn ith(dim: usize, i: usize, v: T) -> Self {
        assert!(i < dim);
        let mut r = Self::zero(dim);
        r[i] = v;
        r
    }
}

impl FVector {
    /// Zero vector with a certain precision.
    pub fn zero_prec(dim: usize, prec: u32) -> Self {
        Self::from_iter(dim, std::iter::repeat(Float::with_val(prec, 0)))
    }
}

impl<T> Borrow<VV<T>> for Vector<T> {
    fn borrow(&self) -> &VV<T> {
        self.view()
    }
}

impl<T> BorrowMut<VV<T>> for Vector<T> {
    fn borrow_mut(&mut self) -> &mut VV<T> {
        self.view_mut()
    }
}

impl<T> Deref for Vector<T> {
    type Target = VV<T>;

    fn deref(&self) -> &Self::Target {
        self.view()
    }
}

impl<T> DerefMut for Vector<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        self.view_mut()
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

        v
    }
}

impl<T: Debug> Debug for Vector<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries(self.as_slice())
            .finish()
    }
}

impl<T> AsRef<VV<T>> for Vector<T> {
    fn as_ref(&self) -> &VV<T> {
        self.view()
    }
}

impl<T> AsMut<VV<T>> for Vector<T> {
    fn as_mut(&mut self) -> &mut VV<T> {
        self.view_mut()
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
    ($t:tt) => {
        impl_addsub!(impl, $t, add, Add);
        impl_addsub!(impl, $t, sub, Sub);
    };
    (impl, $t:tt, $op:tt, $class:tt) => {
        impl_addsub!(impl, &Vector<$t>,   &Vector<$t>,    $t, $op, $class);
        impl_addsub!(impl, &Vector<$t>,   &VV<$t>,        $t, $op, $class);
        impl_addsub!(impl, &VV<$t>,       &Vector<$t>,    $t, $op, $class);
        impl_addsub!(impl, &VV<$t>,       &VV<$t>,        $t, $op, $class);
    };
    (impl, $t:ty, $u:ty, $v:tt, $op:tt, $class:tt) => {
        impl $class<$u> for $t {
            type Output = Vector<$v>;
            fn $op(self, rhs: $u) -> Self::Output {
                assert!(self.dim() == rhs.dim(), "Can not perform operation \
                    for vectors of incompatible sizes");
                macro_rules! prec {
                    (Float) => {
                        {
                            let prec = self.precision();
                            assert!(prec == rhs.precision(),
                                "Can not add vectors of different precision.");
                            prec
                        }
                    };
                    ($_:tt) => { 0u32 };
                };
                let prec = prec!($v);
                macro_rules! op_impl {
                    ($x:tt, $opa:tt) => { op_impl!($x, $x, $opa) };
                    (default, $x:ty, $opa:tt) => {
                        |(a, b): (&$x, &$x)| <&$x>::$opa(a, b).complete()
                    };
                    (Integer, $x:ty, $opa:tt) => { op_impl!(default, $x, $opa) };
                    (Rational, $x:ty, $opa:tt) => { op_impl!(default, $x, $opa) };
                    (Float, $x:ty, $opa:tt) => {
                        |(a, b): (&$x, &$x)| Float::with_val(prec, <&$x>::$opa(a, b))
                    };
                }
                Vector::from_iter(self.dim(),
                    self.iter().zip(rhs.iter()).map(op_impl!($v, $op)))
            }
        }
    };
    ($t:tt, $($o:tt),+) => {
        impl_addsub!($t);
        impl_addsub!($($o),+);
    };
}

impl_addsub!(Integer, Rational, Float);

macro_rules! impl_addsub_reuse {
    (impl, $t:ty, $u:ty) => {
        impl Add<$t> for Vector<$u> {
            type Output = Vector<$u>;
            fn add(mut self, rhs: $t) -> Self::Output {
                self += rhs;
                self
            }
        }

        impl Add<Vector<$u>> for $t {
            type Output = Vector<$u>;
            fn add(self, mut rhs: Vector<$u>) -> Self::Output {
                rhs += self;
                rhs
            }
        }

        impl Sub<$t> for Vector<$u> {
            type Output = Vector<$u>;
            fn sub(mut self, rhs: $t) -> Self::Output {
                self -= rhs;
                self
            }
        }

        impl Sub<Vector<$u>> for $t {
            type Output = Vector<$u>;
            fn sub(self, mut rhs: Vector<$u>) -> Self::Output {
                // If we were to do this in one step,
                // we would need to allocate a temporary object.
                // Might be worth it, not too sure.
                rhs = -rhs;
                rhs + self
            }
        }
    };
    ($t:ty) => {
        impl_addsub_reuse!(impl, &Vector<$t>, $t);
        impl_addsub_reuse!(impl, &VV<$t>, $t);
    };
    ($t:ty, $($o:ty),+) => {
        impl_addsub_reuse!($t);
        impl_addsub_reuse!($($o),+);
    };
}

impl_addsub_reuse!(Integer, Rational, Float);

macro_rules! impl_muldiv {
    ($t:tt) => {
        impl_muldiv!(impl, $t, mul, Mul);
        impl_muldiv!(impl, $t, div, Div);
        impl_muldiv!(invert, $t, &Vector<$t>);
        impl_muldiv!(invert, $t, &VV<$t>);
    };
    (impl, $t:tt, $op:tt, $class:tt) => {
        impl_muldiv!(impl, &Vector<$t>,   $t, $op, $class);
        impl_muldiv!(impl, &VV<$t>,       $t, $op, $class);
    };
    (invert, $t:ty, $v:ty) => {
        impl Mul<$v> for &$t {
            type Output = Vector<$t>;
            fn mul(self, rhs: $v) -> Self::Output {
                rhs * self
            }
        }
    };
    (impl, $v:ty, $t:tt, $op:tt, $class:tt) => {
        impl $class<&$t> for $v {
            type Output = Vector<$t>;
            fn $op(self, rhs: &$t) -> Self::Output {
                macro_rules! prec {
                    (Float) => {
                        {
                            let prec = self.precision();
                            assert!(prec == rhs.prec(), "Can not \
                                multiply/divide vectors with float of \
                                different precision of different precision.");
                            prec
                        }
                    };
                    ($_:tt) => { 0u32 };
                };
                let prec = prec!($t);
                macro_rules! op_impl {
                    ($x:tt, $opa:tt) => { op_impl!($x, $x, $opa) };
                    (default, $x:ty, $opa:tt) => {
                        |i: &$x| <&$x>::$opa(i, rhs).complete()
                    };
                    (Integer, $x:ty, $opa:tt) => { op_impl!(default, $x, $opa) };
                    (Rational, $x:ty, $opa:tt) => { op_impl!(default, $x, $opa) };
                    (Float, $x:ty, $opa:tt) => {
                        |i: &$x| Float::with_val(prec, <&$x>::$opa(i, rhs))
                    };
                }
                Vector::from_iter(self.dim(),
                    self.iter().map(op_impl!($t, $op)))
            }
        }
    };
    ($t:tt, $($o:tt),+) => {
        impl_muldiv!($t);
        impl_muldiv!($($o),+);
    };
}

impl_muldiv!(Integer, Rational, Float);

macro_rules! check_prec {
    (Float, $l:expr, $r:expr) => {
        assert!($l.precision() == $r.precision(),
            "Can't add subtract vectors of different precision.");
    };
    ($_:tt, $l:expr, $r:expr) => {};
}

macro_rules! impl_assign_addsub {
    (impl, $i:tt, $t:ty, $u:ty) => {
        impl AddAssign<$u> for $t {
            fn add_assign(&mut self, rhs: $u) {
                assert!(self.dim() == rhs.dim(), "Can not add vectors of different dimensions.");
                check_prec!($i, self, rhs);
                for i in 0..self.dim() {
                    self[i] += &rhs[i];
                }
            }
        }

        impl SubAssign<$u> for $t {
            fn sub_assign(&mut self, rhs: $u) {
                assert!(self.dim() == rhs.dim(), "Can not subtract vectors of different dimensions.");
                check_prec!($i, self, rhs);
                for i in 0..self.dim() {
                    self[i] -= &rhs[i];
                }
            }
        }
    };
    ($t:tt) => {
        impl_assign_addsub!(impl, $t, Vector<$t>, &Vector<$t>);
        impl_assign_addsub!(impl, $t, Vector<$t>, &VV<$t>);
        impl_assign_addsub!(impl, $t, VV<$t>, &Vector<$t>);
        impl_assign_addsub!(impl, $t, VV<$t>, &VV<$t>);
    };
    ($t:tt, $($o:tt),+) => {
        impl_assign_addsub!($t);
        impl_assign_addsub!($($o),+);
    };
}

impl_assign_addsub!(Integer, Rational, Float);

macro_rules! impl_assign_muldiv {
    (impl, $t:ty, $v:ty) => {
        impl MulAssign<&$t> for $v {
            fn mul_assign(&mut self, rhs: &$t) {
                check_prec!($t, self, rhs);
                self.map_mut(|i| *i *= rhs);
            }
        }

        impl DivAssign<&$t> for $v {
            fn div_assign(&mut self, rhs: &$t) {
                check_prec!($t, self, rhs);
                self.map_mut(|i| *i /= rhs);
            }
        }
    };
    ($t:ty, $($o:ty),+) => {
        impl_assign_muldiv!($($o),+);
    };
    ($t:ty) => {
        impl_assign_muldiv!(impl, $t, Vector<$t>);
        impl_assign_muldiv!(impl, $t, VV<$t>);
    };
}

impl_assign_muldiv!(Integer, Rational, Float);

macro_rules! impl_neg {
    ($t:ty) => {
        impl Neg for Vector<$t> {
            type Output = Vector<$t>;
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

impl_neg!(Integer, Rational, Float);