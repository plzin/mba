//! Vectors.
//!
//! The reason we have a [`VectorStorage`] trait and a [`Vector`] type that just
//! wraps the storage is so that we can implement traits for all vectors, i.e.
//! [`std::ops::Index`]. Otherwise we would try to write something like this
//! ```ignore
//! impl<T, V: Vector<T>> std::ops::Index<usize> for V {
//!     // ...
//! }
//! ```
//! which doesn't work because `V` need not be defined in this crate.
//! So we are essentially using the newtype pattern to be able to do this.
//!
//! Another thing to note is that the [`Vector`]s are aware of what kind of ring
//! the elements come from, i.e. the generic is the [`Ring`] and not the ring
//! element. Which means you can't have vectors of things that are not in rings.
//! The reason is that we otherwise wouldn't be able to implement the any
//! functions that depend on what the ring is, e.g. arithmetic operations.
//! Suppose the [`Vector`] didn't take the ring as its argument but the ring
//! element itself. We would try to implement (e.g.) `add_assign` like this:
//! ```ignore
//! impl<R: Ring, S: VectorStorage<R::Element>> Vector<R::Element, S> {
//!     fn add_assign<U>(&mut self, rhs: &Vector<R::Element, U>, r: &R)
//!         where U: VectorStorage<R::Element>
//!     {
//!         // ...
//!     }
//! }
//! ```
//! The problem is that the type parameter `R` is unconstrained. In this case it
//! feels like a limitation of rust because when we call `add_assign` we know
//! which ring type `R` will be because we are passing an argument of that type.
//! But if there were methods that do not take `R` as an argument in the impl
//! block, then the compiler genuinely couldn't know what `R` is.

use rand::Rng;
use std::borrow::{Borrow, BorrowMut};
use std::fmt::Debug;
use std::marker::PhantomData;
use std::ops::{Index, IndexMut, Range};

use crate::matrix::{ColumnVector, Matrix, MatrixStorage, RowVector};
use crate::rings::{Field, IntDivRing, Ring, RingElement, SqrtRing};
use crate::{CustomMetadata, CustomMetadataSlice, Half};

/// Manages how the entries of the vector are stored.
/// See the module documentation [`crate::vector`].
pub trait VectorStorage<R: Ring>: 'static {
    type Iter<'a>: DoubleEndedIterator<Item = &'a R::Element>;
    type IterMut<'a>: DoubleEndedIterator<Item = &'a mut R::Element>;

    /// The dimension of the vector.
    fn dim(&self) -> usize;

    /// Reference to an entry at a given index.
    fn entry(&self, idx: usize) -> &R::Element;

    /// Mutable reference to an entry at a given index.
    fn entry_mut(&mut self, idx: usize) -> &mut R::Element;

    /// Returns an iterator over the elements.
    fn iter(&self) -> Self::Iter<'_>;

    /// Returns an iterator over the mutable elements.
    fn iter_mut(&mut self) -> Self::IterMut<'_>;
}

/// [`VectorStorage`]s that lay out the elements contiguously in memory
/// implement this.
pub trait ContiguousVectorStorage<R: Ring>: VectorStorage<R> {
    /// Returns a slice of the vector.
    fn as_slice(&self) -> &[R::Element];

    /// Returns a mutable slice of the vector.
    fn as_slice_mut(&mut self) -> &mut [R::Element];
}

/// A vector. See the module documentation [`crate::vector`].
pub struct Vector<R: Ring, S: VectorStorage<R> + ?Sized> {
    phantom: PhantomData<R>,
    storage: S,
}

impl<R: Ring, S: VectorStorage<R> + ?Sized> Vector<R, S> {
    /// Construct the vector from its storage.
    pub fn from_storage(storage: S) -> Self
    where
        S: Sized,
    {
        Self {
            phantom: std::marker::PhantomData,
            storage,
        }
    }

    /// Creates a vector reference from a storage reference.
    pub fn from_storage_ref(storage: &S) -> &Self {
        unsafe { &*(storage as *const _ as *const _) }
    }

    /// Creates a mutable vector reference from a mutable storage reference.
    pub fn from_storage_ref_mut(storage: &mut S) -> &mut Self {
        unsafe { &mut *(storage as *mut _ as *mut _) }
    }

    /// Returns the dimension of the vector.
    pub fn dim(&self) -> usize {
        self.storage.dim()
    }

    /// Returns a reference to an entry at an index.
    pub fn entry(&self, idx: usize) -> &R::Element {
        self.storage.entry(idx)
    }

    /// Returns a mutable reference to an entry at an index.
    pub fn entry_mut(&mut self, idx: usize) -> &mut R::Element {
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

    /// Converts the vector into a row vector view, i.e.
    /// a matrix with one row.
    pub fn row_vector(&self) -> &RowVector<R, S> {
        unsafe { &*(self as *const _ as *const _) }
    }

    /// Converts the vector into a mutable row vector view, i.e.
    /// a matrix with one row.
    pub fn row_vector_mut(&mut self) -> &mut RowVector<R, S> {
        unsafe { &mut *(self as *mut _ as *mut _) }
    }

    /// Converts the vector into a column vector view, i.e.
    /// a matrix with one column.
    pub fn column_vector(&self) -> &ColumnVector<R, S> {
        unsafe { &*(self as *const _ as *const _) }
    }

    /// Converts the vector into a mutable column vector view, i.e.
    /// a matrix with one column.
    pub fn column_vector_mut(&mut self) -> &mut ColumnVector<R, S> {
        unsafe { &mut *(self as *mut _ as *mut _) }
    }

    /// Apply a function to each entry.
    pub fn transform<U, F>(&self, f: F) -> OwnedVector<U>
    where
        U: Ring,
        F: FnMut(&R::Element) -> U::Element,
    {
        Vector::from_iter(self.dim(), self.iter().map(f))
    }

    /// Apply a function to each entry.
    pub fn map_mut(&mut self, mut f: impl FnMut(&mut R::Element)) {
        for e in self.iter_mut() {
            f(e);
        }
    }

    /// Swap the two entries at index i and j.
    pub fn swap(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        unsafe {
            std::ptr::swap(self.entry_mut(i), self.entry_mut(j));
        }
    }

    /// Returns mutable references to two entries by index or panics if the
    /// indices are the same.
    pub fn get_mut_entries(&mut self, i: usize, j: usize) -> (&mut R::Element, &mut R::Element) {
        // See `Matrix::get_mut_rows`.
        assert_ne!(i, j, "Tried to mutably borrow the same entry twice.");
        let i_ptr = self.entry_mut(i) as *mut _;
        let j_ptr = self.entry_mut(j) as *mut _;
        (unsafe { &mut *i_ptr }, unsafe { &mut *j_ptr })
    }

    /// Convert the entries of this vector to a different type.
    /// This function is a bit restrictive because
    /// (e.g.) u32 does not implement From<&u32>.
    /// So what we do is clone the entries in the vector
    /// and then call into() on them.
    /// So because this clone is almost always unnecessary,
    /// we only implement this for Copy types.
    pub fn to<U: Ring>(&self) -> OwnedVector<U>
    where
        R::Element: Copy + Into<U::Element>,
    {
        #[allow(clippy::clone_on_copy)]
        self.transform(|e| e.clone().into())
    }

    /// Are all entries in the vector zero?
    pub fn is_zero(&self) -> bool {
        self.iter().all(RingElement::is_zero)
    }

    /// Negate the vector.
    pub fn neg_assign(&mut self, r: &R) {
        for e in self.iter_mut() {
            r.neg_assign(e);
        }
    }

    /// Negate the vector.
    pub fn neg(mut self, r: &R) -> Self
    where
        S: Sized,
    {
        self.neg_assign(r);
        self
    }

    /// Add a vector to this one.
    pub fn add_assign<T>(&mut self, rhs: &Vector<R, T>, r: &R)
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(
            self.dim(),
            rhs.dim(),
            "Can not add vectors of different dimensions."
        );
        for i in 0..self.dim() {
            r.add_assign(&mut self[i], &rhs[i]);
        }
    }

    /// Allocate a new vector that is the sum of two vectors.
    pub fn add_refs<T>(&self, rhs: &Vector<R, T>, r: &R) -> OwnedVector<R>
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(
            self.dim(),
            rhs.dim(),
            "Can not add vectors of incompatible dimensions."
        );
        OwnedVector::from_iter(
            self.dim(),
            self.iter()
                .zip(rhs.iter())
                .map(|(a, b)| r.add(a.clone(), b)),
        )
    }

    /// Consumes the current vector, adds another vector to it, and returns it.
    pub fn add<T>(mut self, rhs: &Vector<R, T>, r: &R) -> Self
    where
        S: Sized,
        T: VectorStorage<R> + ?Sized,
    {
        self.add_assign(rhs, r);
        self
    }

    /// Subtract a vector from this one.
    pub fn sub_assign<T>(&mut self, rhs: &Vector<R, T>, r: &R)
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(
            self.dim(),
            rhs.dim(),
            "Can not subtract vectors of different dimensions."
        );
        for i in 0..self.dim() {
            r.sub_assign(&mut self[i], &rhs[i]);
        }
    }

    /// Allocate a new vector that is the result of subtracting one vector from
    /// another.
    pub fn sub_refs<T>(&self, rhs: &Vector<R, T>, r: &R) -> OwnedVector<R>
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(
            self.dim(),
            rhs.dim(),
            "Can not subtract vectors of incompatible dimensions."
        );
        OwnedVector::from_iter(
            self.dim(),
            self.iter()
                .zip(rhs.iter())
                .map(|(a, b)| r.sub(a.clone(), b)),
        )
    }

    /// Consumes the current vector, subtracts another vector from it,
    /// and returns it.
    pub fn sub<T>(mut self, rhs: &Vector<R, T>, r: &R) -> Self
    where
        S: Sized,
        T: VectorStorage<R> + ?Sized,
    {
        self.sub_assign(rhs, r);
        self
    }

    /// Subtract a vector from this one and store the result in the given
    /// vector. See also [`Ring::sub_assign_rhs`].
    pub fn sub_assign_rhs<T>(&self, rhs: &mut Vector<R, T>, r: &R)
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(
            self.dim(),
            rhs.dim(),
            "Can not subtract vectors of incompatible dimensions."
        );
        for i in 0..self.dim() {
            r.sub_assign_rhs(&self[i], &mut rhs[i]);
        }
    }

    /// Subtract a vector given by value from this one, store the result in the
    /// argument, and return it.
    /// See [`Self::sub_assign_rhs`], [`Ring::sub_rhs`], and
    /// [`Ring::sub_assign_rhs`].
    pub fn sub_rhs<T>(&self, mut rhs: Vector<R, T>, r: &R) -> Vector<R, T>
    where
        T: VectorStorage<R>,
    {
        self.sub_assign_rhs(&mut rhs, r);
        rhs
    }

    /// Multiply the vector by a scalar.
    pub fn mul_assign(&mut self, c: &R::Element, r: &R) {
        for e in self.iter_mut() {
            r.mul_assign(e, c);
        }
    }

    /// Multiply the vector by a scalar.
    pub fn mul(mut self, c: &R::Element, r: &R) -> Self
    where
        S: Sized,
    {
        self.mul_assign(c, r);
        self
    }

    /// Allocate a new vector that is this vector scaled by the scalar.
    pub fn mul_ref(&self, c: &R::Element, r: &R) -> OwnedVector<R> {
        OwnedVector::from_iter(self.dim(), self.iter().map(|e| r.mul(e.clone(), c)))
    }

    /// Multiply a vector by a scalar and add the result to this vector.
    /// This exists as an optimization.
    pub fn mul_add_assign<T>(&mut self, c: &R::Element, v: &Vector<R, T>, r: &R)
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(self.dim(), v.dim());
        for i in 0..self.dim() {
            r.mul_add_assign(self.entry_mut(i), c, v.entry(i));
        }
    }

    /// Multiply a vector by a scalar and add the result to this vector.
    /// This exists as an optimization.
    pub fn mul_add<T>(mut self, c: &R::Element, v: &Vector<R, T>, r: &R) -> Self
    where
        S: Sized,
        T: VectorStorage<R> + ?Sized,
    {
        self.mul_add_assign(c, v, r);
        self
    }

    /// Compute the dot product of two vectors.
    pub fn dot<T>(&self, other: &Vector<R, T>, r: &R) -> R::Element
    where
        T: VectorStorage<R> + ?Sized,
    {
        assert_eq!(self.dim(), other.dim());
        self.iter()
            .zip(other.iter())
            .fold(R::zero(), |acc, (c, d)| r.mul_add(acc, c, d))
    }

    /// Computes the dot product of the vector with itself.
    pub fn norm_sqr(&self, r: &R) -> R::Element {
        self.iter().fold(R::zero(), |acc, e| r.mul_add(acc, e, e))
    }

    /// Compute square euclidean distance to another vector.
    pub fn dist_sqr<T>(&self, other: &Vector<R, T>, r: &R) -> R::Element
    where
        T: VectorStorage<R> + ?Sized,
    {
        self.iter()
            .zip(other.iter())
            .fold(R::zero(), |acc, (a, b)| {
                // TODO: We should probably try to avoid allocating a new element
                // for the diff each time. We can just reuse one.
                let diff = r.sub(a.clone(), b);
                r.mul_add(acc, &diff, &diff)
            })
    }

    /// Returns the norm of the vector.
    pub fn norm(&self, r: &R) -> R::Element
    where
        R: SqrtRing,
    {
        R::sqrt(&self.norm_sqr(r))
    }

    /// Computes the euclidean distance to another vector.
    pub fn dist<T>(&self, other: &Vector<R, T>, r: &R) -> R::Element
    where
        R: SqrtRing,
        T: VectorStorage<R> + ?Sized,
    {
        R::sqrt(&self.dist_sqr(other, r))
    }

    /// Normalize a vector, i.e. divide it by its [`Self::norm`].
    pub fn normalize(&mut self, r: &R)
    where
        R: SqrtRing + Field,
    {
        self.div_assign(&self.norm(r), r);
    }

    /// Reduce a vector using the rows a upper triangular square matrix, i.e.
    /// the first non-zero element in each row is strictly to the right of the
    /// first non-zero element in the previous row.
    pub fn reduce<T>(&mut self, m: &Matrix<R, T>, r: &R)
    where
        R: IntDivRing,
        T: MatrixStorage<R> + ?Sized,
    {
        if m.num_rows() == 0 {
            return;
        }

        let mut row = 0;
        for i in 0..self.dim() {
            let e = &m[(row, i)];
            if e.is_zero() {
                continue;
            }
            let c = r.neg(R::euclidean_div(&self[i], e));
            self.mul_add_assign(&c, m.row(row), r);
            row += 1;
            if row == m.num_rows() {
                break;
            }
        }
    }
}

impl<R: Field, S: VectorStorage<R> + ?Sized> Vector<R, S> {
    /// Divide the vector by a scalar.
    pub fn div_assign(&mut self, c: &R::Element, r: &R) {
        for e in self.iter_mut() {
            r.div_assign(e, c);
        }
    }

    /// Divide the vector by a scalar.
    pub fn div(mut self, c: &R::Element, r: &R) -> Self
    where
        S: Sized,
    {
        self.div_assign(c, r);
        self
    }

    /// Allocate a new vector that is this vector divided by the scalar.
    pub fn div_ref(&self, c: &R::Element, r: &R) -> OwnedVector<R> {
        OwnedVector::from_iter(self.dim(), self.iter().map(|e| r.div(e.clone(), c)))
    }
}

impl<R: Ring, S: ContiguousVectorStorage<R> + ?Sized> Vector<R, S> {
    /// Returns a slice of the vector.
    pub fn as_slice(&self) -> &[R::Element] {
        self.storage.as_slice()
    }

    /// Returns a mutable slice of the vector.
    pub fn as_slice_mut(&mut self) -> &mut [R::Element] {
        self.storage.as_slice_mut()
    }
}

impl<R, S> AsRef<[R::Element]> for Vector<R, S>
where
    R: Ring,
    S: ContiguousVectorStorage<R> + ?Sized,
{
    fn as_ref(&self) -> &[R::Element] {
        self.as_slice()
    }
}

impl<R, S> AsMut<[R::Element]> for Vector<R, S>
where
    R: Ring,
    S: ContiguousVectorStorage<R> + ?Sized,
{
    fn as_mut(&mut self) -> &mut [R::Element] {
        self.as_slice_mut()
    }
}

impl<R, S> Index<usize> for Vector<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    type Output = R::Element;

    fn index(&self, idx: usize) -> &Self::Output {
        self.entry(idx)
    }
}

impl<R, S> IndexMut<usize> for Vector<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    fn index_mut(&mut self, idx: usize) -> &mut Self::Output {
        self.entry_mut(idx)
    }
}

impl<R, S> Index<Range<usize>> for Vector<R, S>
where
    R: Ring,
    S: ContiguousVectorStorage<R> + ?Sized,
{
    type Output = VectorView<R>;

    fn index(&self, index: Range<usize>) -> &Self::Output {
        self.as_slice().index(index).into()
    }
}

impl<R, S> IndexMut<Range<usize>> for Vector<R, S>
where
    R: Ring,
    S: ContiguousVectorStorage<R> + ?Sized,
{
    fn index_mut(&mut self, index: Range<usize>) -> &mut Self::Output {
        self.as_slice_mut().index_mut(index).into()
    }
}

impl<'a, R, S> IntoIterator for &'a Vector<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    type Item = &'a R::Element;
    type IntoIter = S::Iter<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

impl<'a, R, S> IntoIterator for &'a mut Vector<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    type Item = &'a mut R::Element;
    type IntoIter = S::IterMut<'a>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter_mut()
    }
}

impl<R, S, I: ?Sized, U> PartialEq<I> for Vector<R, S>
where
    R: Ring,
    R::Element: PartialEq<U>,
    S: VectorStorage<R> + ?Sized,
    for<'a> &'a I: IntoIterator<Item = &'a U>,
{
    fn eq(&self, other: &I) -> bool {
        self.iter().eq(other)
    }
}

impl<R, S> std::fmt::Debug for Vector<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list().entries(self.iter()).finish()
    }
}

impl<R: Ring, S: VectorStorage<R> + ?Sized> Eq for Vector<R, S> where R::Element: Eq {}

/// Technically this is not a vector view, because it owns
/// the elements. However, you can only get references
/// to this type, which are the views.
pub type VectorView<R> = Vector<R, SliceVectorStorage<R>>;

impl<R: Ring> VectorView<R> {
    /// Create a vector view of a slice.
    pub fn from_slice(s: &[R::Element]) -> &Self {
        Vector::from_storage_ref(SliceVectorStorage::from_slice(s))
    }

    /// Create a mutable vector view of a slice.
    pub fn from_slice_mut(s: &mut [R::Element]) -> &mut Self {
        Vector::from_storage_ref_mut(SliceVectorStorage::from_slice_mut(s))
    }

    /// Create a vector view of a single element.
    pub fn from_elem(s: &R::Element) -> &Self {
        Self::from_slice(std::slice::from_ref(s))
    }

    /// Create a mutable vector view of a single element.
    pub fn from_elem_mut(s: &mut R::Element) -> &mut Self {
        Self::from_slice_mut(std::slice::from_mut(s))
    }
}

impl<'a, R: Ring> From<&'a [R::Element]> for &'a VectorView<R> {
    fn from(s: &'a [R::Element]) -> Self {
        VectorView::from_slice(s)
    }
}

impl<'a, R: Ring> From<&'a mut [R::Element]> for &'a mut VectorView<R> {
    fn from(s: &'a mut [R::Element]) -> Self {
        VectorView::from_slice_mut(s)
    }
}

impl<R: Ring> ToOwned for VectorView<R> {
    type Owned = OwnedVector<R>;

    fn to_owned(&self) -> Self::Owned {
        OwnedVector::from_iter(self.dim(), self.iter().cloned())
    }
}

#[derive(Debug)]
#[repr(transparent)]
pub struct SliceVectorStorage<R: Ring>([R::Element]);

impl<R: Ring> SliceVectorStorage<R> {
    pub fn from_slice(s: &[R::Element]) -> &Self {
        unsafe { &*std::ptr::from_raw_parts(s.as_ptr() as _, s.len()) }
    }

    pub fn from_slice_mut(s: &mut [R::Element]) -> &mut Self {
        unsafe { &mut *std::ptr::from_raw_parts_mut(s.as_mut_ptr() as _, s.len()) }
    }
}

impl<R: Ring> VectorStorage<R> for SliceVectorStorage<R> {
    type Iter<'a>
        = std::slice::Iter<'a, R::Element>
    where
        R: 'a;
    type IterMut<'a>
        = std::slice::IterMut<'a, R::Element>
    where
        R: 'a;

    fn dim(&self) -> usize {
        self.0.len()
    }

    fn entry(&self, idx: usize) -> &R::Element {
        &self.0[idx]
    }

    fn entry_mut(&mut self, idx: usize) -> &mut R::Element {
        &mut self.0[idx]
    }

    fn iter(&self) -> Self::Iter<'_> {
        self.0.iter()
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        self.0.iter_mut()
    }
}

impl<R: Ring> ContiguousVectorStorage<R> for SliceVectorStorage<R> {
    fn as_slice(&self) -> &[R::Element] {
        &self.0
    }

    fn as_slice_mut(&mut self) -> &mut [R::Element] {
        &mut self.0
    }
}

pub type OwnedVector<R> = Vector<R, OwnedVectorStorage<R>>;

impl<R: Ring> OwnedVector<R> {
    /// Returns a view of the owned vector.
    /// Currently, this has the same layout as the owned vector,
    /// so this is pretty much a no-op.
    pub fn view(&self) -> &VectorView<R> {
        VectorView::from_slice(self.as_slice())
    }

    /// Returns a mutable view of the owned vector.
    /// Currently, this has the same layout as the owned vector,
    /// so this is pretty much a no-op.
    pub fn view_mut(&mut self) -> &mut VectorView<R> {
        VectorView::from_slice_mut(self.as_slice_mut())
    }

    /// Returns an empty vector.
    pub fn empty() -> Self {
        Self::from_storage(OwnedVectorStorage {
            entries: Vec::new(),
        })
    }

    /// Returns a zero vector.
    pub fn zero(dim: usize) -> Self {
        Self::from_iter(dim, core::iter::repeat_with(R::zero))
    }

    /// Creates a vector from a slice.
    pub fn from_entries<U, V>(a: U) -> Self
    where
        U: AsRef<[V]>,
        V: Into<R::Element> + Clone,
    {
        let a = a.as_ref();
        Self::from_iter(a.len(), a.iter().cloned().map(Into::into))
    }

    /// Creates a vector from an iterator.
    pub fn from_iter<U: Iterator<Item = R::Element>>(dim: usize, mut iter: U) -> Self {
        let mut entries = Vec::with_capacity(dim);
        for _ in 0..dim {
            let e = iter.next().expect("Iter needs to return `dim` elements.");
            entries.push(e);
        }

        Self::from_storage(OwnedVectorStorage { entries })
    }

    /// Try to create a vector from an iterator.
    pub fn try_from_iter<E, U: Iterator<Item = Result<R::Element, E>>>(
        dim: usize,
        mut iter: U,
    ) -> Result<Self, E> {
        let mut entries = Vec::with_capacity(dim);
        for _ in 0..dim {
            let e = iter.next().expect("Iter needs to return `dim` elements.")?;
            entries.push(e);
        }

        Ok(Self::from_storage(OwnedVectorStorage { entries }))
    }

    /// Creates a vector from an array.
    pub fn from_array<U: Into<R::Element>, const D: usize>(a: [U; D]) -> Self {
        Self::from_iter(D, a.into_iter().map(U::into))
    }

    /// Creates a random vector.
    pub fn random<Rand: Rng>(dim: usize, ring: &R, rng: &mut Rand) -> Self {
        Self::from_iter(dim, std::iter::repeat_with(|| ring.random(rng)))
    }

    /// Applies a function to each entry.
    pub fn map<F: FnMut(R::Element) -> R::Element>(self, mut f: F) -> OwnedVector<R> {
        let mut entries = Vec::with_capacity(self.dim());
        for e in self.storage.entries {
            entries.push(f(e));
        }

        Self::from_storage(OwnedVectorStorage { entries })
    }

    /// Appends an element to the end.
    pub fn append(&mut self, e: R::Element) {
        self.storage.append(e)
    }

    /// Returns an owned vector from a pointer and dimension.
    pub fn from_raw_entries(entries: Vec<R::Element>) -> Self {
        Self::from_storage(OwnedVectorStorage { entries })
    }
}

impl<R: Ring> Borrow<VectorView<R>> for OwnedVector<R> {
    fn borrow(&self) -> &VectorView<R> {
        self.view()
    }
}

impl<R: Ring> BorrowMut<VectorView<R>> for OwnedVector<R> {
    fn borrow_mut(&mut self) -> &mut VectorView<R> {
        self.view_mut()
    }
}

impl<R: Ring> Clone for OwnedVector<R> {
    fn clone(&self) -> Self {
        Self::from_iter(self.dim(), self.iter().cloned())
    }
}

pub struct OwnedVectorStorage<R: Ring> {
    /// Memory that holds the entries.
    pub entries: Vec<R::Element>,
}

impl<R: Ring> OwnedVectorStorage<R> {
    /// Appends an element to the end.
    pub fn append(&mut self, e: R::Element) {
        self.entries.push(e);
    }
}

impl<R: Ring> VectorStorage<R> for OwnedVectorStorage<R> {
    type Iter<'a>
        = std::slice::Iter<'a, R::Element>
    where
        R: 'a;
    type IterMut<'a>
        = std::slice::IterMut<'a, R::Element>
    where
        R: 'a;

    fn dim(&self) -> usize {
        self.entries.len()
    }

    fn entry(&self, idx: usize) -> &R::Element {
        &self.entries[idx]
    }

    fn entry_mut(&mut self, idx: usize) -> &mut R::Element {
        &mut self.entries[idx]
    }

    fn iter(&self) -> Self::Iter<'_> {
        self.as_slice().iter()
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        self.as_slice_mut().iter_mut()
    }
}

impl<R: Ring> ContiguousVectorStorage<R> for OwnedVectorStorage<R> {
    fn as_slice(&self) -> &[R::Element] {
        self.entries.as_slice()
    }

    fn as_slice_mut(&mut self) -> &mut [R::Element] {
        self.entries.as_mut_slice()
    }
}

/// A vector view where elements are skipped.
/// The stride is the amount of elements to skip
/// not the number of bytes to skip.
pub type StrideVectorView<R> = Vector<R, StrideStorage<R>>;

impl<R: Ring> StrideVectorView<R> {
    /// Returns a vector view with stride from a raw pointer and dimension.
    unsafe fn from_raw_parts<'a>(ptr: *const R::Element, dim: usize, stride: usize) -> &'a Self {
        Self::from_storage_ref(unsafe { StrideStorage::from_raw_parts(ptr, dim, stride) })
    }

    /// Returns a mutable vector view with stride from a raw pointer and dimension.
    unsafe fn from_raw_parts_mut<'a>(
        ptr: *mut R::Element,
        dim: usize,
        stride: usize,
    ) -> &'a mut Self {
        Self::from_storage_ref_mut(unsafe { StrideStorage::from_raw_parts_mut(ptr, dim, stride) })
    }

    /// Returns a vector view with stride from a slice.
    /// The vector view will have the maximum possible dimension
    /// that fits into the slice.
    pub fn from_stride_slice(s: &[R::Element], stride: usize) -> &Self {
        assert!(stride > 0);
        let dim = s.len().div_ceil(stride);
        unsafe { Self::from_raw_parts(s.as_ptr(), dim, stride) }
    }

    /// Returns a vector view with stride from a slice.
    /// The vector view will have the maximum possible dimension
    /// that fits into the slice.
    pub fn from_stride_slice_mut(s: &mut [R::Element], stride: usize) -> &mut Self {
        assert!(stride > 0);
        let dim = s.len() + (stride - 1) / stride;
        unsafe { Self::from_raw_parts_mut(s.as_mut_ptr(), dim, stride) }
    }
}

/// This is extremely hacky, because rust currently does not support
/// custom dynamically sized types (DST). We make this a DST by having
/// a slice as the member, just as for [ContiguousVectorStorage].
/// We abuse the metadata of the slice to store the dimension and stride!
/// (See [std::ptr::Pointee], for what metadata is). Since we cannot
/// implement `Pointee` ourselves, we have to use the usize metadata of the
/// slice. The first 32 bits are the dimension, the second 32 bits are the
/// stride. As I said, very hacky, but in practice the both the dimension and
/// stride should fit into 32 bits.
#[repr(transparent)]
pub struct StrideStorage<R: Ring>(CustomMetadataSlice<R::Element, StrideStorageMetadata>);

struct StrideStorageMetadata {
    dim: Half,
    stride: Half,
}

impl StrideStorageMetadata {
    pub fn new(dim: usize, stride: usize) -> Self {
        Self {
            dim: dim.try_into().unwrap(),
            stride: stride.try_into().unwrap(),
        }
    }
}

impl CustomMetadata for StrideStorageMetadata {
    fn size(&self) -> usize {
        self.dim as usize * self.stride as usize
    }
}

impl<R: Ring> StrideStorage<R> {
    pub(crate) unsafe fn from_raw_parts<'a>(
        data: *const R::Element,
        dim: usize,
        stride: usize,
    ) -> &'a Self {
        let metadata = StrideStorageMetadata::new(dim, stride);
        unsafe { std::mem::transmute(CustomMetadataSlice::new(data, metadata)) }
    }

    pub(crate) unsafe fn from_raw_parts_mut<'a>(
        data: *mut R::Element,
        dim: usize,
        stride: usize,
    ) -> &'a mut Self {
        let metadata = StrideStorageMetadata::new(dim, stride);
        unsafe { std::mem::transmute(CustomMetadataSlice::new_mut(data, metadata)) }
    }

    /// The stride of the vector.
    pub fn stride(&self) -> usize {
        self.0.metadata().stride as usize
    }
}

impl<R: Ring> VectorStorage<R> for StrideStorage<R> {
    type Iter<'a>
        = StrideStorageIter<'a, R::Element>
    where
        R: 'a;
    type IterMut<'a>
        = StrideStorageIterMut<'a, R::Element>
    where
        R: 'a;

    fn dim(&self) -> usize {
        self.0.metadata().dim as usize
    }

    fn entry(&self, idx: usize) -> &R::Element {
        self.0.slice().index(idx * self.stride())
    }

    fn entry_mut(&mut self, idx: usize) -> &mut R::Element {
        let stride = self.stride();
        self.0.slice_mut().index_mut(idx * stride)
    }

    fn iter(&self) -> Self::Iter<'_> {
        let ptr = self.0.slice().as_ptr();
        StrideStorageIter {
            start: ptr,
            end: unsafe { ptr.add(self.dim() * self.stride()) },
            stride: self.stride(),
            marker: PhantomData,
        }
    }

    fn iter_mut(&mut self) -> Self::IterMut<'_> {
        let ptr = self.0.slice_mut().as_mut_ptr();
        StrideStorageIterMut {
            start: ptr,
            end: unsafe { ptr.add(self.dim() * self.stride()) },
            stride: self.stride(),
            marker: PhantomData,
        }
    }
}

pub struct StrideStorageIter<'a, T> {
    pub(self) start: *const T,
    pub(self) end: *const T,
    pub(self) stride: usize,
    pub(self) marker: PhantomData<&'a T>,
}

impl<'a, T> Iterator for StrideStorageIter<'a, T> {
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.end {
            None
        } else {
            let ptr = self.start;
            self.start = unsafe { ptr.add(self.stride) };
            Some(unsafe { &*ptr })
        }
    }
}

impl<'a, T> DoubleEndedIterator for StrideStorageIter<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start >= self.end {
            None
        } else {
            self.end = unsafe { self.end.sub(self.stride) };
            Some(unsafe { &*self.end })
        }
    }
}

pub struct StrideStorageIterMut<'a, T> {
    pub(self) start: *mut T,
    pub(self) end: *mut T,
    pub(self) stride: usize,
    pub(self) marker: PhantomData<&'a mut T>,
}

impl<'a, T> Iterator for StrideStorageIterMut<'a, T> {
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start >= self.end {
            None
        } else {
            let ptr = self.start;
            self.start = unsafe { ptr.add(self.stride) };
            Some(unsafe { &mut *ptr })
        }
    }
}

impl<'a, T> DoubleEndedIterator for StrideStorageIterMut<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start >= self.end {
            None
        } else {
            self.end = unsafe { self.end.sub(self.stride) };
            Some(unsafe { &mut *self.end })
        }
    }
}

#[test]
fn reduce_test() {
    use crate::rings::U8;
    let mut v = Vector::from_entries([9, 9, 9, 9]);
    let m = Matrix::from_rows(&[&[2, 2, 0, 0], &[0, 0, 3, 2], &[0, 0, 0, 1]]);
    v.reduce(&m, &U8);
    println!("{m:3?}");
    assert_eq!(v, [1, 1, 0, 0]);
}
