//! Matrices.
//! See the module documentation for [`crate::vector`].
//! The same thing applies here.

use std::{
    fmt::Debug,
    marker::PhantomData,
    ops::{Index, IndexMut},
};

use itertools::iproduct;
use rand::Rng;

use crate::{
    CustomMetadata, CustomMetadataSlice, Half, rings::Ring, vector::*,
};

/// How are the entries of a matrix stored?
pub trait MatrixStorage<R: Ring>: 'static {
    /// Row/Column vector type.
    type RowVecStorage: VectorStorage<R> + ?Sized;
    type ColVecStorage: VectorStorage<R> + ?Sized;

    /// Row/column iter type.
    type RowIter<'a>: DoubleEndedIterator<
        Item = &'a Vector<R, Self::RowVecStorage>,
    >;
    type RowIterMut<'a>: DoubleEndedIterator<
        Item = &'a mut Vector<R, Self::RowVecStorage>,
    >;
    type ColIter<'a>: DoubleEndedIterator<
        Item = &'a Vector<R, Self::ColVecStorage>,
    >;
    type ColIterMut<'a>: DoubleEndedIterator<
        Item = &'a mut Vector<R, Self::ColVecStorage>,
    >;

    /// Returns the number of rows.
    fn rows(&self) -> usize;

    /// Returns the number of columns.
    fn cols(&self) -> usize;

    /// Returns a view of the row `r`.
    fn row(&self, r: usize) -> &Self::RowVecStorage;

    /// Returns a mutable view of the row `r`.
    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage;

    /// Returns an iterator over the rows.
    fn row_iter<'a>(&'a self) -> Self::RowIter<'a>;

    /// Returns an iterator over the rows.
    fn row_iter_mut<'a>(&'a mut self) -> Self::RowIterMut<'a>;

    /// Returns a view of the column `c`.
    fn col(&self, c: usize) -> &Self::ColVecStorage;

    /// Returns a mutable view of the column `c`.
    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage;

    /// Returns an iterator over the columns.
    fn col_iter<'a>(&'a self) -> Self::ColIter<'a>;

    /// Returns an iterator over the columns.
    fn col_iter_mut<'a>(&'a mut self) -> Self::ColIterMut<'a>;

    /// Returns a reference to the entry at row `r` and column `c`.
    fn entry(&self, r: usize, c: usize) -> &R::Element {
        self.row(r).entry(c)
    }

    /// Returns a mutable reference to the entry at row `r` and column `c`.
    fn entry_mut(&mut self, r: usize, c: usize) -> &mut R::Element {
        self.row_mut(r).entry_mut(c)
    }
}

/// A matrix.
pub struct Matrix<R: Ring, S: MatrixStorage<R> + ?Sized> {
    marker: PhantomData<R>,
    storage: S,
}

impl<R: Ring, S: MatrixStorage<R> + ?Sized> Matrix<R, S> {
    fn from_storage(s: S) -> Self
    where
        S: Sized,
    {
        Self { marker: PhantomData, storage: s }
    }

    /// The number of rows of the matrix.
    pub fn num_rows(&self) -> usize {
        self.storage.rows()
    }

    /// The number of columns of the matrix.
    pub fn num_cols(&self) -> usize {
        self.storage.cols()
    }

    /// Returns the smaller of the two dimensions.
    pub fn min_dim(&self) -> usize {
        std::cmp::min(self.num_rows(), self.num_cols())
    }

    /// Is the matrix empty, i.e. has it zero rows or columns?
    pub fn is_empty(&self) -> bool {
        self.num_rows() == 0 || self.num_cols() == 0
    }

    /// Returns a view of the row `r`.
    pub fn row(&self, r: usize) -> &Vector<R, S::RowVecStorage> {
        Vector::from_storage_ref(self.storage.row(r))
    }

    /// Returns a mutable view of the row `r`.
    pub fn row_mut(&mut self, r: usize) -> &mut Vector<R, S::RowVecStorage> {
        Vector::from_storage_ref_mut(self.storage.row_mut(r))
    }

    /// Returns a view of the column `c`.
    pub fn col(&self, c: usize) -> &Vector<R, S::ColVecStorage> {
        Vector::from_storage_ref(self.storage.col(c))
    }

    /// Returns a mutable view of the column `c`.
    pub fn col_mut(&mut self, c: usize) -> &mut Vector<R, S::ColVecStorage> {
        Vector::from_storage_ref_mut(self.storage.col_mut(c))
    }

    /// Returns a reference to the entry at row `r` and column `c`.
    pub fn entry(&self, r: usize, c: usize) -> &R::Element {
        self.storage.entry(r, c)
    }

    /// Returns a mutable reference to the entry at row `r` and column `c`.
    pub fn entry_mut(&mut self, r: usize, c: usize) -> &mut R::Element {
        self.storage.entry_mut(r, c)
    }

    /// Returns an iterator over the rows.
    pub fn rows(&self) -> S::RowIter<'_> {
        self.storage.row_iter()
    }

    /// Returns an iterator over the mutable rows.
    pub fn rows_mut(&mut self) -> S::RowIterMut<'_> {
        self.storage.row_iter_mut()
    }

    /// Returns an iterator over the columns.
    pub fn cols(&self) -> S::ColIter<'_> {
        self.storage.col_iter()
    }

    /// Returns an iterator over the mutable columns.
    pub fn cols_mut(&mut self) -> S::ColIterMut<'_> {
        self.storage.col_iter_mut()
    }

    /// Returns an iterator over the entries in row-major order.
    pub fn entries_row_major(&self) -> impl Iterator<Item = &R::Element> {
        self.rows().flat_map(|r| r.iter())
    }

    /// Returns an iterator over the mutable entries in row-major order.
    pub fn entries_row_major_mut(
        &mut self,
    ) -> impl Iterator<Item = &mut R::Element> {
        self.rows_mut().flat_map(|r| r.iter_mut())
    }

    /// Returns an iterator over the entries in column-major order.
    pub fn entries_col_major(&self) -> impl Iterator<Item = &R::Element> {
        self.cols().flat_map(|c| c.iter())
    }

    /// Returns an iterator over the mutable entries in column-major order.
    pub fn entries_col_major_mut(
        &mut self,
    ) -> impl Iterator<Item = &mut R::Element> {
        self.cols_mut().flat_map(|c| c.iter_mut())
    }

    /// Call a function on each entry and return the resulting matrix.
    pub fn transform<U: Ring, F: FnMut(&R::Element) -> U::Element>(
        &self,
        f: F,
    ) -> OwnedMatrix<U> {
        OwnedMatrix::from_iter(
            self.num_rows(),
            self.num_cols(),
            self.entries_row_major().map(f),
        )
    }

    /// Call a function on each entry.
    pub fn map_mut<F: FnMut(&mut R::Element)>(&mut self, f: F) {
        self.entries_row_major_mut().for_each(f)
    }

    /// Swap two rows.
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        unsafe {
            core::ptr::swap_nonoverlapping(
                self.entry_mut(i, 0),
                self.entry_mut(j, 0),
                self.num_cols(),
            )
        }
    }

    /// Swap two columns.
    pub fn swap_columns(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        for k in 0..self.num_rows() {
            unsafe {
                core::ptr::swap(self.entry_mut(k, i), self.entry_mut(k, j));
            }
        }
    }

    /// Returns a mutable references to two rows of the matrix. Panics if the
    /// indices are the same.
    #[allow(clippy::type_complexity)]
    pub fn get_rows_mut(
        &mut self,
        i: usize,
        j: usize,
    ) -> (&mut Vector<R, S::RowVecStorage>, &mut Vector<R, S::RowVecStorage>)
    {
        // I actually think this is not UB, although it is a bit ugly. It could
        // maybe be less ugly if the `MatrixStorage` returned pointers to the
        // row storages instead of references, but I don't think it matters.
        assert_ne!(i, j, "Tried to mutably borrow the same row twice.");
        let i_ptr = self.row_mut(i) as *mut _;
        let j_ptr = self.row_mut(j) as *mut _;
        (unsafe { &mut *i_ptr }, unsafe { &mut *j_ptr })
    }

    /// Returns a copy of the matrix.
    pub fn to_owned(&self) -> OwnedMatrix<R> {
        OwnedMatrix::from_iter(
            self.num_rows(),
            self.num_cols(),
            self.entries_row_major().cloned(),
        )
    }

    /// Creates an owned matrix that is the transpose of the current matrix.
    pub fn transposed(&self) -> OwnedMatrix<R> {
        OwnedMatrix::from_iter(
            self.num_cols(),
            self.num_rows(),
            self.entries_col_major().cloned(),
        )
    }

    /// Negates all elements of a row.
    pub fn negate_row(&mut self, row: usize, r: &R) {
        for e in self.row_mut(row) {
            r.neg_assign(e);
        }
    }

    /// Negates all elements of a column.
    pub fn negate_col(&mut self, col: usize, r: &R) {
        for e in self.col_mut(col) {
            r.neg_assign(e);
        }
    }

    /// Multiply a row by an element.
    /// Note that this is only an invertible operation if `c` is a unit.
    pub fn row_multiply(&mut self, row: usize, c: &R::Element, r: &R) {
        for e in self.row_mut(row) {
            r.mul_assign(e, c);
        }
    }

    /// Multiply a column by an element.
    /// Note that this is only an invertible operation if `c` is a unit.
    pub fn col_multiply(&mut self, col: usize, c: &R::Element, r: &R) {
        for e in self.col_mut(col) {
            r.mul_assign(e, c);
        }
    }

    /// Add a scaled row to another row. N += M * c.
    /// `m` and `n` can not be equal.
    pub fn row_multiply_add(
        &mut self,
        n: usize,
        m: usize,
        c: &R::Element,
        r: &R,
    ) {
        assert_ne!(m, n);
        assert!(n < self.num_rows() && m < self.num_rows());
        for i in 0..self.num_cols() {
            // We need unsafe to assert that these are not the same element,
            // which they are not because `m != n`.
            let ne_ptr = self.entry_mut(n, i) as *mut _;
            let me_ptr = self.entry(m, i) as *const _;
            r.mul_add_assign(unsafe { &mut *ne_ptr }, unsafe { &*me_ptr }, c);
        }
    }

    /// Add a scaled column to another column. N += M * c.
    /// `m` and `n` can not be equal.
    pub fn col_multiply_add(
        &mut self,
        n: usize,
        m: usize,
        c: &R::Element,
        r: &R,
    ) {
        assert_ne!(m, n);
        assert!(n < self.num_cols() && m < self.num_cols());
        for i in 0..self.num_rows() {
            // We need unsafe to assert that these are not the same element,
            // which they are not because `m != n`.
            let ne_ptr = self.entry_mut(i, n) as *mut _;
            let me_ptr = self.entry(i, m) as *const _;
            r.mul_add_assign(unsafe { &mut *ne_ptr }, unsafe { &*me_ptr }, c);
        }
    }

    /// Multiply two matrices.
    pub fn mul<T>(&self, rhs: &Matrix<R, T>, ring: &R) -> OwnedMatrix<R>
    where
        T: MatrixStorage<R> + ?Sized,
    {
        let r = self.num_rows();
        let c = rhs.num_cols();
        let iter = iproduct!(0..r, 0..c)
            .map(|(r, c)| self.row(r).dot(rhs.col(c), ring));
        OwnedMatrix::from_iter(r, c, iter)
    }

    /// Post-multiply a matrix with a vector.
    pub fn mul_vec_post<T>(&self, rhs: &Vector<R, T>, r: &R) -> OwnedVector<R>
    where
        T: VectorStorage<R> + ?Sized,
    {
        self.mul(rhs.column_vector(), r).into_vector()
    }

    /// Pre-multiply a matrix with a vector.
    pub fn mul_vec_pre<T>(&self, lhs: &Vector<R, T>, r: &R) -> OwnedVector<R>
    where
        T: VectorStorage<R> + ?Sized,
    {
        lhs.row_vector().mul(self, r).into_vector()
    }

    /// Prints each row of the matrix in its own line using the debug formatter.
    pub fn print_rows(&self) {
        for r in self.rows() {
            println!("{r:?}");
        }
    }
}

impl<R, S> Index<usize> for Matrix<R, S>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
{
    type Output = Vector<R, S::RowVecStorage>;

    fn index(&self, index: usize) -> &Self::Output {
        self.row(index)
    }
}

impl<R, S> IndexMut<usize> for Matrix<R, S>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
{
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.row_mut(index)
    }
}

impl<R, S> Index<(usize, usize)> for Matrix<R, S>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
{
    type Output = R::Element;

    fn index(&self, (r, c): (usize, usize)) -> &Self::Output {
        self.entry(r, c)
    }
}

impl<R, S> IndexMut<(usize, usize)> for Matrix<R, S>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
{
    fn index_mut(&mut self, (r, c): (usize, usize)) -> &mut Self::Output {
        self.entry_mut(r, c)
    }
}

impl<R, S, T> PartialEq<Matrix<R, S>> for Matrix<R, T>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
    T: MatrixStorage<R> + ?Sized,
{
    fn eq(&self, other: &Matrix<R, S>) -> bool {
        self.num_rows() == other.num_rows()
            && self.num_cols() == other.num_cols()
            && self.entries_row_major().eq(other.entries_row_major())
    }
}

impl<R, S> Debug for Matrix<R, S>
where
    R: Ring,
    S: MatrixStorage<R> + ?Sized,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        //write!(f, "({}×{})", self.nrows(), self.ncols())?;
        f.debug_list().entries(self.rows()).finish()
    }
}

pub type OwnedMatrix<R> = Matrix<R, OwnedMatrixStorage<R>>;

impl<R: Ring> OwnedMatrix<R> {
    /// Return an empty matrix.
    pub fn empty() -> Self {
        Self::from_storage(OwnedMatrixStorage::empty())
    }

    /// Returns the number of entries.
    ///
    /// This function exists because `OwnedMatrix` stores the number of entries
    /// directly, so it doesn't have to be computed.
    pub fn num_entries(&self) -> usize {
        self.storage.entries.len()
    }

    /// Returns a slice containing all elements.
    pub fn as_slice(&self) -> &[R::Element] {
        self.storage.as_slice()
    }

    /// Returns a slice containing all elements.
    pub fn as_slice_mut(&mut self) -> &mut [R::Element] {
        self.storage.as_slice_mut()
    }

    /// Returns a view of this matrix.
    pub fn view(&self) -> &MatrixView<R> {
        unsafe {
            MatrixView::from_raw_parts(
                self.storage.entries.as_ptr(),
                self.num_rows(),
                self.num_cols(),
            )
        }
    }

    /// Returns a mutable view of this matrix.
    pub fn view_mut(&mut self) -> &mut MatrixView<R> {
        unsafe {
            MatrixView::from_raw_parts_mut(
                self.storage.entries.as_mut_ptr(),
                self.num_rows(),
                self.num_cols(),
            )
        }
    }

    /// Turns this matrix into a vector of dimension r*c.
    pub fn into_vector(self) -> OwnedVector<R> {
        // Transfer the memory to the vector.
        Vector::from_raw_entries(self.storage.entries)
    }

    /// Creates a matrix from an array of rows.
    pub fn from_array<U: Into<R::Element>, const RS: usize, const CS: usize>(
        a: [[U; CS]; RS],
    ) -> Self {
        let mut entries = Vec::with_capacity(RS * CS);
        for row in a {
            for e in row {
                entries.push(e.into());
            }
        }

        Self::from_storage(OwnedMatrixStorage { entries, rows: RS, cols: CS })
    }

    /// Creates a matrix from an iterator in row-major order.
    pub fn from_iter<I: Iterator<Item = R::Element>>(
        r: usize,
        c: usize,
        mut iter: I,
    ) -> Self {
        let num_entries = r * c;
        let mut entries = Vec::with_capacity(num_entries);
        for _ in 0..num_entries {
            let e = iter
                .next()
                .expect("The iterator needs to return at least r * c items.");
            entries.push(e);
        }

        Self::from_storage(OwnedMatrixStorage { entries, rows: r, cols: c })
    }

    /// Tries to create a matrix from an iterator in row-major order.
    pub fn try_from_iter<E, I: Iterator<Item = Result<R::Element, E>>>(
        r: usize,
        c: usize,
        mut iter: I,
    ) -> Result<Self, E> {
        let num_entries = r * c;
        let mut entries = Vec::with_capacity(num_entries);
        for _ in 0..num_entries {
            let e = iter
                .next()
                .expect("The iterator needs to return at least r * c items.")?;
            entries.push(e);
        }

        Ok(Self::from_storage(OwnedMatrixStorage { entries, rows: r, cols: c }))
    }

    /// Creates a matrix from slice of rows.
    pub fn from_rows<U, V>(rows: &[U]) -> Self
    where
        U: AsRef<[V]>,
        V: Into<R::Element> + Clone,
    {
        if rows.is_empty() {
            return Self::empty();
        }

        let r = rows.len();
        let c = rows[0].as_ref().len();
        assert!(rows.iter().all(|r| r.as_ref().len() == c));

        Matrix::from_iter(
            r,
            c,
            rows.iter()
                .flat_map(|r| r.as_ref().iter().map(|e| e.clone().into())),
        )
    }

    /// Creates a random matrix.
    pub fn random<Rand: Rng>(
        r: usize,
        c: usize,
        ring: &R,
        rng: &mut Rand,
    ) -> Self {
        Self::from_iter(r, c, std::iter::repeat_with(|| ring.random(rng)))
    }

    /// View of the transpose.
    pub fn transpose(&self) -> &TransposedMatrixView<R> {
        self.view().transpose()
    }

    /// Mutable view of the transpose.
    pub fn transpose_mut(&mut self) -> &mut TransposedMatrixView<R> {
        self.view_mut().transpose_mut()
    }

    /// Removes rows at the end of the matrix.
    pub fn shrink(&mut self, new_nrows: usize) {
        self.storage.shrink(new_nrows);
    }

    /// Returns an r×c zero matrix.
    pub fn zero(r: usize, c: usize) -> Self {
        Self::from_storage(OwnedMatrixStorage::zero(r, c))
    }

    /// Returns an nxn identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zero(n, n);
        for i in 0..n {
            m[(i, i)] = R::one();
        }
        m
    }

    /// Remove rows of zeros at the end of the matrix.
    pub fn remove_zero_rows(&mut self) {
        let num = self.rows().rev().take_while(|r| r.is_zero()).count();

        self.shrink(self.num_rows() - num);
    }

    /// Append `n` zero rows to the end of the matrix.
    pub fn append_zero_rows(&mut self, n: usize) {
        self.storage.append_zero_rows(n);
    }
}

impl<R: Ring> Clone for OwnedMatrix<R> {
    fn clone(&self) -> Self {
        OwnedMatrix::from_iter(
            self.num_rows(),
            self.num_cols(),
            self.entries_row_major().cloned(),
        )
    }
}

pub struct OwnedMatrixStorage<R: Ring> {
    /// Memory that holds the entries.
    pub(self) entries: Vec<R::Element>,

    /// The number of rows.
    pub(self) rows: usize,

    /// The number of columns.
    pub(self) cols: usize,
}

impl<R: Ring> OwnedMatrixStorage<R> {
    /// Returns a 0×0 matrix.
    pub fn empty() -> Self {
        Self { entries: Vec::new(), rows: 0, cols: 0 }
    }

    /// Returns a zero-initialized r×c matrix.
    pub fn zero(r: usize, c: usize) -> Self {
        Self { entries: vec![R::zero(); r * c], rows: r, cols: c }
    }

    /// Returns a slice containing all elements.
    pub fn as_slice(&self) -> &[R::Element] {
        self.entries.as_slice()
    }

    /// Returns a slice containing all elements.
    pub fn as_slice_mut(&mut self) -> &mut [R::Element] {
        self.entries.as_mut_slice()
    }

    /// Removes rows at the end of the matrix.
    pub fn shrink(&mut self, new_nrows: usize) {
        if new_nrows == self.rows {
            return;
        }
        assert!(new_nrows < self.rows);

        self.entries.truncate(new_nrows * self.cols);
        self.rows = new_nrows;
    }

    /// Append zero rows to the end of the matrix.
    pub fn append_zero_rows(&mut self, n: usize) {
        if n == 0 {
            return;
        }

        if self.cols == 0 {
            self.rows += n;
            return;
        }

        self.rows += n;
        self.entries.resize(self.rows * self.cols, R::zero());
    }
}

impl<R: Ring> MatrixStorage<R> for OwnedMatrixStorage<R> {
    type ColIter<'a> = MatrixStridedIter<'a, R>;
    type ColIterMut<'a> = MatrixStridedIterMut<'a, R>;
    type ColVecStorage = StrideStorage<R>;
    type RowIter<'a> = MatrixSliceIter<'a, R>;
    type RowIterMut<'a> = MatrixSliceIterMut<'a, R>;
    type RowVecStorage = SliceVectorStorage<R>;

    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn row(&self, r: usize) -> &Self::RowVecStorage {
        assert!(r < self.rows);
        let start = r * self.cols;
        SliceVectorStorage::from_slice(&self.entries[start..start + self.cols])
    }

    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage {
        assert!(r < self.rows);
        let start = r * self.cols;
        SliceVectorStorage::from_slice_mut(
            &mut self.entries[start..start + self.cols],
        )
    }

    fn row_iter(&self) -> Self::RowIter<'_> {
        unsafe {
            MatrixSliceIter::new(self.entries.as_ptr(), self.cols, self.rows)
        }
    }

    fn row_iter_mut(&mut self) -> Self::RowIterMut<'_> {
        unsafe {
            MatrixSliceIterMut::new(
                self.entries.as_mut_ptr(),
                self.cols,
                self.rows,
            )
        }
    }

    fn col(&self, c: usize) -> &Self::ColVecStorage {
        assert!(c < self.cols);
        unsafe {
            StrideStorage::from_raw_parts(
                &self.entries[c],
                self.rows,
                self.cols,
            )
        }
    }

    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage {
        assert!(c < self.cols);
        unsafe {
            StrideStorage::from_raw_parts_mut(
                &mut self.entries[c],
                self.rows,
                self.cols,
            )
        }
    }

    fn col_iter(&self) -> Self::ColIter<'_> {
        unsafe {
            MatrixStridedIter::new(self.entries.as_ptr(), self.rows, self.cols)
        }
    }

    fn col_iter_mut(&mut self) -> Self::ColIterMut<'_> {
        unsafe {
            MatrixStridedIterMut::new(
                self.entries.as_mut_ptr(),
                self.rows,
                self.cols,
            )
        }
    }
}

pub type MatrixView<R> = Matrix<R, SliceMatrixStorage<R>>;

impl<R: Ring> MatrixView<R> {
    /// Returns a view of the matrix transposed.
    pub fn transpose(&self) -> &TransposedMatrixView<R> {
        unsafe { &*(self as *const _ as *const _) }
    }

    /// Returns a mutable view of the matrix transposed.
    pub fn transpose_mut(&mut self) -> &mut TransposedMatrixView<R> {
        unsafe { &mut *(self as *mut _ as *mut _) }
    }

    unsafe fn from_raw_parts<'a>(
        entries: *const R::Element,
        rows: usize,
        cols: usize,
    ) -> &'a Self {
        let metadata = SliceMatrixStorageMetadata::new(rows, cols);
        unsafe {
            std::mem::transmute(CustomMetadataSlice::new(entries, metadata))
        }
    }

    unsafe fn from_raw_parts_mut<'a>(
        entries: *mut R::Element,
        rows: usize,
        cols: usize,
    ) -> &'a mut Self {
        let metadata = SliceMatrixStorageMetadata::new(rows, cols);
        unsafe {
            std::mem::transmute(CustomMetadataSlice::new_mut(entries, metadata))
        }
    }
}

/// Very hacky, see [StrideStorage] for explanation.
/// Metadata stores the number of rows (in the least significant 4 bytes)
/// and columns (in the most significant 4 bytes).
pub struct SliceMatrixStorage<R: Ring>(
    CustomMetadataSlice<R::Element, SliceMatrixStorageMetadata>,
);

struct SliceMatrixStorageMetadata {
    rows: Half,
    cols: Half,
}

impl SliceMatrixStorageMetadata {
    pub fn new(rows: usize, cols: usize) -> Self {
        Self { rows: rows.try_into().unwrap(), cols: cols.try_into().unwrap() }
    }
}

impl CustomMetadata for SliceMatrixStorageMetadata {
    fn size(&self) -> usize {
        self.rows as usize * self.cols as usize
    }
}

impl<R: Ring> MatrixStorage<R> for SliceMatrixStorage<R> {
    type ColIter<'a> = MatrixStridedIter<'a, R>;
    type ColIterMut<'a> = MatrixStridedIterMut<'a, R>;
    type ColVecStorage = StrideStorage<R>;
    type RowIter<'a> = MatrixSliceIter<'a, R>;
    type RowIterMut<'a> = MatrixSliceIterMut<'a, R>;
    type RowVecStorage = SliceVectorStorage<R>;

    fn rows(&self) -> usize {
        self.0.metadata().rows as usize
    }

    fn cols(&self) -> usize {
        self.0.metadata().cols as usize
    }

    fn row(&self, r: usize) -> &Self::RowVecStorage {
        assert!(r < self.rows());
        unsafe {
            SliceVectorStorage::from_slice(std::slice::from_raw_parts(
                self.0.as_ptr().add(r * self.cols()),
                self.cols(),
            ))
        }
    }

    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage {
        assert!(r < self.rows());
        unsafe {
            SliceVectorStorage::from_slice_mut(std::slice::from_raw_parts_mut(
                self.0.as_mut_ptr().add(r * self.cols()),
                self.cols(),
            ))
        }
    }

    fn row_iter(&self) -> Self::RowIter<'_> {
        unsafe {
            MatrixSliceIter::new(self.0.as_ptr(), self.cols(), self.rows())
        }
    }

    fn row_iter_mut(&mut self) -> Self::RowIterMut<'_> {
        unsafe {
            MatrixSliceIterMut::new(
                self.0.as_mut_ptr(),
                self.cols(),
                self.rows(),
            )
        }
    }

    fn col(&self, c: usize) -> &Self::ColVecStorage {
        assert!(c < self.cols());
        unsafe {
            StrideStorage::from_raw_parts(
                self.0.as_ptr().add(c),
                self.rows(),
                self.cols(),
            )
        }
    }

    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage {
        assert!(c < self.cols());
        unsafe {
            StrideStorage::from_raw_parts_mut(
                self.0.as_mut_ptr().add(c),
                self.rows(),
                self.cols(),
            )
        }
    }

    fn col_iter(&self) -> Self::ColIter<'_> {
        unsafe {
            MatrixStridedIter::new(self.0.as_ptr(), self.rows(), self.cols())
        }
    }

    fn col_iter_mut(&mut self) -> Self::ColIterMut<'_> {
        unsafe {
            MatrixStridedIterMut::new(
                self.0.as_mut_ptr(),
                self.rows(),
                self.cols(),
            )
        }
    }
}

pub type TransposedMatrixView<R> = Matrix<R, TransposedMatrixStorage<R>>;

impl<R: Ring> TransposedMatrixView<R> {
    /// Returns a view of the matrix transposed.
    pub fn transpose(&self) -> &MatrixView<R> {
        unsafe { &*(self as *const _ as *const _) }
    }

    /// Returns a mutable view of the matrix transposed.
    pub fn transpose_mut(&mut self) -> &mut MatrixView<R> {
        unsafe { &mut *(self as *mut _ as *mut _) }
    }
}

/// Very hacky, see [StrideStorage] for explanation.
/// Metadata stores the number of columns (in the least significant 4 bytes)
/// and rows (in the most significant 4 bytes).
/// This allows you to transmute a [SliceMatrixStorage] reference into a
/// [TransposedMatrixStorage] reference while transposing the view.
pub struct TransposedMatrixStorage<R: Ring>(
    CustomMetadataSlice<R::Element, TransposedMatrixStorageMetadata>,
);

struct TransposedMatrixStorageMetadata {
    cols: Half,
    rows: Half,
}

impl CustomMetadata for TransposedMatrixStorageMetadata {
    fn size(&self) -> usize {
        self.rows as usize * self.cols as usize
    }
}

impl<R: Ring> MatrixStorage<R> for TransposedMatrixStorage<R> {
    type ColIter<'a> = MatrixSliceIter<'a, R>;
    type ColIterMut<'a> = MatrixSliceIterMut<'a, R>;
    type ColVecStorage = SliceVectorStorage<R>;
    type RowIter<'a> = MatrixStridedIter<'a, R>;
    type RowIterMut<'a> = MatrixStridedIterMut<'a, R>;
    type RowVecStorage = StrideStorage<R>;

    fn rows(&self) -> usize {
        self.0.metadata().rows as usize
    }

    fn cols(&self) -> usize {
        self.0.metadata().cols as usize
    }

    fn row(&self, r: usize) -> &Self::RowVecStorage {
        assert!(r < self.rows());
        unsafe {
            StrideStorage::from_raw_parts(
                self.0.as_ptr().add(r),
                self.cols(),
                self.rows(),
            )
        }
    }

    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage {
        assert!(r < self.rows());
        unsafe {
            StrideStorage::from_raw_parts_mut(
                self.0.as_mut_ptr().add(r),
                self.cols(),
                self.rows(),
            )
        }
    }

    fn row_iter(&self) -> Self::RowIter<'_> {
        unsafe {
            MatrixStridedIter::new(self.0.as_ptr(), self.rows(), self.cols())
        }
    }

    fn row_iter_mut(&mut self) -> Self::RowIterMut<'_> {
        unsafe {
            MatrixStridedIterMut::new(
                self.0.as_mut_ptr(),
                self.rows(),
                self.cols(),
            )
        }
    }

    fn col(&self, c: usize) -> &Self::ColVecStorage {
        assert!(c < self.cols());
        unsafe {
            SliceVectorStorage::from_slice(std::slice::from_raw_parts(
                self.0.as_ptr().add(c * self.cols()),
                self.rows(),
            ))
        }
    }

    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage {
        assert!(c < self.cols());
        unsafe {
            SliceVectorStorage::from_slice_mut(std::slice::from_raw_parts_mut(
                self.0.as_mut_ptr().add(c * self.cols()),
                self.rows(),
            ))
        }
    }

    fn col_iter(&self) -> Self::ColIter<'_> {
        unsafe {
            MatrixSliceIter::new(self.0.as_ptr(), self.cols(), self.rows())
        }
    }

    fn col_iter_mut(&mut self) -> Self::ColIterMut<'_> {
        unsafe {
            MatrixSliceIterMut::new(
                self.0.as_mut_ptr(),
                self.cols(),
                self.rows(),
            )
        }
    }
}

/// A matrix with a single row.
pub type RowVector<R, S> = Matrix<R, RowVectorStorage<R, S>>;

pub struct RowVectorStorage<R: Ring, S: VectorStorage<R> + ?Sized>(
    PhantomData<R>,
    S,
);

impl<R, S> MatrixStorage<R> for RowVectorStorage<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    type ColIter<'a> = std::iter::Map<
        S::Iter<'a>,
        fn(&'a R::Element) -> &'a Vector<R, SliceVectorStorage<R>>,
    >;
    type ColIterMut<'a> = std::iter::Map<
        S::IterMut<'a>,
        fn(&'a mut R::Element) -> &'a mut Vector<R, SliceVectorStorage<R>>,
    >;
    type ColVecStorage = SliceVectorStorage<R>;
    type RowIter<'a> = std::iter::Once<&'a Vector<R, S>>;
    type RowIterMut<'a> = std::iter::Once<&'a mut Vector<R, S>>;
    type RowVecStorage = S;

    fn rows(&self) -> usize {
        1
    }

    fn cols(&self) -> usize {
        self.1.dim()
    }

    fn row(&self, r: usize) -> &Self::RowVecStorage {
        assert!(r < self.rows());
        &self.1
    }

    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage {
        assert!(r < self.rows());
        &mut self.1
    }

    fn row_iter(&self) -> Self::RowIter<'_> {
        std::iter::once(Vector::from_storage_ref(self.row(0)))
    }

    fn row_iter_mut(&mut self) -> Self::RowIterMut<'_> {
        std::iter::once(Vector::from_storage_ref_mut(self.row_mut(0)))
    }

    fn col(&self, c: usize) -> &Self::ColVecStorage {
        assert!(c < self.cols());
        SliceVectorStorage::from_slice(std::slice::from_ref(self.1.entry(c)))
    }

    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage {
        assert!(c < self.cols());
        SliceVectorStorage::from_slice_mut(std::slice::from_mut(
            self.1.entry_mut(c),
        ))
    }

    fn col_iter(&self) -> Self::ColIter<'_> {
        self.1.iter().map(VectorView::from_elem)
    }

    fn col_iter_mut(&mut self) -> Self::ColIterMut<'_> {
        self.1.iter_mut().map(VectorView::from_elem_mut)
    }
}

/// A matrix with a single column.
pub type ColumnVector<R, S> = Matrix<R, ColumnVectorStorage<R, S>>;

pub struct ColumnVectorStorage<R: Ring, S: VectorStorage<R> + ?Sized>(
    PhantomData<R>,
    S,
);

impl<R, S> MatrixStorage<R> for ColumnVectorStorage<R, S>
where
    R: Ring,
    S: VectorStorage<R> + ?Sized,
{
    type ColIter<'a> = std::iter::Once<&'a Vector<R, S>>;
    type ColIterMut<'a> = std::iter::Once<&'a mut Vector<R, S>>;
    type ColVecStorage = S;
    type RowIter<'a> = std::iter::Map<
        S::Iter<'a>,
        fn(&'a R::Element) -> &'a Vector<R, SliceVectorStorage<R>>,
    >;
    type RowIterMut<'a> = std::iter::Map<
        S::IterMut<'a>,
        fn(&'a mut R::Element) -> &'a mut Vector<R, SliceVectorStorage<R>>,
    >;
    type RowVecStorage = SliceVectorStorage<R>;

    fn rows(&self) -> usize {
        self.1.dim()
    }

    fn cols(&self) -> usize {
        1
    }

    fn row(&self, r: usize) -> &Self::RowVecStorage {
        assert!(r < self.rows());
        SliceVectorStorage::from_slice(std::slice::from_ref(self.1.entry(r)))
    }

    fn row_mut(&mut self, r: usize) -> &mut Self::RowVecStorage {
        assert!(r < self.rows());
        SliceVectorStorage::from_slice_mut(std::slice::from_mut(
            self.1.entry_mut(r),
        ))
    }

    fn row_iter(&self) -> Self::RowIter<'_> {
        self.1.iter().map(VectorView::from_elem)
    }

    fn row_iter_mut(&mut self) -> Self::RowIterMut<'_> {
        self.1.iter_mut().map(VectorView::from_elem_mut)
    }

    fn col(&self, c: usize) -> &Self::ColVecStorage {
        assert!(c < self.cols());
        &self.1
    }

    fn col_mut(&mut self, c: usize) -> &mut Self::ColVecStorage {
        assert!(c < self.cols());
        &mut self.1
    }

    fn col_iter(&self) -> Self::ColIter<'_> {
        std::iter::once(Vector::from_storage_ref(self.col(0)))
    }

    fn col_iter_mut(&mut self) -> Self::ColIterMut<'_> {
        std::iter::once(Vector::from_storage_ref_mut(self.col_mut(0)))
    }
}

/// Iterate over vectors of the matrix that are contiguous in memory.
/// This is usually rows for [`OwnedMatrix`] and [`MatrixView`] and columns for
/// [`TransposedMatrixView`].
pub struct MatrixSliceIter<'a, R: Ring> {
    phantom: PhantomData<&'a R::Element>,
    dim: usize,
    start_ptr: *const R::Element,
    end_ptr: *const R::Element,
}

impl<'a, R: Ring> MatrixSliceIter<'a, R> {
    /// Create a new iterator.
    /// - `base` is the base address of the slice of matrix elements.
    /// - `vdim` is the dimension of the returned vectors, which is either the
    ///   number of columns or rows. For [`OwnedMatrix`] and [`MatrixView`] it
    ///   is the number of columns, for [`TransposedMatrixView`] it is the
    ///   number of rows.
    /// - `odim` is the other dimension.
    ///
    /// The slice of elements should have at least `vdim * odim` elements.
    unsafe fn new(base: *const R::Element, vdim: usize, odim: usize) -> Self {
        Self {
            phantom: PhantomData,
            dim: vdim,
            start_ptr: base,
            end_ptr: unsafe { base.add(vdim * odim) },
        }
    }
}

impl<'a, R: Ring> Iterator for MatrixSliceIter<'a, R> {
    type Item = &'a Vector<R, SliceVectorStorage<R>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        let slice =
            unsafe { std::slice::from_raw_parts(self.start_ptr, self.dim) };

        self.start_ptr = unsafe { self.start_ptr.add(self.dim) };

        Some(Vector::from_storage_ref(SliceVectorStorage::from_slice(slice)))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.start_ptr = unsafe { self.start_ptr.add(self.dim * n) };
        self.next()
    }
}

impl<'a, R: Ring> DoubleEndedIterator for MatrixSliceIter<'a, R> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        self.end_ptr = unsafe { self.end_ptr.sub(self.dim) };
        let slice =
            unsafe { std::slice::from_raw_parts(self.end_ptr, self.dim) };

        Some(Vector::from_storage_ref(SliceVectorStorage::from_slice(slice)))
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.end_ptr = unsafe { self.end_ptr.sub(self.dim * n) };
        self.next()
    }
}

/// Mutable version of [`MatrixSliceIter`].
pub struct MatrixSliceIterMut<'a, R: Ring> {
    phantom: PhantomData<&'a mut R::Element>,
    dim: usize,
    start_ptr: *mut R::Element,
    end_ptr: *mut R::Element,
}

impl<'a, R: Ring> MatrixSliceIterMut<'a, R> {
    /// See [`MatrixSliceIter::new`].
    unsafe fn new(base: *mut R::Element, vdim: usize, odim: usize) -> Self {
        Self {
            phantom: PhantomData,
            dim: vdim,
            start_ptr: base,
            end_ptr: unsafe { base.add(vdim * odim) },
        }
    }
}

impl<'a, R: Ring> Iterator for MatrixSliceIterMut<'a, R> {
    type Item = &'a mut Vector<R, SliceVectorStorage<R>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        let slice =
            unsafe { std::slice::from_raw_parts_mut(self.start_ptr, self.dim) };

        self.start_ptr = unsafe { self.start_ptr.add(self.dim) };

        Some(Vector::from_storage_ref_mut(SliceVectorStorage::from_slice_mut(
            slice,
        )))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.start_ptr = unsafe { self.start_ptr.add(self.dim * n) };
        self.next()
    }
}

impl<'a, R: Ring> DoubleEndedIterator for MatrixSliceIterMut<'a, R> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        self.end_ptr = unsafe { self.end_ptr.sub(self.dim) };
        let slice =
            unsafe { std::slice::from_raw_parts_mut(self.end_ptr, self.dim) };

        Some(Vector::from_storage_ref_mut(SliceVectorStorage::from_slice_mut(
            slice,
        )))
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.end_ptr = unsafe { self.end_ptr.sub(self.dim * n) };
        self.next()
    }
}

/// Iterate over vectors of the matrix that are strided in memory.
/// This is usually columns for [`OwnedMatrix`] and [`MatrixView`] and rows for
/// [`TransposedMatrixView`].
pub struct MatrixStridedIter<'a, R: Ring> {
    phantom: PhantomData<&'a R::Element>,
    vdim: usize,
    odim: usize,
    start_ptr: *const R::Element,
    end_ptr: *const R::Element,
}

impl<'a, R: Ring> MatrixStridedIter<'a, R> {
    /// Create a new iterator.
    /// - `base` is the base address of the slice of matrix elements.
    /// - `vdim` is the dimension of the returned vectors, which is either the
    ///   number of columns or rows. For [`OwnedMatrix`] and [`MatrixView`] it
    ///   is the number of rows, for [`TransposedMatrixView`] it is the number
    ///   of columns.
    /// - `odim` is the other dimension.
    ///
    /// The slice of elements should have at least `vdim * odim` elements.
    unsafe fn new(base: *const R::Element, vdim: usize, odim: usize) -> Self {
        Self {
            phantom: PhantomData,
            vdim,
            odim,
            start_ptr: base,
            end_ptr: unsafe { base.add(odim) },
        }
    }
}

impl<'a, R: Ring> Iterator for MatrixStridedIter<'a, R> {
    type Item = &'a Vector<R, StrideStorage<R>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        let s = unsafe {
            StrideStorage::from_raw_parts(self.start_ptr, self.vdim, self.odim)
        };

        self.start_ptr = unsafe { self.start_ptr.add(1) };

        Some(Vector::from_storage_ref(s))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.start_ptr = unsafe { self.start_ptr.add(n) };
        self.next()
    }
}

impl<'a, R: Ring> DoubleEndedIterator for MatrixStridedIter<'a, R> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        self.end_ptr = unsafe { self.end_ptr.sub(1) };

        let s = unsafe {
            StrideStorage::from_raw_parts(self.end_ptr, self.vdim, self.odim)
        };

        Some(Vector::from_storage_ref(s))
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.end_ptr = unsafe { self.end_ptr.sub(n) };
        self.next()
    }
}

/// See [`MatrixStridedIter`].
pub struct MatrixStridedIterMut<'a, R: Ring> {
    phantom: PhantomData<&'a mut R::Element>,
    vdim: usize,
    odim: usize,
    start_ptr: *mut R::Element,
    end_ptr: *mut R::Element,
}

impl<'a, R: Ring> MatrixStridedIterMut<'a, R> {
    /// See [`MatrixStridedIter::new`].
    unsafe fn new(base: *mut R::Element, vdim: usize, odim: usize) -> Self {
        Self {
            phantom: PhantomData,
            vdim,
            odim,
            start_ptr: base,
            end_ptr: unsafe { base.add(odim) },
        }
    }
}

impl<'a, R: Ring> Iterator for MatrixStridedIterMut<'a, R> {
    type Item = &'a mut Vector<R, StrideStorage<R>>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        let s = unsafe {
            StrideStorage::from_raw_parts_mut(
                self.start_ptr,
                self.vdim,
                self.odim,
            )
        };

        self.start_ptr = unsafe { self.start_ptr.add(1) };

        Some(Vector::from_storage_ref_mut(s))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.start_ptr = unsafe { self.start_ptr.add(n) };
        self.next()
    }
}

impl<'a, R: Ring> DoubleEndedIterator for MatrixStridedIterMut<'a, R> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.start_ptr >= self.end_ptr {
            return None;
        }

        self.end_ptr = unsafe { self.end_ptr.sub(1) };

        let s = unsafe {
            StrideStorage::from_raw_parts_mut(
                self.end_ptr,
                self.vdim,
                self.odim,
            )
        };

        Some(Vector::from_storage_ref_mut(s))
    }

    fn nth_back(&mut self, n: usize) -> Option<Self::Item> {
        self.end_ptr = unsafe { self.end_ptr.sub(n) };
        self.next()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::rings::U32;
    #[test]
    fn transpose_test() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let t = m.transpose();
        assert_eq!(t.num_rows(), 2);
        assert_eq!(t.num_cols(), 2);
        assert_eq!(t[0], [2, 4]);
        assert_eq!(t[1], [3, 5]);
        assert_eq!(t[(0, 0)], 2);
        assert_eq!(t[(0, 1)], 4);
        assert_eq!(t[(1, 0)], 3);
        assert_eq!(t[(1, 1)], 5);
    }

    #[test]
    fn row() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        assert_eq!(m.row(0), &[2, 3]);
        assert_eq!(m.row(1), &[4, 5]);
    }

    #[test]
    fn col() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        assert_eq!(m.col(0), &[2, 4]);
        assert_eq!(m.col(1), &[3, 5]);
    }

    #[test]
    fn rows_iter() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let mut r = m.rows();
        assert_eq!(r.next().unwrap(), &[2, 3]);
        assert_eq!(r.next().unwrap(), &[4, 5]);
        assert!(r.next().is_none());
    }

    #[test]
    fn col_iter() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let mut c0 = m.col(0).iter();
        assert_eq!(c0.next(), Some(&2));
        assert_eq!(c0.next(), Some(&4));
        assert!(c0.next().is_none());
        let mut c1 = m.col(1).iter();
        assert_eq!(c1.next(), Some(&3));
        assert_eq!(c1.next(), Some(&5));
        assert!(c0.next().is_none());
    }

    #[test]
    fn row_iter_rev_test() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let mut r = m.rows();
        assert_eq!(r.next_back().unwrap(), &[4, 5]);
        assert_eq!(r.next_back().unwrap(), &[2, 3]);
        assert!(r.next_back().is_none());
    }

    #[test]
    fn col_iter_rev() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let mut r = m.cols();
        assert_eq!(r.next_back().unwrap(), &[3, 5]);
        assert_eq!(r.next_back().unwrap(), &[2, 4]);
        assert!(r.next_back().is_none());
    }

    #[test]
    fn entry_iterators() {
        let m = OwnedMatrix::<U32>::from_rows(&[[2u32, 3u32], [4u32, 5u32]]);
        let mut r = m.entries_row_major();
        assert_eq!(r.next(), Some(&2));
        assert_eq!(r.next(), Some(&3));
        assert_eq!(r.next(), Some(&4));
        assert_eq!(r.next(), Some(&5));
        assert!(r.next().is_none());
        let mut r = m.entries_col_major();
        assert_eq!(r.next(), Some(&2));
        assert_eq!(r.next(), Some(&4));
        assert_eq!(r.next(), Some(&3));
        assert_eq!(r.next(), Some(&5));
        assert!(r.next().is_none());
    }
}
