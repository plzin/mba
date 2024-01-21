//! Matrices

use std::ops::{Index, IndexMut};
use std::{marker::PhantomData, fmt::Debug, ops::Mul};
use itertools::iproduct;
use num_traits::{Zero, One};
use num_bigint::BigInt;
use num_rational::BigRational;
use crate::{vector::*, CustomMetadataSlice, CustomMetadata};

/// How are the entries of a matrix stored?
pub trait MatrixStorage<T> {
    /// Row/Column vector type.
    type RowVecStorage: VectorStorage<T> + ?Sized;
    type ColVecStorage: VectorStorage<T> + ?Sized;

    /// Returns the number of rows.
    fn rows(&self) -> usize;

    /// Returns the number of columns.
    fn cols(&self) -> usize;

    /// Returns a view of the row `r`.
    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage>;

    /// Returns a mutable view of the row `r`.
    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage>;

    /// Returns a view of the column `c`.
    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage>;

    /// Returns a mutable view of the column `c`.
    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage>;

    /// Returns a reference to the entry at row `r` and column `c`.
    fn entry(&self, r: usize, c: usize) -> &T {
        self.row(r).entry(c)
    }

    /// Returns a mutable reference to the entry at row `r` and column `c`.
    fn entry_mut(&mut self, r: usize, c: usize) -> &mut T {
        self.row_mut(r).entry_mut(c)
    }
}

/// A matrix.
pub struct Matrix<T, S: MatrixStorage<T> + ?Sized> {
    marker: PhantomData<T>,
    storage: S,
}

impl<T, S: MatrixStorage<T> + ?Sized> Matrix<T, S> {
    /// The number of rows of the matrix.
    pub fn nrows(&self) -> usize {
        self.storage.rows()
    }

    /// The number of columns of the matrix.
    pub fn ncols(&self) -> usize {
        self.storage.cols()
    }

    /// Is the matrix empty, i.e. has it zero rows or columns?
    pub fn is_empty(&self) -> bool {
        self.nrows() == 0 || self.ncols() == 0
    }

    /// Returns a view of the row `r`.
    pub fn row(&self, r: usize) -> &Vector<T, S::RowVecStorage> {
        self.storage.row(r)
    }

    /// Returns a mutable view of the row `r`.
    pub fn row_mut(&mut self, r: usize) -> &mut Vector<T, S::RowVecStorage> {
        self.storage.row_mut(r)
    }

    /// Returns a view of the column `c`.
    pub fn col(&self, c: usize) -> &Vector<T, S::ColVecStorage> {
        self.storage.col(c)
    }

    /// Returns a mutable view of the column `c`.
    pub fn col_mut(&mut self, c: usize) -> &mut Vector<T, S::ColVecStorage> {
        self.storage.col_mut(c)
    }

    /// Returns a reference to the entry at row `r` and column `c`.
    pub fn entry(&self, r: usize, c: usize) -> &T {
        self.storage.entry(r, c)
    }

    /// Returns a mutable reference to the entry at row `r` and column `c`.
    pub fn entry_mut(&mut self, r: usize, c: usize) -> &mut T {
        self.storage.entry_mut(r, c)
    }

    /// Returns an iterator over the rows.
    pub fn rows(&self) -> impl DoubleEndedIterator<Item = &Vector<T, S::RowVecStorage>> {
        (0..self.nrows()).map(|r| self.row(r))
    }

    /// Returns an iterator over the mutable rows.
    pub fn rows_mut(&mut self) -> impl DoubleEndedIterator<Item = &mut Vector<T, S::RowVecStorage>> {
        (0..self.nrows()).map(|r| unsafe { &mut *(self as *mut Self) }.row_mut(r))
    }

    /// Returns an iterator over the columns.
    pub fn cols(&self) -> impl DoubleEndedIterator<Item = &Vector<T, S::ColVecStorage>> {
        (0..self.ncols()).map(|c| self.col(c))
    }

    /// Returns an iterator over the mutable columns.
    pub fn cols_mut(&mut self) -> impl DoubleEndedIterator<Item = &mut Vector<T, S::ColVecStorage>> {
        (0..self.ncols()).map(|c| unsafe { &mut *(self as *mut Self) }.col_mut(c))
    }

    /// Returns an iterator over the entries in row-major order.
    pub fn entries_row_major(&self) -> impl Iterator<Item = &T> {
        self.rows().flat_map(|r| r.iter())
    }

    /// Returns an iterator over the mutable entries in row-major order.
    pub fn entries_row_major_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.rows_mut().flat_map(|r| r.iter_mut())
    }


    /// Returns an iterator over the entries in column-major order.
    pub fn entries_col_major(&self) -> impl Iterator<Item = &T> {
        self.cols().flat_map(|c| c.iter())
    }

    /// Returns an iterator over the mutable entries in column-major order.
    pub fn entries_col_major_mut(&mut self) -> impl Iterator<Item = &mut T> {
        self.cols_mut().flat_map(|c| c.iter_mut())
    }

    /// Call a function on each entry and return the resulting matrix.
    pub fn transform<U, F: FnMut(&T) -> U>(&self, f: F) -> OwnedMatrix<U> {
        OwnedMatrix::from_iter(
            self.nrows(),
            self.ncols(),
            self.entries_row_major().map(f),
        )
    }
}

impl<T, S: MatrixStorage<T>> Matrix<T, S> {
    fn from_storage(s: S) -> Self {
        Self {
            marker: PhantomData,
            storage: s,
        }
    }
}

impl<T: Clone, S: MatrixStorage<T> + ?Sized> Matrix<T, S> {
    /// Returns a copy of the matrix.
    pub fn to_owned(&self) -> OwnedMatrix<T> {
        OwnedMatrix::from_iter(
            self.nrows(),
            self.ncols(),
            self.entries_row_major().cloned()
        )
    }

    /// Creates an owned matrix that is the transpose of the current matrix.
    pub fn transposed(&self) -> OwnedMatrix<T> {
        OwnedMatrix::from_iter(
            self.ncols(),
            self.nrows(),
            self.entries_col_major().cloned()
        )
    }
}

impl<S: MatrixStorage<BigInt> + ?Sized> Matrix<BigInt, S> {
    /// Flip the sign of a row.
    pub fn flip_sign_row(&mut self, row: usize) {
        for e in self.row_mut(row) {
            crate::neg_mut(e);
        }
    }

    /// Flip the sign of a column.
    pub fn flip_sign_column(&mut self, col: usize) {
        for e in self.col_mut(col) {
            crate::neg_mut(e);
        }
    }

    /// Add a scaled row to another row. N = c * M.
    pub fn row_multiply_add(&mut self, n: usize, m: usize, c: &BigInt) {
        debug_assert!(n < self.nrows() && m < self.nrows());
        for i in 0..self.ncols() {
            let c = &self[(n, i)] * c;
            self[(m, i)] += c;
        }
    }

    /// Add a scaled column to another column. N = c * M.
    pub fn col_multiply_add(&mut self, n: usize, m: usize, c: &BigInt) {
        debug_assert!(n < self.ncols() && m < self.ncols());
        for i in 0..self.nrows() {
            let c = &self[(i, n)] * c;
            self[(i, m)] += c;
        }
    }
}

impl<T, S: MatrixStorage<T> + ?Sized> Index<usize> for Matrix<T, S> {
    type Output = Vector<T, S::RowVecStorage>;

    fn index(&self, index: usize) -> &Self::Output {
        self.row(index)
    }
}

impl<T, S: MatrixStorage<T> + ?Sized> IndexMut<usize> for Matrix<T, S> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.row_mut(index)
    }
}

impl<T, S: MatrixStorage<T> + ?Sized> Index<(usize, usize)> for Matrix<T, S> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        self.entry(index.0, index.1)
    }

}

impl<T, S: MatrixStorage<T> + ?Sized> IndexMut<(usize, usize)> for Matrix<T, S> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        self.entry_mut(index.0, index.1)
    }
}

pub type OwnedMatrix<T> = Matrix<T, OwnedMatrixStorage<T>>;
pub type IOwnedMatrix = OwnedMatrix<BigInt>;
pub type QOwnedMatrix = OwnedMatrix<BigRational>;

impl<T> OwnedMatrix<T> {
    /// Return an empty matrix.
    pub fn empty() -> Self {
        Self::from_raw_parts(core::ptr::null_mut(), 0, 0)
    }

    /// Returns a slice containing all elements.
    pub fn as_slice(&self) -> &[T] {
        self.storage.as_slice()
    }

    /// Returns a slice containing all elements.
    pub fn as_slice_mut(&mut self) -> &mut [T] {
        self.storage.as_slice_mut()
    }

    /// Returns a view of this matrix.
    pub fn view(&self) -> &MatrixView<T> {
        MatrixView::from_raw_parts(
            self.storage.entries, self.nrows(), self.ncols()
        )
    }

    /// Returns a mutable view of this matrix.
    pub fn view_mut(&mut self) -> &mut MatrixView<T> {
        MatrixView::from_raw_parts_mut(
            self.storage.entries, self.nrows(), self.ncols()
        )
    }

    /// Turns this matrix into a vector of dimension r*c.
    pub fn into_vector(self) -> Vector<T, OwnedVectorStorage<T>> {
        let entries = self.storage.entries;
        let dim = self.nrows() * self.ncols();

        // Don't call the destructor on self.
        std::mem::forget(self);

        // Transfer the memory to the vector.
        // It has a matching destructor.
        OwnedVector::from_raw_parts(entries, dim)
    }

    fn from_raw_parts(entries: *mut T, rows: usize, cols: usize) -> Self {
        Self::from_storage(OwnedMatrixStorage {
            entries,
            rows,
            cols,
        })
    }

    /// Returns an uninitialized r×c matrix.
    pub(self) fn uninit(r: usize, c: usize) -> Self {
        if r == 0 || c == 0 {
            return Self::empty();
        }

        // The memory layout of the elements.
        let layout = std::alloc::Layout::from_size_align(
            r * c * core::mem::size_of::<T>(),
            core::mem::align_of::<T>()
        ).unwrap();

        // Allocate memory for the entries.
        let entries = unsafe {
            std::alloc::alloc(layout) as *mut T
        };

        Self::from_raw_parts(entries, r, c)
    }

    /// Creates a matrix from an array of rows.
    pub fn from_array<U: Into<T>, const R: usize, const C: usize>(
        a: [[U; C]; R]
    ) -> Self {
        let m = Self::uninit(R, C);

        // Copy the entries.
        let mut ptr = m.storage.entries;
        for row in a {
            for e in row {
                unsafe {
                    ptr.write(e.into());
                    ptr = ptr.add(1);
                }
            }
        }

        m
    }

    /// Creates a matrix from an iterator in row-major order.
    pub fn from_iter<I: Iterator<Item = T>>(
        r: usize, c: usize, mut iter: I
    ) -> Self {
        let m = Self::uninit(r, c);
        for i in 0..r*c {
            let e = iter.next()
                .expect("The iterator needs to return at least r * c items.");
            unsafe {
                m.storage.entries.add(i).write(e);
            }
        }

        m
    }

    /// Creates a matrix from slice of rows.
    pub fn from_rows<U, V>(rows: &[U]) -> Self
    where
        U: AsRef<[V]>,
        V: Into<T> + Clone,
    {
        assert!(!rows.is_empty());
        let r = rows.len();
        let c = rows[0].as_ref().len();
        assert!(rows.iter().all(|r| r.as_ref().len() == c));

        Matrix::from_iter(r, c, rows.iter()
            .flat_map(|r| r.as_ref().iter().map(|e| e.clone().into()))
        )
    }

    /// View of the transpose.
    pub fn transpose(&self) -> &TransposedMatrixView<T> {
        self.view().transpose()
    }

    /// Mutable view of the transpose.
    pub fn transpose_mut(&mut self) -> &mut TransposedMatrixView<T> {
        self.view_mut().transpose_mut()
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
                self.ncols()
            )
        }
    }

    /// Swap two columns.
    pub fn swap_columns(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        for k in 0..self.ncols() {
            unsafe {
                core::ptr::swap_nonoverlapping(
                    self.entry_mut(i, k), self.entry_mut(j, k), 1
                );
            }
        }
    }

    /// Removes rows at the end of the matrix.
    pub fn shrink(&mut self, new_nrows: usize) {
        self.storage.shrink(new_nrows);
    }
}

impl<T: Zero + One> OwnedMatrix<T> {
    /// Returns an r×c zero matrix.
    pub fn zero(r: usize, c: usize) -> Self {
        Self::from_iter(r, c, std::iter::repeat_with(|| T::zero()))
    }

    /// Returns an nxn identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zero(n, n);
        for i in 0..n {
            m[(i, i)] = T::one();
        }
        m
    }
}

pub struct OwnedMatrixStorage<T> {
    /// Memory that holds the entries.
    pub(self) entries: *mut T,

    /// The number of rows.
    pub(self) rows: usize,

    /// The number of columns.
    pub(self) cols: usize,
}

impl<T> OwnedMatrixStorage<T> {
    /// Returns a slice containing all elements.
    pub fn as_slice(&self) -> &[T] {
        unsafe {
            std::slice::from_raw_parts(self.entries, self.rows*self.cols)
        }
    }

    /// Returns a slice containing all elements.
    pub fn as_slice_mut(&mut self) -> &mut [T] {
        unsafe {
            std::slice::from_raw_parts_mut(self.entries, self.rows*self.cols)
        }
    }

    /// Removes rows at the end of the matrix.
    pub fn shrink(&mut self, new_nrows: usize) {
        if new_nrows == self.rows {
            return;
        }
        assert!(new_nrows < self.rows);

        // Call the destructor on all elements.
        for i in new_nrows*self.cols..self.rows*self.cols {
            unsafe {
                core::ptr::drop_in_place(self.entries.add(i));
            }
        }

        let esize = core::mem::size_of::<T>();
        let layout = std::alloc::Layout::from_size_align(
            self.rows * self.cols * esize,
            core::mem::align_of::<T>()
        ).unwrap();

        // This should hopefully not actually allocate.
        let new_ptr = unsafe {
            std::alloc::realloc(
                self.entries as _,
                layout,
                new_nrows * self.cols * esize
            )
        };
        self.entries = new_ptr as _;
        self.rows = new_nrows;
    }
}

impl<T> Drop for OwnedMatrixStorage<T> {
    fn drop(&mut self) {
        if self.entries.is_null() {
            return;
        }

        unsafe {
            // Call the destructor on all elements.
            for i in 0..self.rows*self.cols {
                core::ptr::drop_in_place(self.entries.add(i));
            }

            let layout = std::alloc::Layout::from_size_align(
                self.rows * self.cols * core::mem::size_of::<T>(),
                core::mem::align_of::<T>()
            ).unwrap();

            // Free the memory.
            std::alloc::dealloc(self.entries as _, layout);
        }
    }
}

impl<T> MatrixStorage<T> for OwnedMatrixStorage<T> {
    type RowVecStorage = SliceVectorStorage<T>;
    type ColVecStorage = StrideStorage<T>;

    fn rows(&self) -> usize {
        self.rows
    }

    fn cols(&self) -> usize {
        self.cols
    }

    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows);
        unsafe {
            VectorView::from_slice(std::slice::from_raw_parts(
                self.entries.add(r * self.cols),
                self.cols
            ))
        }
    }

    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows);
        unsafe {
            VectorView::from_slice_mut(std::slice::from_raw_parts_mut(
                self.entries.add(r * self.cols),
                self.cols
            ))
        }
    }

    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols);
        StrideVectorView::from_raw_parts(
            unsafe { self.entries.add(c) },
            self.rows,
            self.cols
        )
    }

    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols);
        StrideVectorView::from_raw_parts_mut(
            unsafe { self.entries.add(c) },
            self.rows,
            self.cols
        )
    }
}

pub type MatrixView<T> = Matrix<T, SliceMatrixStorage<T>>;
pub type IMatrixView = MatrixView<BigInt>;
pub type RMatrixView = MatrixView<BigRational>;

impl<T> MatrixView<T> {
    /// Returns a view of the matrix transposed.
    pub fn transpose(&self) -> &TransposedMatrixView<T> {
        unsafe {
            &*(self as *const _ as *const _)
        }
    }

    /// Returns a mutable view of the matrix transposed.
    pub fn transpose_mut(&mut self) -> &mut TransposedMatrixView<T> {
        unsafe {
            &mut *(self as *mut _ as *mut _)
        }
    }

    pub fn from_raw_parts<'a>(entries: *const T, rows: usize, cols: usize) -> &'a Self {
        let metadata = SliceMatrixStorageMetadata { rows, cols };
        unsafe {
            std::mem::transmute(CustomMetadataSlice::new(entries, metadata))
        }
    }

    pub fn from_raw_parts_mut<'a>(entries: *mut T, rows: usize, cols: usize) -> &'a mut Self {
        let metadata = SliceMatrixStorageMetadata { rows, cols };
        unsafe {
            std::mem::transmute(CustomMetadataSlice::new_mut(entries, metadata))
        }
    }
}

/// Very hacky, see [StrideStorage] for explanation.
/// Metadata stores the number of rows (in the least significant 4 bytes)
/// and columns (in the most significant 4 bytes).
pub struct SliceMatrixStorage<T>(
    CustomMetadataSlice<T, SliceMatrixStorageMetadata>
);

struct SliceMatrixStorageMetadata {
    rows: usize,
    cols: usize,
}

impl CustomMetadata for SliceMatrixStorageMetadata {
    fn size(&self) -> usize {
        self.rows * self.cols
    }
}

impl<T> MatrixStorage<T> for SliceMatrixStorage<T> {
    type RowVecStorage = SliceVectorStorage<T>;
    type ColVecStorage = StrideStorage<T>;

    fn rows(&self) -> usize {
        self.0.metadata().rows
    }

    fn cols(&self) -> usize {
        self.0.metadata().cols
    }

    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        unsafe {
            VectorView::from_slice(std::slice::from_raw_parts(
                self.0.as_ptr().add(r * self.cols()),
                self.cols()
            ))
        }
    }

    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        unsafe {
            VectorView::from_slice_mut(std::slice::from_raw_parts_mut(
                self.0.as_mut_ptr().add(r * self.cols()),
                self.cols()
            ))
        }
    }

    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        StrideVectorView::from_raw_parts(
            unsafe { self.0.as_ptr().add(c) },
            self.rows(),
            self.cols()
        )
    }

    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        StrideVectorView::from_raw_parts_mut(
            unsafe { self.0.as_mut_ptr().add(c) },
            self.rows(),
            self.cols()
        )
    }
}

pub type TransposedMatrixView<T> = Matrix<T, TransposedMatrixStorage<T>>;

impl<T> TransposedMatrixView<T> {
    /// Returns a view of the matrix transposed.
    pub fn transpose(&self) -> &MatrixView<T> {
        unsafe {
            &*(self as *const _ as *const _)
        }
    }

    /// Returns a mutable view of the matrix transposed.
    pub fn transpose_mut(&mut self) -> &mut MatrixView<T> {
        unsafe {
            &mut *(self as *mut _ as *mut _)
        }
    }
}

/// Very hacky, see [StrideStorage] for explanation.
/// Metadata stores the number of columns (in the least significant 4 bytes)
/// and rows (in the most significant 4 bytes).
/// This allows you to transmute a [SliceMatrixStorage] reference into a
/// [TransposedMatrixStorage] reference while transposing the view.
pub struct TransposedMatrixStorage<T>(
    CustomMetadataSlice<T, TransposedMatrixStorageMetadata>
);

struct TransposedMatrixStorageMetadata {
    cols: usize,
    rows: usize,
}

impl CustomMetadata for TransposedMatrixStorageMetadata {
    fn size(&self) -> usize {
        self.rows * self.cols
    }
}

impl<T> MatrixStorage<T> for TransposedMatrixStorage<T> {
    type RowVecStorage = StrideStorage<T>;
    type ColVecStorage = SliceVectorStorage<T>;
    fn rows(&self) -> usize {
        self.0.metadata().rows
    }

    fn cols(&self) -> usize {
        self.0.metadata().cols
    }

    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        StrideVectorView::from_raw_parts(
            unsafe { self.0.as_ptr().add(r) },
            self.cols(),
            self.rows()
        )
    }

    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        StrideVectorView::from_raw_parts_mut(
            unsafe { self.0.as_mut_ptr().add(r) },
            self.cols(),
            self.rows()
        )
    }

    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        unsafe {
            VectorView::from_slice(std::slice::from_raw_parts(
                self.0.as_ptr().add(c * self.cols()),
                self.rows()
            ))
        }
    }

    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        unsafe {
            VectorView::from_slice_mut(std::slice::from_raw_parts_mut(
                self.0.as_mut_ptr().add(c * self.cols()),
                self.rows()
            ))
        }
    }
}

/// A matrix with a single row.
pub type RowVector<T, S> = Matrix<T, RowVectorStorage<T, S>>;

#[repr(transparent)]
pub struct RowVectorStorage<T, S: VectorStorage<T> + ?Sized>(PhantomData<T>, S);

impl<T, S: VectorStorage<T> + ?Sized> MatrixStorage<T> for RowVectorStorage<T, S> {
    type RowVecStorage = S;
    type ColVecStorage = SliceVectorStorage<T>;

    fn rows(&self) -> usize {
        1
    }

    fn cols(&self) -> usize {
        self.1.dim()
    }

    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        Vector::from_storage_ref(&self.1)
    }

    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        Vector::from_storage_ref_mut(&mut self.1)
    }

    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        VectorView::from_slice(std::slice::from_ref(self.1.entry(c)))
    }

    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        VectorView::from_slice_mut(std::slice::from_mut(self.1.entry_mut(c)))
    }
}

/// A matrix with a single column.
pub type ColumnVector<T, S> = Matrix<T, ColumnVectorStorage<T, S>>;

#[repr(transparent)]
pub struct ColumnVectorStorage<T, S: VectorStorage<T> + ?Sized>(PhantomData<T>, S);

impl<T, S: VectorStorage<T> + ?Sized> MatrixStorage<T> for ColumnVectorStorage<T, S> {
    type RowVecStorage = SliceVectorStorage<T>;
    type ColVecStorage = S;

    fn rows(&self) -> usize {
        self.1.dim()
    }

    fn cols(&self) -> usize {
        1
    }

    fn row(&self, r: usize) -> &Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        VectorView::from_slice(std::slice::from_ref(self.1.entry(r)))
    }

    fn row_mut(&mut self, r: usize) -> &mut Vector<T, Self::RowVecStorage> {
        assert!(r < self.rows());
        VectorView::from_slice_mut(std::slice::from_mut(self.1.entry_mut(r)))
    }

    fn col(&self, c: usize) -> &Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        Vector::from_storage_ref(&self.1)
    }

    fn col_mut(&mut self, c: usize) -> &mut Vector<T, Self::ColVecStorage> {
        assert!(c < self.cols());
        Vector::from_storage_ref_mut(&mut self.1)
    }
}

impl<T, S, R> Mul<&Matrix<T, R>> for &Matrix<T, S>
where
    T: InnerProduct,
    S: MatrixStorage<T> + ?Sized,
    R: MatrixStorage<T> + ?Sized,
{
    type Output = OwnedMatrix<T>;
    fn mul(self, rhs: &Matrix<T, R>) -> Self::Output {
        let r = self.nrows();
        let c = rhs.ncols();
        let iter = iproduct!(0..r, 0..c)
            .map(|(r, c)| self.row(r).dot(rhs.col(c)));
        OwnedMatrix::from_iter(r, c, iter)
    }
}

impl<T, S, R> Mul<&Vector<T, R>> for &Matrix<T, S>
where
    S: MatrixStorage<T> + ?Sized,
    R: VectorStorage<T> + ?Sized,
    for<'a> &'a Matrix<T, S>: Mul<&'a ColumnVector<T, R>, Output = OwnedMatrix<T>>
{
    type Output = OwnedVector<T>;

    fn mul(self, rhs: &Vector<T, R>) -> Self::Output {
        let r = self * rhs.column_vector();
        r.into_vector()
    }
}

impl<T, S, R> Mul<&Matrix<T, R>> for &Vector<T, S>
where
    T: InnerProduct,
    S: VectorStorage<T> + ?Sized,
    R: MatrixStorage<T> + ?Sized,
{
    type Output = OwnedVector<T>;

    fn mul(self, rhs: &Matrix<T, R>) -> Self::Output {
        let r = self.row_vector() * rhs;
        r.into_vector()
    }
}

impl<T: PartialEq, S, R> PartialEq<Matrix<T, R>> for Matrix<T, S>
where
    T: InnerProduct,
    S: MatrixStorage<T> + ?Sized,
    R: MatrixStorage<T> + ?Sized,
{
    fn eq(&self, other: &Matrix<T, R>) -> bool {
        self.nrows() == other.nrows() && self.ncols() == other.ncols()
            && self.entries_row_major().eq(other.entries_row_major())
    }
}

impl<T: Debug, S: MatrixStorage<T> + ?Sized> Debug for Matrix<T, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries(self.rows())
            .finish()
    }
}

impl<T: Clone> Clone for OwnedMatrix<T> {
    fn clone(&self) -> Self {
        OwnedMatrix::from_iter(
            self.nrows(),
            self.ncols(),
            self.entries_row_major().cloned()
        )
    }
}

#[test]
fn transpose_test() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    let t = m.transpose();
    assert_eq!(t.nrows(), 2);
    assert_eq!(t.ncols(), 2);
    assert_eq!(t[0], [2, 4]);
    assert_eq!(t[1], [3, 5]);
    assert_eq!(t[(0, 0)], 2);
    assert_eq!(t[(0, 1)], 4);
    assert_eq!(t[(1, 0)], 3);
    assert_eq!(t[(1, 1)], 5);
}

#[test]
fn row() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    assert_eq!(m.row(0), &[2, 3]);
    assert_eq!(m.row(1), &[4, 5]);
}

#[test]
fn col() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    assert_eq!(m.col(0), &[2, 4]);
    assert_eq!(m.col(1), &[3, 5]);
}

#[test]
fn rows_iter() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    let mut r = m.rows();
    assert_eq!(r.next().unwrap(), &[2, 3]);
    assert_eq!(r.next().unwrap(), &[4, 5]);
    assert!(r.next().is_none());
}

#[test]
fn col_iter() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
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
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    let mut r = m.rows();
    assert_eq!(r.next_back().unwrap(), &[4, 5]);
    assert_eq!(r.next_back().unwrap(), &[2, 3]);
    assert!(r.next_back().is_none());
}

#[test]
fn col_iter_rev() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
    let mut r = m.cols();
    assert_eq!(r.next_back().unwrap(), &[3, 5]);
    assert_eq!(r.next_back().unwrap(), &[2, 4]);
    assert!(r.next_back().is_none());
}

#[test]
fn entry_iterators() {
    let m = OwnedMatrix::<i32>::from_rows(&[
        [2, 3],
        [4, 5],
    ]);
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