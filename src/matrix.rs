//! Owned matrix with constant runtime dimension.

use std::{marker::PhantomData, fmt::Debug};
use rug::{Integer, Complete, Float};
use crate::vector::*;

/// Matrix.
pub struct Matrix<T> {
    /// Memory that holds the entries.
    pub(self) entries: *mut T,

    /// The number of rows.
    pub rows: usize,

    /// The number of columns.
    pub cols: usize,
}

pub type IMatrix = Matrix<Integer>;
pub type FMatrix = Matrix<Float>;

impl<T> Matrix<T> {
    /// Return an empty matrix.
    pub fn empty() -> Self {
        Self {
            rows: 0,
            cols: 0,
            entries: core::ptr::null_mut(),
        }
    }

    /// Is the matrix empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_null()
    }

    /// Returns an uninitialized mxn matrix.
    pub(self) fn uninit(m: usize, n: usize) -> Self {
        if m == 0 || n == 0 {
            return Self::empty();
        }

        // The memory layout of the elements.
        let layout = std::alloc::Layout::from_size_align(
            m * n * core::mem::size_of::<T>(),
            core::mem::align_of::<T>()
        ).unwrap();

        // Allocate memory for the entries.
        let entries = unsafe {
            std::alloc::alloc(layout) as *mut T
        };

        Self {
            rows: m,
            cols: n,
            entries
        }
    }

    /// Creates a matrix from an array of rows.
    pub fn from_array<U: Into<T>, const R: usize, const C: usize>(
        a: [[U; C]; R]
    ) -> Self {
        let m = Self::uninit(R, C);

        // Copy the entries.
        let mut ptr = m.entries;
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

    /// Creates a matrix from an iterator.
    pub fn from_iter<I: Iterator<Item = T>>(r: usize, c: usize, mut iter: I) -> Self {
        let mut m = Self::uninit(r, c);
        for i in 0..r*c {
            let e = iter.next()
                .expect("The iterator needs to return at least r * c items.");
            unsafe {
                m.entries.add(i).write(e);
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

    /// Returns an iterator over the rows as Integer slices.
    pub fn rows(&self) -> RowIter<T> {
        RowIter::from_matrix(self)
    }

    /// Returns an iterator over the rows as mutable Integer slices.
    pub fn rows_mut(&mut self) -> RowIterMut<T> {
        RowIterMut::from_matrix(self)
    }

    /// Returns a slice of a row.
    pub fn row(&self, r: usize) -> &VV<T> {
        debug_assert!(r < self.rows);
        unsafe {
            VV::from_slice(core::slice::from_raw_parts(
                self.entries.add(r * self.cols),
                self.cols
            ))
        }
    }

    /// Returns a mutable slice of a row.
    pub fn row_mut(&mut self, r: usize) -> &mut VV<T> {
        debug_assert!(r < self.rows);
        unsafe {
            VV::from_slice_mut(core::slice::from_raw_parts_mut(
                self.entries.add(r * self.cols),
                self.cols
            ))
        }
    }

    /// Returns an iterator over the column.
    pub fn column(&self, c: usize) -> Column<T> {
        Column::from_matrix(self, c)
    }

    /// Returns an iterator over the column.
    pub fn column_mut(&mut self, c: usize) -> ColumnMut<T> {
        ColumnMut::from_matrix(self, c)
    }

    /// Returns an immutable reference to an element.
    pub fn entry(&self, r: usize, c: usize) -> &T {
        debug_assert!(r < self.rows && c < self.cols, "Matrix index out of range");

        unsafe {
            &*self.entries.add(r * self.cols + c)
        }
    }

    /// Returns a mutable reference to an element.
    pub fn entry_mut(&mut self, r: usize, c: usize) -> &mut T {
        debug_assert!(r < self.rows && c < self.cols, "Matrix index out of range");

        unsafe {
            &mut *self.entries.add(r * self.cols + c)
        }
    }

    /// Returns a pointer to an entry.
    pub fn entry_ptr(&mut self, r: usize, c: usize) -> *mut T {
        debug_assert!(r < self.rows && c < self.cols, "Matrix index out of range");
        unsafe {
            self.entries.add(r * self.cols + c)
        }
    }

    /// Swap two rows.
    pub fn swap_rows(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        unsafe {
            core::ptr::swap_nonoverlapping(
                self.entry_ptr(i, 0),
                self.entry_ptr(j, 0),
                self.cols
            )
        }
    }

    /// Swap two columns.
    pub fn swap_columns(&mut self, i: usize, j: usize) {
        if i == j {
            return;
        }

        for k in 0..self.cols {
            unsafe {
                core::ptr::swap_nonoverlapping(
                    self.entry_ptr(i, k), self.entry_ptr(j, k), 1
                );
            }
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

        let esize = core::mem::size_of::<Integer>(); 
        let layout = std::alloc::Layout::from_size_align(
            self.rows * self.cols * esize,
            core::mem::align_of::<Integer>()
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

impl Matrix<Integer> {
    /// Returns an mxn zero matrix.
    pub fn zero(m: usize, n: usize) -> Self {
        let r = Self::uninit(m, n);
        for i in 0..m*n {
            unsafe {
                r.entries.add(i).write(Integer::new());
            }
        }
        r
    }

    /// Returns an nxn identity matrix.
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zero(n, n);
        for i in 0..n {
            m[(i, i)] = Integer::from(1);
        }
        m
    }

    /// Flip the sign of a row.
    pub fn flip_sign_row(&mut self, row: usize) {
        for e in self.row_mut(row) {
            *e *= -1;
        }
    }

    /// Flip the sign of a column.
    pub fn flip_sign_column(&mut self, col: usize) {
        for e in self.column_mut(col) {
            *e *= -1;
        }
    }

    /// Add a scaled row to another row. N = c * M.
    pub fn row_multiply_add(&mut self, n: usize, m: usize, c: &Integer) {
        debug_assert!(n < self.rows && m < self.rows);
        for i in 0..self.cols {
            let c = (&self[(n, i)] * c).complete();
            self[(m, i)] += c;
        }
    }

    /// Add a scaled column to another column. N = c * M.
    pub fn col_multiply_add(&mut self, n: usize, m: usize, c: &Integer) {
        debug_assert!(n < self.cols && m < self.cols);
        for i in 0..self.rows {
            let c = (&self[(i, n)] * c).complete();
            self[(i, m)] += c;
        }
    }
}

impl<T> std::ops::Index<usize> for Matrix<T> {
    type Output = VV<T>;

    fn index(&self, index: usize) -> &Self::Output {
        self.row(index)
    }
}

impl<T> std::ops::IndexMut<usize> for Matrix<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        self.row_mut(index)
    }
}

impl<T> std::ops::Index<(usize, usize)> for Matrix<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        self.entry(index.0, index.1)
    }

}

impl<T> std::ops::IndexMut<(usize, usize)> for Matrix<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        self.entry_mut(index.0, index.1)
    }
}

impl std::ops::Mul for &IMatrix {
    type Output = IMatrix;

    fn mul(self, rhs: Self) -> Self::Output {
        assert!(self.cols == rhs.rows, "Can't multiply matrices because of incompatible dimensions");

        let mut m = Matrix::zero(self.rows, rhs.cols);

        for i in 0..m.rows {
            for j in 0..m.cols {
                m[(i, j)] = self.row(i).iter()
                    .zip(rhs.column(i))
                    .map(|(l, r)| (l * r).complete())
                    .sum();
            }
        }

        m
    }
}

impl std::ops::Mul<&IVector> for &IMatrix {
    type Output = IVector;

    fn mul(self, rhs: &IVector) -> Self::Output {
        assert!(self.cols == rhs.dim, "Can't multiply matrix and vector because of incompatible dimensions");

        let mut v = Vector::zero(self.rows);

        for i in 0..self.rows {
            v[i] = self.row(i).iter()
                .zip(rhs.iter())
                .map(|(l, r)| (l * r).complete())
                .sum();
        }

        v
    }

}

impl<T: Debug> Debug for Matrix<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_list()
            .entries(self.rows().map(|r| r.as_slice()))
            .finish()
    }
}

impl<T> Drop for Matrix<T> {
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
                self.rows * self.cols * core::mem::size_of::<Integer>(),
                core::mem::align_of::<Integer>()
            ).unwrap();

            // Free the memory.
            std::alloc::dealloc(self.entries as _, layout);
        }
    }
}

pub struct RowIter<'a, T> {
    /// Front iterator.
    front: *mut T,

    /// Back iterator.
    back: *mut T,

    /// The number of elements in a column.
    count: usize,
    lifetime: PhantomData<&'a T>,
}

impl<'a, T> RowIter<'a, T> {
    pub fn from_matrix(a: &'a Matrix<T>) -> Self {
        Self {
            front: a.entries,
            count: a.cols,
            back: unsafe { a.entries.add(a.cols * a.rows).sub(a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a, T> Iterator for RowIter<'a, T> {
    type Item = &'a VV<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.front as usize > self.back as usize {
            return None;
        }

        unsafe {
            let row = std::slice::from_raw_parts(self.front, self.count);
            self.front = self.front.add(self.count);
            Some(VV::from_slice(row))
        }
    }
}

impl<'a, T> DoubleEndedIterator for RowIter<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.front as usize > self.back as usize {
            return None;
        }

        unsafe {
            let row = std::slice::from_raw_parts(self.back, self.count);
            self.back = self.back.sub(self.count);
            Some(VV::from_slice(row))
        }
    }
}

pub struct RowIterMut<'a, T> {
    /// Front iterator.
    front: *mut T,

    /// Back iterator.
    back: *mut T,

    /// The number of elements in a column.
    count: usize,
    lifetime: PhantomData<&'a T>,
}

impl<'a, T> RowIterMut<'a, T> {
    pub fn from_matrix(a: &'a mut Matrix<T>) -> Self {
        Self {
            front: a.entries,
            count: a.cols,
            back: unsafe { a.entries.add(a.cols * a.rows).sub(a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a, T> Iterator for RowIterMut<'a, T> {
    type Item = &'a mut VV<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.front as usize > self.back as usize {
            return None;
        }

        unsafe {
            let row = std::slice::from_raw_parts_mut(self.front, self.count);
            self.front = self.front.add(self.count);
            Some(VV::from_slice_mut(row))
        }
    }
}

impl<'a, T> DoubleEndedIterator for RowIterMut<'a, T> {
    fn next_back(&mut self) -> Option<Self::Item> {
        if self.front as usize > self.back as usize {
            return None;
        }

        unsafe {
            let row = std::slice::from_raw_parts_mut(self.back, self.count);
            self.back = self.back.sub(self.count);
            Some(VV::from_slice_mut(row))
        }
    }
}

pub struct Column<'a, T> {
    ptr: *mut T,
    off: usize,
    last: *mut T,
    lifetime: PhantomData<&'a T>,
}

impl<'a, T> Column<'a, T> {
    pub fn from_matrix(a: &'a Matrix<T>, col: usize) -> Self {
        debug_assert!(col < a.cols);
        Self {
            ptr: unsafe { a.entries.add(col) },
            off: a.cols,
            last: unsafe { a.entries.add(col + (a.rows - 1) * a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a, T> Iterator for Column<'a, T> {
    type Item = &'a T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr > self.last {
            None
        } else {
            unsafe {
                let v = &*self.ptr;
                self.ptr = self.ptr.add(self.off);
                Some(v)
            }
        }
    }
}

pub struct ColumnMut<'a, T> {
    ptr: *mut T,
    off: usize,
    last: *mut T,
    lifetime: PhantomData<&'a mut T>,
}

impl<'a, T> ColumnMut<'a, T> {
    pub fn from_matrix(a: &'a mut Matrix<T>, col: usize) -> Self {
        debug_assert!(col < a.cols);
        Self {
            ptr: unsafe { a.entries.add(col) },
            off: a.cols,
            last: unsafe { a.entries.add(col + (a.rows - 1) * a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a, T> Iterator for ColumnMut<'a, T> {
    type Item = &'a mut T;
    fn next(&mut self) -> Option<Self::Item> {
        if self.ptr > self.last {
            None
        } else {
            unsafe {
                let v = &mut *self.ptr;
                self.ptr = self.ptr.add(self.off);
                Some(v)
            }
        }
    }
}
