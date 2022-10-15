use std::marker::PhantomData;
use rug::{Integer, Complete};
use crate::vector::Vector;

pub struct Matrix {
    /// The number of rows.
    pub rows: usize,

    /// The number of columns.
    pub cols: usize,

    /// Memory that holds the entries.
    pub(self) entries: *mut Integer,
}

impl Matrix {
    /// Return an empty matrix.
    pub fn empty() -> Self {
        Self {
            rows: 0,
            cols: 0,
            entries: core::ptr::null_mut(),
        }
    }

    /// Returns an uninitialized mxn matrix.
    pub(self) fn uninit(m: usize, n: usize) -> Self {
        if m == 0 || n == 0 {
            return Self::empty();
        }

        // The memory layout of the elements.
        let layout = std::alloc::Layout::from_size_align(
            m * n * core::mem::size_of::<Integer>(),
            core::mem::align_of::<Integer>()
        ).unwrap();

        // Allocate memory for the entries.
        let entries = unsafe {
            std::alloc::alloc(layout) as *mut Integer
        };

        Self {
            rows: m,
            cols: n,
            entries
        }
    }

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

    /// Creates a matrix from an array.
    pub fn from_array<T: Into<Integer>, const R: usize, const C: usize>(
        a: [[T; C]; R]
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

    /// Creates a matrix from a vector.
    pub fn from_vec<T: Into<Integer>>(v: Vec<Vec<T>>) -> Self {
        let m = v.len();
        if m == 0 {
            return Self::empty();
        }

        let n = v[0].len();
        if n == 0 {
            return Self::empty();
        }

        let r = Self::uninit(m, n);

        let mut ptr = r.entries;
        for row in v {
            assert!(row.len() == n,
                "Can't create matrix from vector \
                    because not every row has the same number of elements");
            for e in row {
                unsafe {
                    ptr.write(e.into());
                    ptr = ptr.add(1);
                }
            }
        }

        r
    }

    /// Returns a slice of a row.
    pub fn row(&self, r: usize) -> &[Integer] {
        debug_assert!(r < self.rows);
        unsafe {
            core::slice::from_raw_parts(
                self.entries.add(r * self.cols),
                self.cols
            )
        }
    }

    /// Returns a mutable slice of a row.
    pub fn row_mut(&mut self, r: usize) -> &mut [Integer] {
        debug_assert!(r < self.rows);
        unsafe {
            core::slice::from_raw_parts_mut(
                self.entries.add(r * self.cols),
                self.cols
            )
        }
    }

    /// Returns an iterator over the column.
    pub fn column(&self, c: usize) -> Column<'_> {
        Column::from_matrix(self, c)
    }

    /// Returns an iterator over the column.
    pub fn column_mut(&mut self, c: usize) -> ColumnMut<'_> {
        ColumnMut::from_matrix(self, c)
    }

    /// Returns an immutable reference to an element.
    pub fn entry(&self, r: usize, c: usize) -> &Integer {
        debug_assert!(r < self.rows && c < self.cols, "Matrix index out of range");

        unsafe {
            &*self.entries.add(r * self.cols + c)
        }
    }

    /// Returns a mutable reference to an element.
    pub fn entry_mut(&mut self, r: usize, c: usize) -> &mut Integer {
        debug_assert!(r < self.rows && c < self.cols, "Matrix index out of range");

        unsafe {
            &mut *self.entries.add(r * self.cols + c)
        }
    }

    /// Returns a pointer to an entry.
    pub fn entry_ptr(&mut self, r: usize, c: usize) -> *mut Integer {
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

    /// Swap to columns.
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

    /// Add a scaled row to another row. N = c * M.
    pub fn col_multiply_add(&mut self, n: usize, m: usize, c: &Integer) {
        debug_assert!(n < self.cols && m < self.cols);
        for i in 0..self.rows {
            let c = (&self[(i, n)] * c).complete();
            self[(i, m)] += c;
        }
    }
}

impl std::ops::Index<(usize, usize)> for Matrix {
    type Output = Integer;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        self.entry(index.0, index.1)
    }

}

impl std::ops::IndexMut<(usize, usize)> for Matrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        self.entry_mut(index.0, index.1)
    }

}

impl std::ops::Mul for &Matrix {
    type Output = Matrix;

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

impl std::ops::Mul<&Vector> for &Matrix {
    type Output = Vector;

    fn mul(self, rhs: &Vector) -> Self::Output {
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

impl std::fmt::Debug for Matrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut l = f.debug_list();

        for r in 0..self.rows {
            l.entry(&self.row(r));
        }

        l.finish()
    }
}

impl Drop for Matrix {
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

pub struct Column<'a> {
    ptr: *mut Integer,
    off: usize,
    last: *mut Integer,
    lifetime: PhantomData<&'a Integer>,
}

impl<'a> Column<'a> {
    pub fn from_matrix(a: &'a Matrix, col: usize) -> Self {
        debug_assert!(col < a.cols);
        Self {
            ptr: unsafe { a.entries.add(col) },
            off: a.cols,
            last: unsafe { a.entries.add(col + (a.rows - 1) * a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a> Iterator for Column<'a> {
    type Item = &'a Integer;
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

pub struct ColumnMut<'a> {
    ptr: *mut Integer,
    off: usize,
    last: *mut Integer,
    lifetime: PhantomData<&'a mut Integer>,
}

impl<'a> ColumnMut<'a> {
    pub fn from_matrix(a: &'a mut Matrix, col: usize) -> Self {
        debug_assert!(col < a.cols);
        Self {
            ptr: unsafe { a.entries.add(col) },
            off: a.cols,
            last: unsafe { a.entries.add(col + (a.rows - 1) * a.cols) },
            lifetime: PhantomData,
        }
    }
}

impl<'a> Iterator for ColumnMut<'a> {
    type Item = &'a mut Integer;
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
