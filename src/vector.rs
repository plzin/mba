use rug::{Integer, Complete};

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

    /// Creates a vector from a slice.
    pub fn from_slice(a: &[Integer]) -> Self {
        let v = Self::uninit(a.len());
        for i in 0..v.dim {
            unsafe {
                v.entries.add(i).write(a[i].clone());
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

    /// Returns a slice of the vector entries.
    pub fn as_slice(&self) -> &[Integer] {
        unsafe {
            core::slice::from_raw_parts(self.entries, self.dim)
        }
    }

    /// Returns a mutable slice of the vector entries.
    pub fn as_slice_mut(&mut self) -> &mut [Integer] {
        unsafe {
            core::slice::from_raw_parts_mut(self.entries, self.dim)
        }
    }

    /// Returns an immutable reference to an entry.
    pub fn entry(&self, i: usize) -> &Integer {
        unsafe {
            &*self.entries.add(i)
        }
    }

    /// Returns a mutable reference to an entry.
    pub fn entry_mut(&mut self, i: usize) -> &mut Integer {
        unsafe {
            &mut *self.entries.add(i)
        }
    }

    /// Swap the two entries i and j.
    pub fn swap(&mut self, i : usize, j: usize) {
        unsafe {
            core::ptr::swap(self.entries.add(i), self.entries.add(j));
        }
    }

    /// Returns an iterator over the elements.
    pub fn iter(&self) -> std::slice::Iter<'_, Integer> {
        self.as_slice().iter()
    }

    /// Returns an iterator over the mutable elements.
    pub fn iter_mut(&mut self) -> std::slice::IterMut<'_, Integer> {
        self.as_slice_mut().iter_mut()
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

impl std::ops::Add for &Vector {
    type Output = Vector;

    fn add(self, rhs: Self) -> Self::Output {
        assert!(self.dim == rhs.dim, "Can not add vectors of incompatible sizes");
        let mut v = Vector::zero(self.dim);

        for i in 0..v.dim {
            v[i] = (&self[i] + &rhs[i]).complete();
        }

        return v;
    }
}

impl std::ops::AddAssign for Vector {
    fn add_assign(&mut self, mut rhs: Self) {
        assert!(self.dim == rhs.dim, "Can not add vectors of incompatible sizes");

        for i in 0..self.dim {
            let r = unsafe {
                core::ptr::read(&rhs[i])
            };

            self[i] += r;
        }

        unsafe { rhs.free(); }
        core::mem::forget(rhs);
    }
}

impl std::ops::Mul<&Integer> for &Vector {
    type Output = Vector;

    fn mul(self, rhs: &Integer) -> Self::Output {
        let mut v = Vector::zero(self.dim);

        for i in 0..v.dim {
            v[i] = (&self[i] * rhs).complete();
        }

        return v;
    }
}

impl std::ops::Rem<&Integer> for Vector {
    type Output = Vector;
    fn rem(mut self, rhs: &Integer) -> Self::Output {
        for e in self.iter_mut() {
            *e %= rhs;
        }
        self
    }
}

impl std::ops::RemAssign<&Integer> for Vector {
    fn rem_assign(&mut self, rhs: &Integer) {
        for e in self.iter_mut() {
            *e %= rhs;
        }
    }
}

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
