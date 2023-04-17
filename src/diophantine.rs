use rug::{Integer, Complete};
use rand::random;

use crate::matrix::Matrix;
use crate::vector::Vector;

/// Describes the lattice of solutions.
#[derive(Debug)]
pub struct AffineLattice {
    pub offset: Vector,
    pub basis: Vec<Vector>,
}

impl AffineLattice {
    /// Creates an empty lattice.
    pub fn empty() -> Self {
        Self {
            offset: Vector::empty(),
            basis: Vec::new(),
        }
    }

    /// Is this lattice empty?
    pub fn is_empty(&self) -> bool {
        self.offset.dim == 0
    }

    /// Returns a random point on the lattice.
    pub fn sample_point(&self) -> Vector {
        assert!(!self.is_empty(), "Lattice is empty.");
        let mut s = self.offset.clone();

        for b in &self.basis {
            let f = Integer::from(random::<usize>());
            s += b * &f;
        }

        s
    }
}


/// Computes the hermite normal form of a matrix in place
/// and returns the transformation matrix.
fn hermite_normal_form(a: &mut Matrix) -> Matrix {
    // The transformation matrix.
    let mut u = Matrix::identity(a.rows);

    let mut i = 0;
    let mut j = 0;
    while i < a.rows && j < a.cols {
        // Choose a pivot in the jth column.
        let pivot = a.column(j)
            .enumerate()
            .skip(i)
            .filter(|e| *e.1 != 0)
            .min_by_key(|e| e.1.clone().abs())
            .map(|e| e.0);

        let Some(pivot) = pivot else {
            // If we didn't find a pivot then the column is 0.
            // Continue with the next one.
            j += 1;
            continue;
        };

        // Move the pivot to the beginning.
        a.swap_rows(i, pivot);
        u.swap_rows(i, pivot);

        // Try to eliminate every other entry in the column.
        // This might not work instantly.
        // If there remain non-zero entries in this column,
        // then we will go over this column again.
        for k in i+1..a.rows {
            if a[(k, j)] != 0 {
                let m = -(&a[(k, j)] / &a[(i, j)]).complete();

                a.row_multiply_add(i, k, &m);
                u.row_multiply_add(i, k, &m);
            }
        }

        // If there is any non-zero element then we need to continue in the same column.
        if a.column(j).skip(i + 1).any(|e| *e != 0) {
            continue;
        }

        // Flip sign if necessary.
        if a[(i, j)] < 0 {
            a.flip_sign_row(i);
            u.flip_sign_row(i);
        }

        // Reduce the elements above the pivot
        // (in the column of the pivot and rows above the pivot).
        if a[(i, j)] != 0 {
            for k in 0..i {
                let entry = &a[(k, j)];
                // The Hermite normal form requires the entries above the pivot to be positive,
                // so if the entry is negative we add 1 to the result of the euclidean division.
                // E.g. -3/2 = -1 => m = 2 so we add -3+2*2 = 1.
                // This isn't really needed for mixed-boolean arithmetic.
                let m = match entry.cmp0() {
                    std::cmp::Ordering::Equal => continue,
                    std::cmp::Ordering::Greater => -(entry / &a[(i, j)]).complete(),
                    std::cmp::Ordering::Less => -(entry / &a[(i, j)]).complete() + 1,
                };

                if m != 0 {
                    a.row_multiply_add(i, k, &m);
                    u.row_multiply_add(i, k, &m);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        j += 1;
        i += 1;
    }

    return u;
}

/// Solves a system of linear diophantine equations.
pub fn solve(a: &Matrix, b: &Vector) -> AffineLattice {
    assert!(a.rows == b.dim,
        "Vector must have an entry for each row in the matrix.");

    let mut m = Matrix::zero(a.cols + 1, a.rows + 1);

    // Initialize the matrix m.
    for i in 0..a.rows {
        for (j, e) in a.row(i).iter().enumerate() {
            m[(j, i)] = e.clone();
        }
    }

    for i in 0..b.dim {
        m[(a.cols, i)] = b[i].clone();
    }

    m[(a.cols, a.rows)] = Integer::from(1);

    // Transform it into hermite normal form.
    let u = hermite_normal_form(&mut m);

    // Compute the rank of the matrix.
    // It has a special form that we can take advantage of.
    let rank = (0..m.rows)
        .position(|i| m.row(i).iter().all(|e| *e == 0))
        .unwrap_or_else(|| std::cmp::min(m.rows, m.cols));

    // Make sure the hermite normal form has the correct form,
    // because only then does it have a solution.
    let r = rank - 1;
    let has_solution = m[(r, m.cols - 1)] == 1
        && m.row(r).iter().take(m.cols - 1).all(|e| *e == 0);

    if !has_solution {
        return AffineLattice::empty();
    }

    let mut offset = Vector::from_slice(&u.row(rank - 1)[..u.rows-1]);
    for i in 0..offset.dim {
        offset[i] *= -1;
    }

    let basis: Vec<_> = (rank..u.rows).map(|i|
        Vector::from_slice(&u.row(i)[..u.rows - 1])).collect();

    AffineLattice {
        offset,
        basis,
    }
}

/// Solves a linear system of equations Ax=b mod n.
pub fn solve_modular(a: &Matrix, b: &Vector, n: &Integer) -> AffineLattice {
    //
    // Concatenate an n times the identity matrix to the right of A.
    //
    let mut m = Matrix::zero(a.rows, a.cols + a.rows);

    // Copy the old matrix.
    for i in 0..a.rows {
        for j in 0..a.cols {
            m[(i, j)] = a[(i, j)].clone();
        }
    }

    // Set the identity matrix.
    for i in 0..a.rows {
        m[(i, a.cols + i)] = n.clone();
    }

    // Solve the diophantine system.
    let l = solve(&m, b);
    if l.is_empty() {
        return l;
    }

    // Clean up the solution by taking everything mod n
    // removing the last components that correspond to the multipliers
    // of the n's and then removing (now) linearly dependent basis vectors.

    let offset = Vector::from_slice(&l.offset.as_slice()[..a.cols]) % n;

    let basis = l.basis.iter()
        .map(|e| Vector::from_slice(&e.as_slice()[..a.cols]) % n)
        .filter(|e| e.iter().any(|e| *e != 0))
        .collect();

    AffineLattice {
        offset,
        basis,
    }
}
