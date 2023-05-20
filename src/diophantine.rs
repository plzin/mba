use num_traits::Zero;
use rug::ops::DivRounding;
use rug::{Integer, Complete};
use rand::random;

use crate::matrix::Matrix;
use crate::vector::*;
use crate::lattice::{AffineLattice, Lattice};

/// Computes the (row-style) hermite normal form of a matrix in place
/// and returns the transformation matrix.
pub fn hermite_normal_form(a: &mut Matrix) -> Matrix {
    // The transformation matrix.
    let mut u = Matrix::identity(a.rows);

    let mut r = 0;
    let mut c = 0;
    while r < a.rows && c < a.cols {
        // Choose a pivot in the jth column.
        let pivot = a.column(c)
            .enumerate()
            .skip(r)
            .filter(|e| !e.1.is_zero())
            //.min_by_key(|e| e.1.clone().abs())
            .min_by(|a, b| a.1.cmp_abs(b.1))
            .map(|e| e.0);

        let Some(pivot) = pivot else {
            // If we didn't find a pivot then the column is 0.
            // Continue with the next one.
            c += 1;
            continue;
        };

        // Move the pivot to the beginning.
        a.swap_rows(r, pivot);
        u.swap_rows(r, pivot);

        // Try to eliminate every other entry in the column.
        // This might not work instantly.
        // If there remain non-zero entries in this column,
        // then we will go over this column again.
        for k in r+1..a.rows {
            if a[(k, c)] != 0 {
                let m = -(&a[(k, c)] / &a[(r, c)]).complete();

                a.row_multiply_add(r, k, &m);
                u.row_multiply_add(r, k, &m);
            }
        }

        // If there is any non-zero element then we need to continue in the same column.
        if a.column(c).skip(r + 1).any(|e| *e != 0) {
            continue;
        }

        // Flip sign if necessary.
        if a[(r, c)] < 0 {
            a.flip_sign_row(r);
            u.flip_sign_row(r);
        }

        // Reduce the elements above the pivot
        // (in the column of the pivot and rows above the pivot).
        // The Hermite normal form requires the entries
        // above the pivot to be positive.
        if !a[(r, c)].is_zero() {
            for k in 0..r {
                let entry = &a[(k, c)];
                let m = -Integer::from(entry.div_euc(&a[(r, c)]));

                if m != 0 {
                    a.row_multiply_add(r, k, &m);
                    u.row_multiply_add(r, k, &m);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    return u;
}

/// Solves a system of linear diophantine equations.
pub fn solve(a: &Matrix, b: &IVector) -> AffineLattice {
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
    let rank = m.rows()
        .take_while(|r| r.iter().any(|e| !e.is_zero()))
        .count();

    // Make sure the hermite normal form has the correct form,
    // because only then does it have a solution.
    let r = rank - 1;
    let has_solution = m[(r, m.cols - 1)] == 1
        && m.row(r).iter().take(m.cols - 1).all(|e| *e == 0);

    if !has_solution {
        return AffineLattice::empty();
    }

    let offset = -Vector::from_entries(
        &u.row(r).as_slice()[..u.rows-1]
    );

    let basis = Matrix::from_iter(u.rows - rank, u.rows - 1,
        u.rows().skip(rank).flat_map(|r| r.iter().take(u.rows - 1).cloned()) 
    );

    AffineLattice::from_offset_basis(offset, basis)
}

/// Solves a linear system of equations Ax=b mod n.
pub fn solve_modular(a: &Matrix, b: &IVector, n: &Integer) -> AffineLattice {
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

    let offset = IVector::from_entries(&l.offset.as_slice()[..a.cols])
        .map(|i| *i = i.div_rem_euc_ref(n).complete().1);

    let iter = l.lattice.basis.rows()
        .flat_map(|e| e.iter().take(a.cols).map(|i| i.clone() % n))
        .chain((0..a.cols).flat_map(|i| (0..a.cols).map(
            move |j| if j == i { n.clone() } else { Integer::new() }))
        );
    let mut bm = Matrix::from_iter(l.lattice.rank() + a.cols, a.cols,
        iter
    );

    let lattice = Lattice::from_generating_set(bm);

    AffineLattice {
        offset,
        lattice,
    }
}

#[test]
fn small_test() {
    let a = Matrix::from_rows(&[
        [0, 0, -1, 1],
        [0, 1, 0, 0],
        [0, 1, 0, 0],
        [1, 0, 0, 0],
    ]);

    let b = Vector::from_entries(&[0, 1, 1, 2]);

    let l = solve(&a, &b);
    assert!(l.offset.as_slice() == &[2, 1, 0, 0]);
    assert!(l.lattice.rank() == 1);
    assert!(l.lattice.basis.row(0).as_slice() == &[0, 0, 1, 1]);
}