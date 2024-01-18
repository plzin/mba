//! Solves linear diophantine equations.

use num_bigint::BigInt;
use num_traits::Euclid;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use crate::matrix::*;
use crate::vector::*;
use crate::lattice::{AffineLattice, Lattice};

/// Computes the (row-style) hermite normal form of a matrix in place
/// and returns the transformation matrix.
pub fn hermite_normal_form(a: &mut IOwnedMatrix) -> IOwnedMatrix {
    // The transformation matrix.
    let mut u = Matrix::identity(a.nrows());

    let mut r = 0;
    let mut c = 0;
    while r < a.nrows() && c < a.ncols() {
        // Choose a pivot in the jth column.
        let pivot = a.col(c)
            .iter()
            .enumerate()
            .skip(r)
            .filter(|e| !e.1.is_zero())
            .min_by(|a, b| a.1.magnitude().cmp(b.1.magnitude()))
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
        for k in r+1..a.nrows() {
            if !a[(k, c)].is_zero() {
                let m = -(&a[(k, c)] / &a[(r, c)]);

                a.row_multiply_add(r, k, &m);
                u.row_multiply_add(r, k, &m);
            }
        }

        // If there is any non-zero element then we need to continue in the same column.
        if a.col(c).iter().skip(r + 1).any(|e| !e.is_zero()) {
            continue;
        }

        // Flip sign if necessary.
        if a[(r, c)].is_negative() {
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
                let m = -entry.div_euclid(&a[(r, c)]);

                if !m.is_zero() {
                    a.row_multiply_add(r, k, &m);
                    u.row_multiply_add(r, k, &m);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    u
}

/// Solves a system of linear diophantine equations.
pub fn solve(a: &IMatrixView, b: &IVectorView) -> AffineLattice {
    assert!(a.nrows() == b.dim(),
        "Vector must have an entry for each row in the matrix.");

    let mut m = Matrix::zero(a.ncols() + 1, a.nrows() + 1);

    // Initialize the matrix m.
    for i in 0..a.nrows() {
        for (j, e) in a.row(i).iter().enumerate() {
            m[(j, i)] = e.clone();
        }
    }

    for i in 0..b.dim() {
        m[(a.ncols(), i)] = b[i].clone();
    }

    m[(a.ncols(), a.nrows())] = One::one();

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
    let has_solution = m[(r, m.ncols() - 1)].is_one()
        && m.row(r).iter().take(m.ncols() - 1).all(|e| e.is_zero());

    if !has_solution {
        return AffineLattice::empty();
    }

    let offset = -IOwnedVector::from_entries(
        &u.row(r).as_slice()[..u.nrows()-1]
    );

    let basis = Matrix::from_iter(u.nrows() - rank, u.nrows() - 1,
        u.rows().skip(rank).flat_map(|r| r.iter().take(u.nrows() - 1).cloned())
    );

    AffineLattice::from_offset_basis(offset, basis)
}

/// Solves a linear system of equations Ax=b mod n.
/// The solution lattice consists of all integer solutions to the equations.
pub fn solve_modular(
    a: &IMatrixView,
    b: &IVectorView,
    n: &BigInt
) -> AffineLattice {
    //
    // Concatenate an n times the identity matrix to the right of A.
    //
    let mut m = Matrix::zero(a.nrows(), a.ncols() + a.nrows());

    // Copy the old matrix.
    for i in 0..a.nrows() {
        for j in 0..a.ncols() {
            m[(i, j)] = a[(i, j)].clone();
        }
    }

    // Set the identity matrix.
    for i in 0..a.nrows() {
        m[(i, a.ncols() + i)] = n.clone();
    }

    // Solve the diophantine system.
    let l = solve(m.view(), b);
    if l.is_empty() {
        return l;
    }

    // Clean up the solution by taking everything mod n
    // removing the last components that correspond to the multipliers
    // of the n's and then removing (now) linearly dependent basis vectors.

    let offset = IOwnedVector::from_entries(&l.offset.as_slice()[..a.ncols()])
        .map(|i| i.rem_euclid(n));

    // This might be the worst code in the history of code.
    let iter = l.lattice.basis.rows()
        .flat_map(|e| e.iter().take(a.ncols()).map(|i| i.clone() % n))
        .chain((0..a.ncols()).flat_map(|i| (0..a.ncols()).map(
            move |j| if j == i { n.clone() } else { Zero::zero() }))
        );
    let bm = Matrix::from_iter(l.lattice.rank() + a.ncols(), a.ncols(),
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

    let b = Vector::from_entries([0, 1, 1, 2]);

    let l = solve(a.view(), b.view());
    assert_eq!(l.offset.as_slice(),
        [2.into(), 1.into(), 0.into(), 0.into()]);
    assert_eq!(l.lattice.rank(), 1);
    assert_eq!(l.lattice.basis.row(0).as_slice(),
        [0.into(), 0.into(), 1.into(), 1.into()]);
}