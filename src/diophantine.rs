//! Solves linear diophantine equations.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::Euclid;
use num_traits::One;
use num_traits::Signed;
use num_traits::Zero;

use crate::keep_bits;
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

                a.row_multiply_add(k, r, &m);
                u.row_multiply_add(k, r, &m);
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
                    a.row_multiply_add(k, r, &m);
                    u.row_multiply_add(k, r, &m);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    u
}

/// Computes the (row-style) hermite normal form of a matrix in place
/// and returns the transformation matrix.
pub fn hermite_normal_form_mod(a: &mut IOwnedMatrix, bits: u32) -> IOwnedMatrix {
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

                a.row_multiply_add(k, r, &m);
                u.row_multiply_add(k, r, &m);
                a.row_mut(k).keep_signed_bits(bits);
                u.row_mut(k).keep_signed_bits(bits);
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
                    a.row_multiply_add(k, r, &m);
                    u.row_multiply_add(k, r, &m);
                    a.row_mut(k).keep_signed_bits(bits);
                    u.row_mut(k).keep_signed_bits(bits);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    a.keep_bits(bits);
    u.keep_bits(bits);
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

/// Computes a diagonal matrix D in-place
/// and returns matrices (S, T), such that D=SAT mod n.
pub fn diagonalize(
    a: &mut IMatrixView, bits: u32
) -> (IOwnedMatrix, IOwnedMatrix) {
    // The matrices S and T are initialized to the identity.
    // S/T keeps track of the row/column operations.
    let mut s = IOwnedMatrix::identity(a.nrows());
    let mut t = IOwnedMatrix::identity(a.ncols());

    for i in 0..a.min_dim() {
        //
        // Eliminate row i and column i.
        //
        loop {
            // Is there a non-zero element in the column?
            let col_zero = a.col(i)
                .iter()
                .skip(i+1)
                .all(|e| e.is_zero());

            if !col_zero {
                //
                // Eliminate the column.
                //

                // Find a pivot in the column.
                let pivot = a.col(i)
                    .iter()
                    .enumerate()
                    .skip(i)
                    .filter(|e| !e.1.is_zero())
                    .min_by_key(|e| e.1.magnitude())
                    .map(|e| e.0)
                    .unwrap(); // We know there is a non-zero element.

                // Move the pivot to the beginning.
                a.swap_rows(i, pivot);
                s.swap_rows(i, pivot);

                // Try to eliminate every other entry in the column.
                for k in i+1..a.nrows() {
                    if !a[(k, i)].is_zero() {
                        let m = -(&a[(k, i)] / &a[(i, i)]);
                        a.row_multiply_add(k, i, &m);
                        s.row_multiply_add(k, i, &m);
                        a.row_mut(k).keep_bits(bits);
                        s.row_mut(k).keep_bits(bits);
                    }
                }

                // Keep eliminating the column.
                continue;
            }

            // If we get here, the column is zero.

            // Is there a non-zero element in the row?
            let row_zero = a.row(i)
                .iter()
                .skip(i+1)
                .all(|e| e.is_zero());

            // If the row is zero, then continue with the next row/column.
            if row_zero {
                break;
            }

            //
            // Eliminate the row.
            //

            // Find a pivot in the row.
            let pivot = a.row(i)
                .iter()
                .enumerate()
                .skip(i)
                .filter(|e| !e.1.is_zero())
                .min_by_key(|e| e.1.magnitude())
                .map(|e| e.0)
                .unwrap(); // We know there is a non-zero element.

            // Move the pivot to the beginning.
            a.swap_columns(i, pivot);
            t.swap_columns(i, pivot);

            // Try to eliminate every other entry in the row.
            for k in i+1..a.ncols() {
                if !a[(i, k)].is_zero() {
                    let m = -(&a[(i, k)] / &a[(i, i)]);
                    a.col_multiply_add(k, i, &m);
                    t.col_multiply_add(k, i, &m);
                    a.col_mut(k).keep_bits(bits);
                    t.col_mut(k).keep_bits(bits);
                }
            }
        }
    }

    (s, t)
}

/// Solves a system of linear congruences Ax=b.
pub fn solve_congruences(
    mut a: IOwnedMatrix, b: IOwnedVector, bits: u32,
) -> AffineLattice {
    debug_assert!(a.nrows() == b.dim(), "Invalid system of congruences");

    // Diagonalize the system.
    let (s, t) = diagonalize(a.view_mut(), bits); println!("s: {s:?}\na: {a:?}\nt: {t:?}");

    // We could already do this in diagonalize if we really wanted.
    let mut b = &s * &b;
    b.keep_bits(bits);

    // If there is a non-zero entry in b at index > a.min_dim()
    // then the system has no solution, since the corresponding
    // row in a is zero, so we are solving 0=x.
    if b.iter().skip(a.min_dim()).any(|e| !e.is_zero()) {
        return AffineLattice::empty();
    }

    // Some solution to the system.
    let mut offset = Vector::zero(a.ncols());

    // The basis of the kernel.
    let mut basis = Vec::new();

    // Solve the scalar linear congruences.
    let n = BigInt::one() << bits;
    for i in 0..a.min_dim() {
        let l = keep_bits(&a[(i, i)], bits);
        let r = keep_bits(&b[i], bits);
        let Some((x, kern)) = solve_scalar_congruence(l, r, &n) else {
            // If there is no solution,
            // then the whole system does not have a solution.
            return AffineLattice::empty();
        };

        offset[i] = x;

        if !kern.is_zero() {
            let mut v = Vector::zero(a.ncols());
            v[i] = kern;
            basis.push(v);
        }
    }

    // If there are more variables then equations
    // then there are no restrictions on the variables
    // from index d.rows
    for i in a.nrows()..a.ncols() {
        let mut v = Vector::zero(a.ncols());
        v[i] = One::one();
        basis.push(v);
    }

    offset = &t * &offset;
    for v in &mut basis {
        *v = &t * &*v;
    }

    let lattice = Lattice::from_basis(IOwnedMatrix::from_rows(&basis));

    AffineLattice {
        offset,
        lattice,
    }
}

/// Solves ax=b mod n.
/// Returns None if there is no solution.
/// Otherwise returns all solutions in the form (c, d)
/// where c+di are all solutions.
pub fn solve_scalar_congruence(
    a: BigInt, b: BigInt, n: &BigInt
) -> Option<(BigInt, BigInt)> {
    assert!(!n.is_zero());
    // Handle the case that a is zero, so we don't have to think about it.
    if a.is_zero() {
        return b.is_zero().then_some((Zero::zero(), One::one()));
    }

    let (mut old_r, mut r) = (n.clone(), a);
    let (mut old_t, mut t) = (BigInt::zero(), BigInt::one());
    while !r.is_zero() {
        let q = &old_r / &r;
        let new_r = &old_r - &q * &r;
        let new_t = &old_t - &q * &t;
        (old_r, r) = (r, new_r);
        (old_t, t) = (t, new_t);
    }

    // old_r is gcd(a, n).
    let gcd = old_r;

    // There is a solution iff gcd divides b.
    // old_t is the Bezout coefficient: a*old_t=gcd(a, n) mod n.
    if !b.is_multiple_of(&gcd) {
        return None;
    }
    let x = &b / &gcd * &old_t;
    let x = x.rem_euclid(n);

    if t.is_negative() {
        t = -t;
    }

    Some((x, t))
}

#[test]
fn scalar_congruence_test() {
    let a = BigInt::from(2);
    let b = BigInt::from(8);
    let n = BigInt::from(10);
    let (x, kern) = solve_scalar_congruence(
        a.clone(),
        b.clone(),
        &n
    ).unwrap();
    println!("{x} + {kern}*i");
    assert_eq!((&a * &x).rem_euclid(&n), b.rem_euclid(&n));
    assert!((&a * &kern).rem_euclid(&n).is_zero());
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