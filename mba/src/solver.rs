//! Solves systems of linear equations modulo some number, usually a power of
//! two.

use num_bigint::BigInt;
use num_bigint::BigUint;
use num_traits::Euclid;

use crate::matrix::*;
use crate::rings::Ring as _;
use crate::rings::RingElement;
use crate::rings::IntDivRing;
use crate::rings::Z;
use crate::vector::*;
use crate::lattice::{AffineLattice, Lattice};

/// Computes the (row-style) Hermite normal form of a matrix in place and
/// returns the transformation matrix.
///
/// This is basically only really well defined for the integers ([`Z`]).
///
/// You can use it for the other rings but it is not a normal form in that case,
/// i.e. there could be multiple matrices that fulfill the conditions of a HNF
/// for a given matrix. (The conditions include the sign of elements which is
/// not really well defined for the other rings.)
pub fn hermite_normal_form<R: IntDivRing>(
    a: &mut OwnedMatrix<R>,
    ring: &R,
) -> OwnedMatrix<R> {
    // The transformation matrix.
    let mut u = Matrix::identity(a.num_rows());

    let mut r = 0;
    let mut c = 0;
    while r < a.num_rows() && c < a.num_cols() {
        // Choose a pivot in the jth column.
        let pivot = a.col(c)
            .iter()
            .enumerate()
            .skip(r)
            .filter(|e| !e.1.is_zero())
            .min_by(|a, b| ring.cmp_abs(a.1, b.1))
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
        for k in r+1..a.num_rows() {
            if a[(k, c)].is_zero() {
                continue;
            }

            //let m = ring.neg(R::rounded_div(&a[(k, c)], &a[(r, c)]));
            let m = ring.neg(R::euclidean_div(&a[(k, c)], &a[(r, c)]));

            a.row_multiply_add(k, r, &m, ring);
            u.row_multiply_add(k, r, &m, ring);
        }

        // If there is any non-zero element then we need to continue in the same column.
        if a.col(c).iter().skip(r + 1).any(|e| !e.is_zero()) {
            continue;
        }

        // Flip sign if necessary.
        if ring.is_negative(&a[(r, c)]) {
            a.negate_row(r, ring);
            u.negate_row(r, ring);
        }

        // Reduce the elements above the pivot
        // (in the column of the pivot and rows above the pivot).
        // The Hermite normal form requires the entries
        // above the pivot to be positive.
        if !a[(r, c)].is_zero() {
            for k in 0..r {
                let entry = &a[(k, c)];
                let m = ring.neg(R::euclidean_div(entry, &a[(r, c)]));

                if !m.is_zero() {
                    a.row_multiply_add(k, r, &m, ring);
                    u.row_multiply_add(k, r, &m, ring);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    u
}

/// Computes a matrix that looks similar to the Hermite normal form, but it is
/// not a normal form. This hermite form is basically an upper triangular matrix
/// such that the elements above the diagonal are smaller in absolute value than
/// the diagonal element.
pub fn modular_hermite_form<R: IntDivRing>(
    a: &mut OwnedMatrix<R>,
    ring: &R,
) -> OwnedMatrix<R> {
    // The transformation matrix.
    let mut u = Matrix::identity(a.num_rows());

    let mut r = 0;
    let mut c = 0;
    'outer: while r < a.num_rows() && c < a.num_cols() {
        // We can rank three kinds of pivots from best to worst:
        //
        // (0. The whole (part of the) column is 0 and the pivot is 0, in which
        //    case the whole (part of) the column is already 0, so we can move
        //    on to the next column.)
        // 1. The pivot is 1. In this case we can directly eliminate the other
        //    elements to 0 without having to compute the inverse or doing
        //    division.
        // 2. The pivot is a unit. In this case we can compute the inverse and
        //    multiply the row by it to make the pivot 1 in which case we can
        //    directly eliminate the other elements to 0 without division.
        // 3. The pivot is not a unit. In this case we have to do the division
        //    to reduce the other elements. We can't necessarily directly
        //    eliminate them.
        //
        // In the case of matrices over `Z/2^{bits}Z` where the elements are
        // sampled uniformly at random and the part of the column we have to
        // consider has size `n`, the (approximate) probability to find each
        // case on the first iteration:
        // 0. `(1/2)^(bits*n)` (i.e. very small)
        // 1. `1 - (1 - (1/2)^bits)^n` (i.e. more `bits` -> less likely, bigger
        //    `n` -> more likely, but overall small)
        // 2. `1 - (1/2)^n` (i.e. very likely)
        // 3. `(1/2)^n` (i.e. small)
        //
        // Note that the probabilities for case 2 and 3 are approximate because
        // they overlap with cases 0 and 1, i.e. you would have to subtract 0
        // from 2 and 1 from 3 to get the exact probabilities. In general the
        // probabilities for case 0 and 1 are small, but for specific cases this
        // is not true.
        //
        // Also note that if we find case 0, 1, or 2, then we will eliminate the
        // (part of) the column in this iteration, so the only case where we
        // need multiple iterations is case 3 in which case we will only ever
        // find cases 0 and 3 in the subsequent iterations.
        //
        // The point is, it is very likely that we will find a unit in these
        // rings. Which makes it worth it for this algorithm to exist in the
        // first place.
        //
        // It is less clear wether it is worth it to check for case 1, but in
        //
        // It might be worth considering if it is worth to check if the pivot is
        // -1 as well, because the inverse computation is trivial in that case.

        let (pivot, is_one) = 'b: {
            // Do we have an element that is 1?
            if let Some((pivot, _)) = a.col(c)
                .iter()
                .enumerate()
                .skip(r)
                .find(|(_, e)| e.is_one())
            {
                break 'b (pivot, true);
            }

            // Do we have an element that is a unit?
            if let Some((pivot, inv)) = a.col(c)
                .iter()
                .enumerate()
                .skip(r)
                .find_map(|(i, e)| ring.inverse(e).map(|inv| (i, inv)))
            {
                // Multiply the row by the inverse.
                a.row_multiply(pivot, &inv, ring);
                u.row_multiply(pivot, &inv, ring);
                break 'b (pivot, true);
            }

            // Find the smallest non-zero element by absolute value.
            if let Some((pivot, _)) = a.col(c)
                .iter()
                .enumerate()
                .skip(r)
                .filter(|e| !e.1.is_zero())
                .min_by(|a, b| ring.cmp_abs(a.1, b.1))
            {
                break 'b (pivot, false);
            }


            // The column is 0, continue with the next one.
            c += 1;
            continue 'outer;
        };

        // Move the pivot to the beginning.
        a.swap_rows(r, pivot);
        u.swap_rows(r, pivot);

        // Try to eliminate every other entry in the column.
        // This might not work instantly.
        // If there remain non-zero entries in this column,
        // then we will go over this column again.
        for k in r+1..a.num_rows() {
            let e = &a[(k, c)];
            if e.is_zero() {
                continue;
            }

            let m = if is_one {
                ring.neg(e.clone())
            } else {
                ring.neg(R::euclidean_div(e, &a[(r, c)]))
            };

            a.row_multiply_add(k, r, &m, ring);
            u.row_multiply_add(k, r, &m, ring);
        }

        // If there is any non-zero element then we need to continue in the
        // same column.
        if !is_one && a.col(c).iter().skip(r + 1).any(|e| !e.is_zero()) {
            continue;
        }

        // Reduce the elements above the pivot
        // (in the column of the pivot and rows above the pivot).
        if !a[(r, c)].is_zero() {
            for k in 0..r {
                let e = &a[(k, c)];
                if e.is_zero() {
                    continue;
                }

                let m = if is_one {
                    ring.neg(e.clone())
                } else {
                    ring.neg(R::euclidean_div(e, &a[(r, c)]))
                };

                if !m.is_zero() {
                    a.row_multiply_add(k, r, &m, ring);
                    u.row_multiply_add(k, r, &m, ring);
                }
            }
        }

        // Continue with the bottom right part of the matrix that remains.
        c += 1;
        r += 1;
    }

    u
}

/// Constructs the matrix whose HNF we can compute to solve the system of
/// linear equations Ax=b.
fn construct_hnf_matrix<R: IntDivRing>(
    a: &MatrixView<R>,
    b: &VectorView<R>,
) -> OwnedMatrix<R> {
    assert_eq!(a.num_rows(), b.dim(),
        "Vector must have an entry for each row in the matrix.");

    let mut m = Matrix::zero(a.num_cols() + 1, a.num_rows() + 1);

    // Copy the transposed matrix into the upper left part.
    for i in 0..a.num_rows() {
        for (j, e) in a.row(i).iter().enumerate() {
            m[(j, i)] = e.clone();
        }
    }

    // Copy the vector into the lower left part.
    for i in 0..b.dim() {
        m[(a.num_cols(), i)] = b[i].clone();
    }

    // Set the lower right entry to 1.
    m[(a.num_cols(), a.num_rows())] = R::one();

    m
}

/// Solves a system of linear equations by constructing a certain matrix and
/// computing the Hermite normal form of it.
///
/// I wrote this with the intention of it being used for rings other than [`Z`]
/// but even after changing things (e.g. allowing the bottom right entry to be
/// a unit) it still did not work in general, so it can only correctly solve
/// systems of linear equations over [`Z`].
fn solve_via_integer_hnf<R: IntDivRing>(
    a: &MatrixView<R>,
    b: &VectorView<R>,
    ring: &R,
) -> AffineLattice<R> {
    let mut m = construct_hnf_matrix(a, b);

    // Transform it into Hermite normal form.
    let u = hermite_normal_form(&mut m, ring);

    // Compute the rank of the matrix.
    // It has a special form that we can take advantage of.
    let rank = m.rows()
        .take_while(|r| r.iter().any(|e| !e.is_zero()))
        .count();

    // Make sure the Hermite normal form has the correct form,
    // because only then does it have a solution.
    let r = rank - 1;

    let has_solution =
        m.row(r).iter().take(m.num_cols() - 1).all(|e| e.is_zero()) &&
        m[(r, m.num_cols() - 1)].is_one();

    if !has_solution {
        return AffineLattice::empty(a.num_cols());
    }

    let offset = OwnedVector::from_entries(
        &u.row(r)[0..u.num_rows()-1]
    ).neg(ring);

    let basis = Matrix::from_iter(u.num_rows() - rank, u.num_rows() - 1,
        u.rows().skip(rank).flat_map(|r| r.iter().take(u.num_rows() - 1).cloned())
    );

    AffineLattice::from_offset_basis(offset, basis)
}

/// Solves a system of linear equations Ax=b mod n by translating it into a
/// system of linear diophantine equations which can be solved by the standard
/// integer HNF algorithm.
///
/// You should probably use [`solve_via_modular_diagonalize`] or
/// [`solve_via_integer_diagonalize`] instead. [`solve_via_modular_diagonalize`]
/// is faster for (at least) binary rings.
pub fn solve_modular_via_integer_hnf(
    a: &MatrixView<Z>,
    b: &VectorView<Z>,
    n: &BigUint,
) -> AffineLattice<Z> {
    let n = BigInt::from(n.clone());
    //
    // Concatenate an n times the identity matrix to the right of A.
    //
    let mut m = Matrix::zero(a.num_rows(), a.num_cols() + a.num_rows());

    // Copy the old matrix.
    for i in 0..a.num_rows() {
        for j in 0..a.num_cols() {
            m[(i, j)] = a[(i, j)].clone();
        }
    }

    // Set the identity matrix.
    for i in 0..a.num_rows() {
        m[(i, a.num_cols() + i)] = n.clone();
    }

    // Solve the diophantine system.
    let l = solve_via_integer_hnf(m.view(), b, &Z);
    if l.is_empty() {
        return AffineLattice::empty(a.num_cols());
    }

    // Clean up the solution by taking everything mod n
    // removing the last components that correspond to the multipliers
    // of the n's and then removing (now) linearly dependent basis vectors.

    let n = &n;
    let offset = Vector::from_iter(a.num_cols(), l.offset[0..a.num_cols()]
        .iter()
        .map(|i| i.rem_euclid(n))
    );

    // This might be the worst code in the history of code.
    let iter = l.lattice.basis.rows()
        .flat_map(|e| e.iter().take(a.num_cols()).map(|i| i.rem_euclid(n)))
        .chain((0..a.num_cols()).flat_map(|i| (0..a.num_cols()).map(
            move |j| if j == i { n.clone() } else { Z::zero() }))
        );
    let bm = Matrix::from_iter(l.lattice.rank() + a.num_cols(), a.num_cols(),
        iter
    );

    let lattice = Lattice::from_generating_set(bm, &Z);

    AffineLattice {
        offset,
        lattice,
    }
}

/// Computes a diagonal matrix D in-place and returns matrices (S, T), such
/// that D=SAT.
///
/// The `integer` in the name refers to the fact that this function does not
/// try to compute inverses of elements as opposed to [`modular_diagonalize`].
///
/// So this function is more appropriate for integer matrices.
pub fn integer_diagonalize<R: IntDivRing>(
    a: &mut MatrixView<R>,
    ring: &R,
) -> (OwnedMatrix<R>, OwnedMatrix<R>) {
    // The matrices S and T are initialized to the identity.
    // S/T keeps track of the row/column operations.
    let mut s = Matrix::identity(a.num_rows());
    let mut t = Matrix::identity(a.num_cols());

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
                    .min_by(|l, r| ring.cmp_abs(l.1, r.1))
                    .map(|e| e.0)
                    .unwrap(); // We know there is a non-zero element.

                // Move the pivot to the beginning.
                a.swap_rows(i, pivot);
                s.swap_rows(i, pivot);

                // Try to eliminate every other entry in the column.
                for k in i+1..a.num_rows() {
                    if a[(k, i)].is_zero() {
                        continue;
                    }

                    //let m = ring.neg(R::rounded_div(&a[(k, i)], &a[(i, i)]));
                    let m = ring.neg(R::euclidean_div(&a[(k, i)], &a[(i, i)]));
                    a.row_multiply_add(k, i, &m, ring);
                    s.row_multiply_add(k, i, &m, ring);
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
                .min_by(|l, r| ring.cmp_abs(l.1, r.1))
                .map(|e| e.0)
                .unwrap(); // We know there is a non-zero element.

            // Move the pivot to the beginning.
            a.swap_columns(i, pivot);
            t.swap_columns(i, pivot);

            // Try to eliminate every other entry in the row.
            for k in i+1..a.num_cols() {
                if a[(i, k)].is_zero() {
                    continue;
                }

                //let m = ring.neg(R::rounded_div(&a[(i, k)], &a[(i, i)]));
                let m = ring.neg(R::euclidean_div(&a[(i, k)], &a[(i, i)]));
                a.col_multiply_add(k, i, &m, ring);
                t.col_multiply_add(k, i, &m, ring);
            }
        }
    }

    (s, t)
}

/// Computes a diagonal matrix D in-place and returns matrices (S, T), such
/// that D=SAT.
///
/// The `modular` in the name refers to the fact that this function tries to
/// compute inverses of elements as opposed to [`integer_diagonalize`].
pub fn modular_diagonalize<R: IntDivRing>(
    a: &mut MatrixView<R>,
    ring: &R,
) -> (OwnedMatrix<R>, OwnedMatrix<R>) {
    // The matrices S and T are initialized to the identity.
    // S/T keeps track of the row/column operations.
    let mut s = Matrix::identity(a.num_rows());
    let mut t = Matrix::identity(a.num_cols());

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
                let (pivot, is_one) = 'b: {
                    // Do we have an element that is 1?
                    if let Some((pivot, _)) = a.col(i)
                        .iter()
                        .enumerate()
                        .skip(i)
                        .find(|(_, e)| e.is_one())
                    {
                        break 'b (pivot, true);
                    }

                    // Do we have an element that is a unit?
                    if let Some((pivot, inv)) = a.col(i)
                        .iter()
                        .enumerate()
                        .skip(i)
                        .find_map(|(i, e)| ring.inverse(e).map(|inv| (i, inv)))
                    {
                        // Multiply the row by the inverse.
                        a.row_multiply(pivot, &inv, ring);
                        s.row_multiply(pivot, &inv, ring);
                        break 'b (pivot, true);
                    }

                    // Find the smallest non-zero element by absolute value.
                    let pivot = a.col(i)
                        .iter()
                        .enumerate()
                        .skip(i)
                        .filter(|e| !e.1.is_zero())
                        .min_by(|a, b| ring.cmp_abs(a.1, b.1))
                        .expect("Failed to find pivot in non-zero column").0;

                    (pivot, false)
                };


                // Move the pivot to the beginning.
                a.swap_rows(i, pivot);
                s.swap_rows(i, pivot);

                // Try to eliminate every other entry in the column.
                for k in i+1..a.num_rows() {
                    let e = &a[(k, i)];
                    if e.is_zero() {
                        continue;
                    }

                    let m = if is_one {
                        ring.neg(e.clone())
                    } else {
                        ring.neg(R::euclidean_div(e, &a[(i, i)]))
                    };

                    a.row_multiply_add(k, i, &m, ring);
                    s.row_multiply_add(k, i, &m, ring);
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
            let (pivot, is_one) = 'b: {
                // Do we have an element that is 1?
                if let Some((pivot, _)) = a.row(i)
                    .iter()
                    .enumerate()
                    .skip(i)
                    .find(|(_, e)| e.is_one())
                {
                    break 'b (pivot, true);
                }

                // Do we have an element that is a unit?
                if let Some((pivot, inv)) = a.row(i)
                    .iter()
                    .enumerate()
                    .skip(i)
                    .find_map(|(i, e)| ring.inverse(e).map(|inv| (i, inv)))
                {
                    a.col_multiply(pivot, &inv, ring);
                    t.col_multiply(pivot, &inv, ring);
                    break 'b (pivot, true);
                }

                // Find the smallest non-zero element by absolute value.
                let pivot = a.row(i)
                    .iter()
                    .enumerate()
                    .skip(i)
                    .filter(|e| !e.1.is_zero())
                    .min_by(|l, r| ring.cmp_abs(l.1, r.1))
                    .expect("Failed to find pivot in non-zero row").0;

                (pivot, false)
            };

            // Move the pivot to the beginning.
            a.swap_columns(i, pivot);
            t.swap_columns(i, pivot);

            // Try to eliminate every other entry in the row.
            for k in i+1..a.num_cols() {
                let e = &a[(i, k)];
                if e.is_zero() {
                    continue;
                }

                let m = if is_one {
                    ring.neg(e.clone())
                } else {
                    ring.neg(R::euclidean_div(e, &a[(i, i)]))
                };

                a.col_multiply_add(k, i, &m, ring);
                t.col_multiply_add(k, i, &m, ring);
            }
        }
    }

    (s, t)
}

/// Uses the `diagonalize` function to diagonalize the matrix and then solves
/// the linear system using it.
fn solve_via_diagonalize<R, F>(
    mut a: OwnedMatrix<R>,
    b: OwnedVector<R>,
    ring: &R,
    diagonalize: F,
) -> AffineLattice<R>
where
    R: IntDivRing,
    F: Fn(&mut MatrixView<R>, &R) -> (OwnedMatrix<R>, OwnedMatrix<R>),
{
    assert_eq!(a.num_rows(), b.dim(), "Invalid system of congruences");

    // Diagonalize the system.
    //let mut d = a.clone();
    let (s, t) = diagonalize(a.view_mut(), ring);
    //println!("s: {s:?}\na: {d:?}\nt: {t:?}");
    //assert_eq!(s.mul(&a.mul(&t, ring), ring), d);

    // We could already do this in diagonalize if we really wanted.
    let b = s.mul_vec_post(&b, ring);

    // If there is a non-zero entry in b at index > a.min_dim()
    // then the system has no solution, since the corresponding
    // row in a is zero, so we are solving 0=x.
    if b.iter().skip(a.min_dim()).any(|e| !e.is_zero()) {
        return AffineLattice::empty(a.num_cols());
    }

    // Some solution to the system.
    let mut offset = Vector::zero(a.num_cols());

    // The basis of the kernel.
    let mut basis = OwnedMatrix::zero(0, a.num_cols());

    // Solve the scalar linear congruences.
    for i in 0..a.min_dim() {
        let l = &a[(i, i)];
        let r = &b[i];
        let Some((x, kern)) = solve_scalar_congruence(l, r, ring) else {
            // If there is no solution,
            // then the whole system does not have a solution.
            return AffineLattice::empty(a.num_cols());
        };

        offset[i] = x;

        if !kern.is_zero() {
            let r = basis.num_rows();
            basis.append_zero_rows(1);
            basis[(r, i)] = kern;
        }
    }

    // If there are more variables then equations
    // then there are no restrictions on the variables
    // from index d.rows
    for i in a.num_rows()..a.num_cols() {
        let r = basis.num_rows();
        basis.append_zero_rows(1);
        basis[(r, i)] = R::one();
    }

    offset = t.mul_vec_post(&offset, ring);
    basis = basis.mul(t.transpose(), ring);

    let lattice = Lattice::from_basis(basis);

    AffineLattice {
        offset,
        lattice,
    }
}


/// Solves a system of linear equations Ax=b by "diagonalizing" the matrix
/// with [`integer_diagonalize`]. This should probably not be used as
/// [`solve_via_modular_diagonalize`] is faster.
pub fn solve_via_integer_diagonalize<R: IntDivRing>(
    a: OwnedMatrix<R>,
    b: OwnedVector<R>,
    ring: &R,
) -> AffineLattice<R> {
    solve_via_diagonalize(a, b, ring, integer_diagonalize)
}

/// Solves a system of linear equations Ax=b by "diagonalizing" the matrix
/// with [`modular_diagonalize`].
pub fn solve_via_modular_diagonalize<R: IntDivRing>(
    a: OwnedMatrix<R>,
    b: OwnedVector<R>,
    ring: &R,
) -> AffineLattice<R> {
    solve_via_diagonalize(a, b, ring, modular_diagonalize)
}

/// Solves `ax=b` in the ring. Returns None if there is no solution. Otherwise
/// returns all solutions in the form `(c, d)` where `c+di` are all solutions
/// where `i` ranges over all elements of the ring.
///
/// This is only meant to work for rings that are the integers mod n.
pub fn solve_scalar_congruence<R: IntDivRing>(
    a: &R::Element,
    b: &R::Element,
    ring: &R,
) -> Option<(R::Element, R::Element)> {
    // Handle the case that a is zero, so we don't have to think about it.
    if a.is_zero() {
        return b.is_zero().then_some((R::zero(), R::one()));
    }

    // We are basically going to use the extended euclidean algorithm on
    // the diophantine equation ax+ny=b where n is the number of values
    // (2^8 for u8).
    // But n doesn't fit into T, so we have to do a hack in the first step.
    // Usually we'd divide n by a but instead we divide n-a by a and add 1.
    // This makes the code structurally uglier, but otherwise I'm pretty
    // much just following the pseudo code on wikipedia.
    let (mut old_r, mut r) = (R::zero(), a.clone());
    let (mut old_t, mut t) = (R::zero(), R::one());
    let mut q = ring.inc(R::euclidean_div(&ring.neg(a.clone()), a));

    loop {
        let new_r = ring.sub_rhs(&old_r, ring.mul(q.clone(), &r));
        let new_t = ring.sub_rhs(&old_t, ring.mul(q.clone(), &t));
        (old_r, r) = (r, new_r);
        (old_t, t) = (t, new_t);
        if r.is_zero() {
            break;
        }

        q = R::euclidean_div(&old_r, &r);
    }

    // old_r is gcd(a, n).
    let gcd = old_r;

    // There is a solution iff gcd divides b.
    // old_t is the Bezout coefficient: a*old_t=gcd(a, n) mod n.
    let x = ring.mul(R::euclidean_div(b, &gcd), &old_t);
    if &ring.mul(a.clone(), &x) != b {
        return None;
    }

    // We compute the absolute value of `t` because its equivalence class is
    // represented by the smaller integer. This isn't strictly necessary.
    ring.abs_assign(&mut t);

    Some((x, t))
}

#[cfg(test)]
mod test {
    use rand::{distr::{Distribution as _, Uniform}, rngs::StdRng, SeedableRng};

    use crate::rings::{BinaryRing, VarU8, U128, U16, U32, U64, U8};

    use super::*;

    /// Tests that all methods of solving a system of linear equations give the
    /// same result.
    fn same_result<R: BinaryRing>(
        a: &MatrixView<R>,
        b: &VectorView<R>,
        ring: &R,
    ) {
        let a_int = a.transform(|e| R::to_representative(e).into());
        let b_int = b.transform(|e| R::to_representative(e).into());
        let l1 = solve_modular_via_integer_hnf(a_int.view(), b_int.view(), &(BigUint::one() << ring.bits()));
        let l1_c = {
            let mut l = l1.clone();
            l.offset.reduce(&l1.lattice.basis, &Z);
            l
        };
        let l2 = solve_via_integer_diagonalize(a.to_owned(), b.to_owned(), ring);
        let l2_c = l2.canonicalize(ring);
        let l3 = solve_via_modular_diagonalize(a.to_owned(), b.to_owned(), ring);
        let l3_c = l3.canonicalize(ring);

        assert_eq!(l1_c, l2_c, "a: {a:?}, b: {b:?}, l1: {l1:?}, l2: {l2:?}");
        assert_eq!(l2_c, l3_c, "a: {a:?}, b: {b:?}, l2: {l2:?}, l3: {l3:?}");
    }

    /// Generates a random system of linear equations and tests that all methods
    /// of solving it give the same result.
    fn random_test<R: BinaryRing>(
        ring: &R,
    ) {
        let rng = &mut StdRng::seed_from_u64(0);
        let dist = Uniform::new(0, 20).unwrap();
        for _i in 0..1000 {
            let dim = dist.sample(rng);
            let cols = dist.sample(rng);
            let a = Matrix::random(dim, cols, ring, rng);
            let b = Vector::random(dim, ring, rng);

            same_result(a.view(), b.view(), ring);
        }
    }

    #[test]
    fn random_test_u2() {
        random_test(&VarU8::new(2));
    }

    #[test]
    fn random_test_u8() {
        random_test(&U8);
    }

    #[test]
    fn random_test_u16() {
        random_test(&U16);
    }

    #[test]
    fn random_test_u32() {
        random_test(&U32);
    }

    #[test]
    fn random_test_u64() {
        random_test(&U64);
    }

    #[test]
    fn random_test_u128() {
        random_test(&U128);
    }

    #[test]
    fn solve_modular_via_integer_hnf_test() {
        let a = OwnedMatrix::<Z>::from_rows(&[
            [0, 3],
            [2, 1],
        ]);
        let b = OwnedVector::<Z>::from_entries([0, 3]);
        let l = solve_modular_via_integer_hnf(a.view(), b.view(), &4u32.into());
        assert!(l.is_empty());
    }

    #[test]
    fn scalar_congruence_test() {
        use crate::rings::BigIntModN;
        let ring = BigIntModN::new(10u32.into());
        let a = BigUint::from(2u32);
        let b = BigUint::from(8u32);
        let (x, kern) = solve_scalar_congruence(
            &a,
            &b,
            &ring,
        ).unwrap();
        println!("{x} + {kern}*i");
        assert_eq!((&a * &x).rem_euclid(ring.modulus()),
            b.rem_euclid(ring.modulus()));
        assert!((&a * &kern).rem_euclid(ring.modulus()).is_zero());
    }

    #[test]
    fn small_test() {
        let a = OwnedMatrix::<Z>::from_rows(&[
            [0, 0, -1, 1],
            [0, 1, 0, 0],
            [0, 1, 0, 0],
            [1, 0, 0, 0],
        ]);

        let b = OwnedVector::<Z>::from_entries([0, 1, 1, 2]);

        let l = solve_via_integer_hnf(a.view(), b.view(), &Z);
        assert_eq!(l.offset.as_slice(),
            [2.into(), 1.into(), 0.into(), 0.into()]);
        assert_eq!(l.lattice.rank(), 1);
        assert_eq!(l.lattice.basis.row(0).as_slice(),
            [0.into(), 0.into(), 1.into(), 1.into()]);
    }
}