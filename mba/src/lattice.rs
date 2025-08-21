//! Integer lattices and algorithms.

use crate::rings::{
    BinaryRing, F32, F64, Field, IntDivRing, OrderedRing, Q, Ring, RingElement as _, SqrtRing, Z,
};
use crate::solver::hermite_normal_form;
use crate::{matrix::*, vector::*};
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{FromPrimitive, ToPrimitive};
use rand::Rng;
use std::fmt::Debug;

/// A "lattice", i.e. the span of a set of vectors.
/// This can never be empty because it always contains the zero vector even if
/// [`Lattice::basis`] has zero rows. The matrix always stores the ambient
/// dimension (the number of columns) even if there are zero rows.
#[derive(Clone, Debug, PartialEq)]
pub struct Lattice<R: Ring> {
    /// The basis matrix of the lattice.
    /// The basis vectors are the rows of the matrix.
    pub basis: OwnedMatrix<R>,
}

impl<R: Ring> Lattice<R> {
    /// Creates a lattice that contains only the zero vector.
    pub fn zero(ambient_dim: usize) -> Self {
        Self {
            basis: OwnedMatrix::zero(0, ambient_dim),
        }
    }

    /// The lattice basis are the rows of the matrix.
    pub fn from_basis(basis: OwnedMatrix<R>) -> Self {
        Self { basis }
    }

    /// The rows of the matrix generate the lattice
    /// but are potentially linearly dependent.
    /// This function will compute the Hermite normal form
    /// and remove zero rows.
    pub fn from_generating_set(mut generating_set: OwnedMatrix<R>, r: &R) -> Self
    where
        R: IntDivRing,
    {
        hermite_normal_form(&mut generating_set, r);
        generating_set.remove_zero_rows();
        Self {
            basis: generating_set,
        }
    }

    /// Returns the rank of the lattice, i.e. the number if basis vectors.
    pub fn rank(&self) -> usize {
        self.basis.num_rows()
    }

    /// Returns the dimension of the ambient space or 0 if the lattice is empty.
    pub fn ambient_dim(&self) -> usize {
        self.basis.num_cols()
    }

    /// Returns the vector on the lattice that is the linear combination of the
    /// basis vectors with the given coefficients.
    pub fn at<S: VectorStorage<R> + ?Sized>(
        &self,
        coefficients: &Vector<R, S>,
        r: &R,
    ) -> OwnedVector<R> {
        self.at_impl(coefficients, OwnedVector::zero(self.ambient_dim()), r)
    }

    /// Adds the vector on the lattice that is the linear combination of the
    /// basis vectors with the given `coefficients` to the `initial` vector.
    pub(self) fn at_impl<S: VectorStorage<R> + ?Sized>(
        &self,
        coefficients: &Vector<R, S>,
        initial: OwnedVector<R>,
        r: &R,
    ) -> OwnedVector<R> {
        assert!(coefficients.dim() == self.rank());
        self.basis
            .rows()
            .zip(coefficients)
            .fold(initial, |acc, (e, c)| acc.mul_add(c, e, r))
    }

    /// Returns a random point on the lattice.
    pub fn sample_point<Rand: Rng>(&self, rng: &mut Rand, r: &R) -> OwnedVector<R> {
        self.sample_point_impl(Vector::zero(self.ambient_dim()), rng, r)
    }

    /// Samples a point from the lattice that is added to the initial vector.
    pub(self) fn sample_point_impl<Rand: Rng>(
        &self,
        initial: OwnedVector<R>,
        rng: &mut Rand,
        r: &R,
    ) -> OwnedVector<R> {
        assert!(initial.dim() == self.ambient_dim());

        self.basis
            .rows()
            .fold(initial, |acc, b| acc.mul_add(&r.random(rng), b, r))
    }

    /// Returns an integer lattice such that two lattices are equal (i.e.
    /// contain the same points) iff `l1.canonicalize() == l2.canonicalize()`.
    ///
    /// The basis of the returned lattice is guaranteed to be in Hermite normal
    /// form and of full rank. In particular, the diagonal entries are positive.
    pub fn canonicalize(&self, ring: &R) -> Lattice<Z>
    where
        R: BinaryRing,
    {
        // Transform the matrix and vector to integer matrices and vectors.
        let mut a = self.basis.transform(|e| R::to_representative(e).into());

        // Append rows that contain the modulus.
        a.append_zero_rows(a.num_cols());
        for i in 0..a.num_cols() {
            a[(self.basis.num_rows() + i, i)] = Z::one() << ring.bits();
        }

        // Compute the HNF of the integer matrix.
        let _ = hermite_normal_form(&mut a, &Z);

        a.remove_zero_rows();

        Lattice { basis: a }
    }

    /// Size reduce the basis.
    /// This essentially is Gram-Schmidt
    /// but rounding the coefficients to integers.
    pub fn size_reduce(&mut self, r: &R)
    where
        R: IntDivRing,
    {
        size_reduce(&mut self.basis, r);
    }

    /// Performs LLL basis reduction.
    pub fn lll<W: WorkingTypeFor<R>>(&mut self, delta: &W::Element, wt: &W, r: &R)
    where
        R: IntDivRing,
    {
        lll(&mut self.basis, delta, wt, r);
    }

    pub fn cvp_planes_coeff<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        rad_sqr: Option<W::Element>,
        wt: &W,
        r: &R,
    ) -> Option<OwnedVector<R>> {
        cvp_planes(&self.basis, t, rad_sqr, wt, r)
    }

    pub fn cvp_planes<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        rad_sqr: Option<W::Element>,
        wt: &W,
        r: &R,
    ) -> Option<OwnedVector<R>> {
        Some(self.at(&self.cvp_planes_coeff(t, rad_sqr, wt, r)?, r))
    }

    pub fn cvp_rounding_coeff<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        wt: &W,
        r: &R,
    ) -> OwnedVector<R> {
        cvp_rounding(&self.basis, t, wt, r)
    }

    pub fn cvp_rounding<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        wt: &W,
        r: &R,
    ) -> OwnedVector<R> {
        self.at(&cvp_rounding(&self.basis, t, wt, r), r)
    }

    pub fn cvp_nearest_plane_coeff<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        wt: &W,
        r: &R,
    ) -> OwnedVector<R> {
        cvp_nearest_plane(&self.basis, t, true, wt, r).1
    }

    pub fn cvp_nearest_plane<W: WorkingTypeFor<R>>(
        &self,
        t: &VectorView<R>,
        wt: &W,
        r: &R,
    ) -> OwnedVector<R> {
        cvp_nearest_plane(&self.basis, t, false, wt, r).0
    }
}

/// A lattice that is offset from the origin by a vector.
/// Mathematicians would call this a lattice coset.
///
/// In contrast to [`Lattice`], this can be empty which is represented by an
/// offset vector of dimension 0. The [`AffineLattice::lattice`] still stores
/// the ambient dimension.
#[derive(Clone, Debug, PartialEq)]
pub struct AffineLattice<R: Ring> {
    pub offset: OwnedVector<R>,
    pub lattice: Lattice<R>,
}

impl<R: Ring> AffineLattice<R> {
    /// Creates an empty affine lattice that contains no points.
    pub fn empty(ambient_dim: usize) -> Self {
        Self {
            offset: Vector::empty(),
            lattice: Lattice::zero(ambient_dim),
        }
    }

    /// Creates an affine lattice from an offset and a basis.
    pub fn from_offset_basis(offset: OwnedVector<R>, basis: OwnedMatrix<R>) -> Self {
        Self {
            offset,
            lattice: Lattice::from_basis(basis),
        }
    }

    /// The ambient dimension of the lattice.
    pub fn ambient_dim(&self) -> usize {
        self.lattice.ambient_dim()
    }

    /// Is this lattice empty?
    pub fn is_empty(&self) -> bool {
        self.offset.is_empty()
    }

    /// Returns the vector on the lattice that is the linear combination of the
    /// basis vectors with the given coefficients.
    pub fn at<S: VectorStorage<R> + ?Sized>(
        &self,
        coefficients: &Vector<R, S>,
        r: &R,
    ) -> OwnedVector<R> {
        self.lattice.at_impl(coefficients, self.offset.clone(), r)
    }

    /// Returns a random point on the lattice.
    pub fn sample_point<Rand: Rng>(&self, rng: &mut Rand, r: &R) -> OwnedVector<R> {
        self.lattice.sample_point_impl(self.offset.clone(), rng, r)
    }

    /// Creates a matrix `M` such that affine lattice `L1` and `L2` are equal,
    /// i.e. contain the same points, iff `canonicalize(L1) = canonicalize(L2)`.
    pub fn canonicalize(&self, ring: &R) -> AffineLattice<Z>
    where
        R: BinaryRing,
    {
        if self.is_empty() {
            return AffineLattice::empty(self.ambient_dim());
        }

        let lattice = self.lattice.canonicalize(ring);
        let mut b = self.offset.transform(|e| R::to_representative(e).into());

        // Reduce the entries of the vector.
        b.reduce(&lattice.basis, &Z);

        AffineLattice { offset: b, lattice }
    }

    /// Returns the shortest vector on the **affine(!)** lattice.
    /// This solves the CVP on the normal lattice with
    /// target offset.
    pub fn svp<W: WorkingTypeFor<R>>(
        &self,
        rad: Option<W::Element>,
        wt: &W,
        r: &R,
    ) -> Option<OwnedVector<R>> {
        self.lattice
            .cvp_planes(self.offset.view(), rad, wt, r)
            .map(|v| self.offset.sub_rhs(v, r))
    }
}

/// Size reduce the basis.
/// This essentially is Gram-Schmidt
/// but rounding the coefficients to integers.
pub fn size_reduce<R: IntDivRing>(basis: &mut OwnedMatrix<R>, r: &R) {
    for i in 0..basis.num_rows() {
        for j in 0..i {
            let (b_i, b_j) = basis.get_rows_mut(i, j);
            let dot = b_i.dot(b_j, r);
            let b_j_norm_sqr = b_j.norm_sqr(r);
            let q = R::rounded_div(&dot, &b_j_norm_sqr);

            if q.is_zero() {
                continue;
            }

            let q = r.neg(q);
            b_i.mul_add_assign(&q, b_j, r);
        }
    }
}

/// Performs LLL basis reduction using rational numbers.
pub fn lll<R: IntDivRing, W: WorkingTypeFor<R>>(
    basis: &mut OwnedMatrix<R>,
    delta: &W::Element,
    wt: &W,
    r: &R,
) {
    let n = basis.num_rows();
    let mut did_swap = true;

    while did_swap {
        size_reduce(basis, r);

        // Lovasz condition
        did_swap = false;
        for i in 0..n - 1 {
            let b = &basis[i];
            let c = &basis[i + 1];

            let lhs = wt.from_ring(&c.norm_sqr(r), r);

            let b_norm_sqr = wt.from_ring(&b.norm_sqr(r), r);
            let q = wt.div(wt.from_ring(&c.dot(b, r), r), &b_norm_sqr);
            let rhs = wt.mul(b_norm_sqr, &(wt.add(wt.neg(wt.square(q)), delta)));

            if wt.is_lt(&lhs, &rhs) {
                basis.swap_rows(i, i + 1);
                did_swap = true;
                break;
            }
        }
    }
}

/// Approximates the CVP using Babai's rounding technique.
/// The returns value is a vector of the coefficients
/// of the linear combination of the basis vectors,
/// so if you want the actual point,
/// you still have to matrix multiply with the basis matrix.
/// In practice, it is a good idea to reduce the basis (e.g. using LLL)
/// so that the approximation is good.
pub fn cvp_rounding<R: Ring, W: WorkingTypeFor<R>>(
    basis: &OwnedMatrix<R>,
    t: &VectorView<R>,
    wt: &W,
    ring: &R,
) -> OwnedVector<R> {
    let mut a = OwnedMatrix::zero(basis.num_cols(), basis.num_rows());

    for r in 0..a.num_rows() {
        for c in 0..a.num_cols() {
            a[(r, c)] = wt.from_ring(&basis[(c, r)], ring);
        }
    }

    let b = OwnedVector::from_iter(t.dim(), t.iter().map(|i| wt.from_ring(i, ring)));

    // If the system has full rank,
    // then we can just use an exact solution.
    let x = if t.dim() == basis.num_rows() {
        solve_linear(a, b, wt)
    }
    // Otherwise we can treat it as a Ordinary Least Squares problem,
    // i.e. find the point in the subspace spanned by the vectors
    // that is closest to the given point.
    else {
        let a_t = a.transpose();
        solve_linear(a_t.mul(&a, wt), a_t.mul_vec_post(&b, wt), wt)
    };

    let x = x.expect("Basis vectors were not independent.");

    let mut res = OwnedVector::zero(x.dim());
    for (i, f) in res.iter_mut().zip(x.iter()) {
        *i = W::to_ring(f, ring);
    }

    res
}

/// Gram-Schmidt of the rows of the matrix `a`.
/// The rows will span the same subspace as the original rows,
/// even if they are not linearly independent.
/// The rows of the result are not normalized.
/// Call `gram_schmidt_orthonormal` if you want that.
pub fn gram_schmidt<R: Field>(mut a: OwnedMatrix<R>, r: &R) -> OwnedMatrix<R> {
    for i in 0..a.num_rows() {
        for j in 0..i {
            let (a_i, a_j) = a.get_rows_mut(i, j);
            let f = r.neg(r.div(a_j.dot(a_i, r), &a_j.norm_sqr(r)));
            a_i.mul_add_assign(&f, a_j, r);
        }
    }

    a
}

/// Gram-Schmidt of the rows of the matrix `a`.
/// The rows of the result are normalized.
/// Call `gram_schmidt` if you don't want that.
/// This functions requires that the rows are linearly independent
/// (unlike `gram_schmidt`).
pub fn gram_schmidt_orthonormal<R: Field + SqrtRing>(
    mut a: OwnedMatrix<R>,
    r: &R,
) -> OwnedMatrix<R> {
    a = gram_schmidt(a, r);
    for i in 0..a.num_rows() {
        a[i].normalize(r);
    }
    a
}

/// Computes the RQ-decomposition (row QR-decomposition) of a matrix.
/// Returns a lower triangular matrix `R` and an orthogonal matrix `Q`,
/// such that `R * Q = A`.
pub fn rq_decomposition<R: Field + SqrtRing>(
    a: &OwnedMatrix<R>,
    r: &R,
) -> (OwnedMatrix<R>, OwnedMatrix<R>) {
    let q = gram_schmidt_orthonormal(a.clone(), r);
    let r = a.mul(q.transpose(), r);
    (r, q)
}

/// Approximates the CVP using Babai's nearest plane algorithm.
/// The return value is a pair of the vector and the vector of coefficients
/// of that vector if `compute_coefficients` is true.
/// Otherwise the second element of the pair is an empty vector.
pub fn cvp_nearest_plane<R: Ring, W: WorkingTypeFor<R>>(
    basis: &OwnedMatrix<R>,
    t: &VectorView<R>,
    compute_coefficients: bool,
    wt: &W,
    r: &R,
) -> (OwnedVector<R>, OwnedVector<R>) {
    let q = gram_schmidt(basis.transform(|i| wt.from_ring(i, r)), wt);

    let mut off: OwnedVector<R> = t.to_owned();
    let b = off.transform::<W, _>(|i| wt.from_ring(i, r));

    let mut coeff = if compute_coefficients {
        OwnedVector::zero(basis.num_rows())
    } else {
        OwnedVector::empty()
    };

    for i in (0..basis.num_rows()).rev() {
        let a = &q[i];
        // c = round(⟨a, b⟩ / ⟨a, a⟩)
        let c = b
            .iter()
            .zip(a.iter())
            .fold(W::Element::zero(), |acc, (i, f)| wt.mul_add(acc, i, f));
        let c = wt.div(c, &a.norm_sqr(wt));
        let c = r.neg(W::to_ring(&c, r));
        off.mul_add_assign(&c, &basis[i], r);
        if compute_coefficients {
            coeff[i] = c;
        }
    }

    (t.sub_rhs(off, r), coeff)
}

/// Solves the CVP exactly.
/// - `basis` is the matrix that contains the basis as rows.
/// - `t` is the target vector.
/// - `prec` is the precision of the floating point numbers used.
/// - `rad_sqr` is an optional float that contains the square of the maximum
///   distance to search for.
///
/// The returned vector is the vector of coefficients.
/// This will always return some vector unless no vector is within `r` of the target.
/// This algorithm is a generalization Babai's nearest plane algorithm
/// that searches all planes that could contain the closest vector.
/// It is the simplest one I could think of.
pub fn cvp_planes<R: Ring, W: WorkingTypeFor<R>>(
    basis: &OwnedMatrix<R>,
    t: &VectorView<R>,
    rad_sqr: Option<W::Element>,
    wt: &W,
    ring: &R,
) -> Option<OwnedVector<R>> {
    assert!(
        basis.num_cols() == t.dim(),
        "Mismatch of basis/target vector dimension."
    );
    let bf = basis.transform::<W, _>(|i| wt.from_ring(i, ring));

    // Q is the Gram-Schmidt orthonormalization and
    // R is the basis matrix with respect to the Gram-Schmidt basis.
    let (r, q) = rq_decomposition(&bf, wt);
    assert!(r.num_cols() == r.num_rows());

    // Write the target vector in that basis.
    // If the vector is not in the span of the basis,
    // then this will project it into the span.
    let qt = q.mul_vec_post(&t.transform(|i| wt.from_ring(i, ring)), wt);

    // We just do this to avoid having to pass around an Option.
    let rad = rad_sqr.unwrap_or_else(|| wt.infinity());

    // Multiply the coefficients of the vector.
    return cvp_impl(r.num_cols() - 1, &r, qt.view(), &rad, wt, ring).map(|v| v.0);

    /// This actually finds the closest point.
    /// `rad` is the squared norm.
    /// This returns the coordinates of the closest point in the basis
    /// as the first entry and the actual point as the second entry.
    fn cvp_impl<R: Ring, W: WorkingTypeFor<R>>(
        i: usize,
        r: &OwnedMatrix<W>,
        qt: &VectorView<W>,
        rad: &W::Element,
        wt: &W,
        ring: &R,
    ) -> Option<(OwnedVector<R>, OwnedVector<W>)> {
        // One dimensional lattice.
        if i == 0 {
            let qt = &qt[0];
            let r = &r[(0, 0)];
            // `m` is the index of the closest point.
            // `d` is the distance to it.
            let m = W::round(wt.div(qt.clone(), r));
            let plane = wt.mul(r.clone(), &m);
            let d = wt.square(wt.sub(plane.clone(), qt));
            return wt.is_le(&d, rad).then(|| {
                (
                    OwnedVector::from_array([W::to_ring(&m, ring)]),
                    OwnedVector::from_array([plane]),
                )
            });
        }

        let qtc = &qt[i];
        let rc = &r[(i, i)];

        // Index of the closest plane before rounding.
        let start_index_fl = wt.div(qtc.clone(), rc);
        let start_index = W::round(start_index_fl.clone());

        // Suppose the start index was -0.4.
        // We would want to check the planes in the order
        // 0, -1, 1, -2, 2.
        // But if the start index was 0.4, we would want
        // 0, 1, -1, 2, -2.
        // So depending on whether we round up or down
        // to the integer start index, we will negate
        // the offset from the start index.
        let negate_offset = wt.is_lt(&start_index_fl, &start_index);

        // The current plane's offset.
        let mut offset = 0isize;

        let mut min_dist = rad.clone();
        let mut min = None;

        // Iterate over all possible planes.
        loop {
            // Index of the plane we are considering.
            let index = W::add_isize(start_index.clone(), offset);

            // Compute the index offset of the next plane.
            // Negate the offset first.
            offset = -offset;

            // If we negate the offsets, i.e. 0, -1, 1, -2, ...,
            // then if we are <= 0, we need to subtract 1.
            // E.g. if the offset was 1, then we negated it to -1
            // and subtract 1.
            if negate_offset && offset <= 0 {
                offset -= 1;
            }
            // In the 0, 1, -1, 2, ... case,
            // if the offset is >= 0, we need to add 1.
            // E.g. if the offset was -1, then we negated it to 1
            // and add 1.
            else if !negate_offset && offset >= 0 {
                offset += 1;
            }

            // If this overflows, you probably would have
            // ctrl-c'd before we got here.
            // debug_assert!(offset != 0);

            // Distance to the plane.
            // This is the distance of the target to
            // its orthogonal projection in the plane.
            let d = wt.square(wt.sub(wt.mul(index.clone(), rc), qtc));

            // If the plane is not in the radius,
            // then the next one in the loop definitely is not
            // either, because of the way we iterate over the planes.
            if wt.is_gt(&d, &min_dist) {
                break;
            }

            // We can use a smaller radius inside the plane.
            // The target to its projection to any point
            // in the plane form a right triangle.
            // So by Pythagoras the squared distance in the
            // plane can only be <= rad - d.
            // rad and d are already the square of the distance.
            let plane_dist = wt.sub(min_dist.clone(), &d);

            // Compute the point in the plane we need to be close to now.
            let point_in_plane = qt.sub_rhs(r[i].mul_ref(&index, wt), wt);

            // Recursively find the closest point.
            let Some((mut v, mut w)) =
                cvp_impl(i - 1, r, point_in_plane.view(), &plane_dist, wt, ring)
            else {
                continue;
            };
            assert!(v.dim() == i);

            // w is the closest point.
            // It is the closest point of the previous call,
            // plus the index of the plane times the basis vector
            // of the current iteration.
            w.mul_add_assign(&index, &r[i][0..i], wt);
            w.append(wt.mul(index.clone(), rc));

            // Compute the distance to the point.
            let d = w.dist_sqr(&qt[0..i + 1], wt);

            // If the distance is smaller than the current minimal dist,
            // then we have found the new best point.
            if wt.is_le(&d, &min_dist) {
                // v is the new coordinate vector of the point.
                // It is the coordinates of the closest point of
                // the previous call, concat the index of the plane.
                v.append(W::to_ring(&index, ring));

                min = Some((v, w));
                min_dist = d;
            }
        }

        min
    }
}

/// Solves a square system of linear equations.
fn solve_linear<W: WorkingType>(
    mut a: OwnedMatrix<W>,
    mut b: OwnedVector<W>,
    wt: &W,
) -> Option<OwnedVector<W>> {
    assert!(
        a.num_rows() == a.num_cols(),
        "This function only supports non-singular square systems."
    );
    let rank = a.num_rows();
    for i in 0..rank {
        // Choose a pivot in the c-th column.
        let pivot = a
            .col(i)
            .iter()
            .enumerate()
            .skip(i)
            .filter(|e| !e.1.is_zero())
            .max_by(|e, f| wt.cmp_abs(e.1, f.1))?
            .0;

        // Swap the pivot row with the current row.
        a.swap_rows(pivot, i);
        for r in i + 1..rank {
            let (pivot_row, row) = a.get_rows_mut(i, r);
            let fac = wt.neg(wt.div(row[i].clone(), &pivot_row[i]));
            for c in i + 1..rank {
                wt.mul_add_assign(&mut row[c], &fac, &pivot_row[c]);
            }
            let (b_r, b_i) = b.get_mut_entries(r, i);
            wt.mul_add_assign(b_r, &fac, b_i);
        }
    }

    let mut result = Vector::zero(rank);
    for i in (0..rank).rev() {
        let row = a.row(i);
        let sum = (i + 1..rank).fold(b[i].clone(), |acc, j| wt.mul_sub(acc, &row[j], &result[j]));

        result[i] = wt.div(sum, &row[i]);
    }

    Some(result)
}

/// Lattice algorithms are possibly internally working with a different ring.
/// These rings need to implement this.
pub trait WorkingType: Field + OrderedRing + SqrtRing {
    fn eps(&self) -> Self::Element;

    fn infinity(&self) -> Self::Element {
        panic!("`infinity` is not supported for this type.");
    }

    fn round(s: Self::Element) -> Self::Element;

    fn add_isize(s: Self::Element, i: isize) -> Self::Element;
}

/// This trait facilitates the conversion between the ring of the lattice and
/// the ring used internally.
pub trait WorkingTypeFor<R: Ring>: WorkingType {
    #[allow(clippy::wrong_self_convention)]
    fn from_ring(&self, e: &R::Element, r: &R) -> Self::Element;

    fn to_ring(s: &Self::Element, r: &R) -> R::Element;
}

macro_rules! float_working_type {
    ($ring:ident) => {
        impl WorkingType for $ring {
            fn eps(&self) -> Self::Element {
                <$ring as Ring>::Element::EPSILON
            }

            fn infinity(&self) -> Self::Element {
                <$ring as Ring>::Element::INFINITY
            }

            fn round(s: Self::Element) -> Self::Element {
                s.round()
            }

            fn add_isize(s: Self::Element, i: isize) -> Self::Element {
                s + i as Self::Element
            }
        }
    };
}

float_working_type!(F32);
float_working_type!(F64);

impl WorkingTypeFor<Z> for F32 {
    fn from_ring(&self, e: &BigInt, _: &Z) -> Self::Element {
        e.to_f32().unwrap()
    }

    fn to_ring(s: &Self::Element, _: &Z) -> <Z as Ring>::Element {
        BigInt::from_f32(s.round()).unwrap()
    }
}

impl WorkingTypeFor<Z> for F64 {
    fn from_ring(&self, e: &BigInt, _: &Z) -> Self::Element {
        e.to_f64().unwrap()
    }

    fn to_ring(s: &Self::Element, _: &Z) -> <Z as Ring>::Element {
        BigInt::from_f64(s.round()).unwrap()
    }
}

impl WorkingType for Q {
    fn eps(&self) -> Self::Element {
        Q::zero()
    }

    fn round(s: Self::Element) -> Self::Element {
        s.round()
    }

    fn add_isize(s: Self::Element, i: isize) -> Self::Element {
        s + BigRational::from_isize(i).unwrap()
    }
}

impl WorkingTypeFor<Z> for Q {
    fn from_ring(&self, e: &BigInt, _: &Z) -> Self::Element {
        BigRational::from_integer(e.clone())
    }

    fn to_ring(s: &Self::Element, _: &Z) -> <Z as Ring>::Element {
        s.round().numer().clone()
    }
}

#[test]
fn cvp_exact_dim20() {
    let l = Lattice::<Z>::from_basis(OwnedMatrix::from_rows(&[
        [
            1, -74, 20, 19, -5, 21, -19, 18, -54, -56, -40, -38, -58, 54, -34, -1, -3, 0, -27, 8,
        ],
        [
            -44, 19, -20, 14, -34, -61, 53, -31, 42, 42, 27, 40, -77, -59, 2, -26, 9, 3, -49, 33,
        ],
        [
            -19, -12, 18, -31, 63, 8, 52, 52, -50, -19, 22, -1, 51, -8, -57, 70, -62, -45, 48, -51,
        ],
        [
            -66, -5, 15, 34, 2, -50, 82, 28, 16, 18, 17, 66, 23, -38, 39, -36, -66, 19, -64, -62,
        ],
        [
            -15, 10, -65, -46, -55, -22, -88, -49, 46, -19, -4, -6, -5, -20, 23, -24, -86, -29, -7,
            16,
        ],
        [
            35, -57, 9, 41, 9, -7, 63, 11, -72, -23, -2, -103, 62, -9, 64, 1, 25, 10, 48, -41,
        ],
        [
            18, 11, 31, 2, -51, 48, 11, -33, 69, 93, -12, -52, -3, -100, 14, -15, 85, -61, 6, 27,
        ],
        [
            -58, -54, -29, -28, 93, -88, 3, -63, -3, -9, 26, 22, 89, 26, 11, -2, -21, -5, 37, -36,
        ],
        [
            24, 72, -79, -38, 50, 17, -54, -24, 38, -9, -64, -32, -10, 70, -67, -88, -4, -34, -4,
            -3,
        ],
        [
            -24, -54, 11, 34, -14, -5, -132, 46, 51, 67, 24, -3, -10, -6, 22, 38, -15, -23, 29, -39,
        ],
        [
            72, 9, 49, -19, 66, -6, -53, -77, 40, -9, -52, -20, 90, -12, 58, -107, -47, -64, -66,
            -10,
        ],
        [
            48, -28, 41, -76, 17, -13, -20, -16, -15, 75, -30, 51, -55, 31, -50, 3, -60, -34, 13,
            31,
        ],
        [
            62, -9, 13, -10, 39, 50, 81, 94, -38, 7, -62, -49, 34, 61, 45, 30, -51, -78, -70, -4,
        ],
        [
            34, 92, -16, -3, 113, -24, 1, 40, -30, 91, -57, -6, -28, -7, 13, 64, 18, 24, -33, -10,
        ],
        [
            -1, -20, -45, 44, -75, 31, -9, -47, 74, -7, 64, 77, -41, 9, 52, -33, -83, 118, 63, 1,
        ],
        [
            -9, -37, 11, -106, -13, 1, 74, 11, 89, -10, 61, 43, -17, -45, -7, 5, -103, 43, -36, 46,
        ],
        [
            60, -101, -48, -23, 37, -45, 1, 23, 52, -18, -19, -36, -125, -23, 22, -32, -28, -41,
            16, -34,
        ],
        [
            19, -47, -85, 17, 2, 1, 12, 19, 27, -21, 43, 43, -9, 8, -60, 36, 83, 65, -50, 58,
        ],
        [
            59, -4, 51, 36, -10, -12, -19, -54, 93, 12, 23, 31, 77, 18, 45, 57, 46, 52, -21, -51,
        ],
        [
            -57, 49, 26, 10, -22, 32, -52, -71, -9, 31, 45, 28, -28, -30, 24, 44, 88, 2, -63, -105,
        ],
    ]));
    let t = OwnedVector::<Z>::from_entries([
        -48, 69, -76, 36, -72, 31, -53, -7, 54, 74, 6, -82, -13, -32, 7, 53, -60, -44, 38, -97,
    ]);
    assert_eq!(
        l.cvp_planes(t.view(), None, &F64, &Z).unwrap(),
        [
            -30, 35, -98, 61, -27, 75, -32, -3, 70, 8, 3, -77, -29, -103, 61, 58, -71, 41, 37, -40
        ]
        .iter()
        .map(|i| BigInt::from(*i))
        .collect::<Vec<_>>()
    );
}

#[test]
fn gram_schmidt_test() {
    let a = OwnedMatrix::<F64>::from_rows(&[[1., 2., 3.], [3., 4., 5.]]);
    let (r, q) = rq_decomposition(&a, &F64);
    println!("q: {q:?}\nr: {r:?}");
    println!("{:?}", r.mul(&q, &F64));
    println!(
        "{:?}",
        q.mul_vec_post(VectorView::from_slice(&[5., 6., 7.]), &F64)
    );
}

#[test]
fn nearest_plane_example() {
    let l = Lattice::from_basis(OwnedMatrix::<Z>::from_rows(&[
        [2, 3, 1],
        [4, 1, -3],
        [2, 2, 2],
    ]));

    let t = OwnedVector::from_entries([4, 2, 7]);
    let c = l.cvp_nearest_plane(t.view(), &F64, &Z);
    println!("{c:?}");
}

#[test]
fn babai_rounding_example() {
    use crate::rings::BinaryBigInt;
    let generators = OwnedMatrix::<Z>::from_rows(&[
        [-97, 75, -97, 75, 22],
        [101, 38, 101, 38, 117],
        [256, 0, 0, 0, 0],
        [0, 256, 0, 0, 0],
        [0, 0, 256, 0, 0],
        [0, 0, 0, 256, 0],
        [0, 0, 0, 0, 256],
    ]);

    let mut lattice = Lattice::from_generating_set(generators, &Z);
    println!("{lattice:?}");
    lattice.lll(&BigRational::new(99.into(), 100.into()), &Q, &Z);
    println!("{lattice:?}");
    let lattice = AffineLattice {
        offset: Vector::from_entries([1, 1, 0, 0, 0]),
        lattice,
    };

    let mut rng = rand::rng();
    let ur = BinaryBigInt::new(8);
    let v = OwnedVector::<Z>::from_iter(
        lattice.lattice.rank(),
        std::iter::repeat_with(|| ur.random(&mut rng).into()),
    );

    let sample = lattice.at(&v, &Z);
    println!("{sample:?}");
    let closest = lattice.lattice.cvp_rounding(sample.view(), &F64, &Z);
    println!("{closest:?}");
}

#[test]
fn babai_rounding_identity_dim_2() {
    use rand::random;
    let lattice = Lattice::from_basis(OwnedMatrix::<Z>::from_rows(&[[1, 0], [0, 1]]));

    for _ in 0..256 {
        let v = OwnedVector::<Z>::from_array([random::<u64>(), random::<u64>()]);
        let r = lattice.cvp_rounding(v.view(), &Q, &Z);
        assert_eq!(r, v);
    }
}

#[test]
fn babai_rounding_identity_dim_2_subspace() {
    use rand::random;
    let lattice = Lattice::from_basis(OwnedMatrix::<Z>::from_rows(&[[1, 0, 0], [0, 1, 0]]));

    for _ in 0..256 {
        let v = OwnedVector::<Z>::from_array([random::<u32>(), random::<u32>(), random::<u32>()]);
        let r = lattice.cvp_rounding(v.view(), &F64, &Z);
        assert_eq!(r.as_slice()[..2], v.as_slice()[..2]);
    }
}

#[test]
fn babai_rounding_linear_dim_3() {
    let lattice = Lattice::from_basis(OwnedMatrix::<Z>::from_rows(&[[3, 3, 3]]));

    assert_eq!(
        lattice.cvp_rounding_coeff(Vector::from_entries([2, 2, 2]).view(), &F64, &Z),
        [BigInt::from(1)]
    );
    assert_eq!(
        lattice.cvp_rounding_coeff(Vector::from_entries([2, -2, 0]).view(), &F64, &Z),
        [BigInt::from(0)]
    );
}
