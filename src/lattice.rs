//! Integer lattices and algorithms.

use crate::diophantine::hermite_normal_form;
use crate::{matrix::*, vector::*, keep_signed_bits_mut};
use num_traits::{Zero, ToPrimitive, FromPrimitive};
use num_bigint::{BigInt, RandBigInt};
use num_rational::BigRational;
use std::fmt::Debug;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign, Neg};
use std::cmp::Ordering;


/// An integer lattice.
#[derive(Debug)]
pub struct Lattice {
    /// The basis matrix of the lattice.
    /// The basis vectors are the rows of the matrix.
    pub basis: IOwnedMatrix,
}

impl Lattice {
    /// Creates an empty lattice.
    pub fn empty() -> Self {
        Self {
            basis: OwnedMatrix::empty()
        }
    }

    /// Is the lattice empty?
    pub fn is_empty(&self) -> bool {
        self.basis.is_empty()
    }

    /// The lattice basis are the rows of the matrix.
    pub fn from_basis(basis: IOwnedMatrix) -> Self {
        Self { basis }
    }

    /// The rows of the matrix generate the lattice
    /// but are potentially linearly dependent.
    /// This function will compute the Hermite normal form
    /// and remove zero rows.
    pub fn from_generating_set(mut generating_set: IOwnedMatrix) -> Self {
        hermite_normal_form(&mut generating_set);
        let rank = generating_set.nrows() - generating_set.rows().rev()
            .take_while(|r| r.iter().all(|i| i.is_zero()))
            .count();
        generating_set.shrink(rank);
        Self {
            basis: generating_set
        }
    }

    /// Returns the rank of the lattice, i.e. the number if basis vectors.
    pub fn rank(&self) -> usize {
        self.basis.nrows()
    }

    /// Returns the dimension of the ambient space.
    pub fn ambient_dim(&self) -> usize {
        self.basis.ncols()
    }

    /// Returns the vector on the lattice that is the linear
    /// combination of the basis vectors with the given coefficients.
    pub fn at<S: VectorStorage<BigInt> + ?Sized>(
        &self, coefficients: &Vector<BigInt, S>
    ) -> IOwnedVector {
        assert!(coefficients.dim() == self.rank());

        self.basis.rows().zip(coefficients)
            .fold(Vector::zero(self.ambient_dim()), |acc, (e, c)| acc + &(e * c))
    }

    /// Samples a point from the lattice that is added to the initial vector.
    pub(self) fn sample_point_impl(
        &self, bits: u32, initial: IOwnedVector
    ) -> IOwnedVector {
        assert!(!self.is_empty(), "Lattice is empty.");
        assert!(initial.dim() == self.ambient_dim());

        let rng = &mut rand::thread_rng();

        let mut s = initial;
        for b in self.basis.rows() {
            let f = rng.gen_bigint(bits as u64);
            s += &(b * &f);
        }

        s.map_mut(|i| keep_signed_bits_mut(i, bits));
        s
    }

    /// Returns a random point on the lattice mod 2^bits.
    pub fn sample_point(&self, bits: u32) -> IOwnedVector {
        self.sample_point_impl(bits, Vector::zero(self.ambient_dim()))
    }

    /// Size reduce the basis.
    /// This essentially is Gram-Schmidt
    /// but rounding the coefficients to integers.
    pub fn size_reduce(&mut self) {
        size_reduce(&mut self.basis);
    }

    /// Performs LLL basis reduction.
    pub fn lll<WT: WorkingType>(&mut self, delta: &WT::Scalar, ty: WT)
        where WT::Scalar: InnerProduct
    {
        lll(&mut self.basis, delta, ty);
    }

    pub fn cvp_planes_coeff<WT: WorkingType>(
        &self, t: &IVectorView, rad_sqr: Option<WT::Scalar>, ty: WT
    ) -> Option<IOwnedVector>
        where WT::Scalar: VectorNorm
    {
        cvp_planes(&self.basis, t, rad_sqr, ty)
    }

    pub fn cvp_planes<WT: WorkingType>(
        &self, t: &IVectorView, rad_sqr: Option<WT::Scalar>, ty: WT
    ) -> Option<IOwnedVector>
        where WT::Scalar: VectorNorm
    {
        Some(self.at(&self.cvp_planes_coeff(t, rad_sqr, ty)?))
    }

    pub fn cvp_rounding_coeff<WT: WorkingType>(
        &self, t: &IVectorView, ty: WT
    ) -> IOwnedVector {
        cvp_rounding(&self.basis, t, ty)
    }

    pub fn cvp_rounding<WT: WorkingType>(
        &self, t: &IVectorView, ty: WT
    ) -> IOwnedVector {
        self.at(&cvp_rounding(&self.basis, t, ty))
    }

    pub fn cvp_nearest_plane_coeff<WT: WorkingType>(
        &self, t: &IVectorView, ty: WT
    ) -> IOwnedVector {
        cvp_nearest_plane(&self.basis, t, true, ty).1
    }

    pub fn cvp_nearest_plane<WT: WorkingType>(
        &self, t: &IVectorView, ty: WT
    ) -> IOwnedVector {
        cvp_nearest_plane(&self.basis, t, false, ty).0
    }
}

/// A lattice that is offset from the origin by a vector.
/// Mathematicians would call this a lattice coset.
#[derive(Debug)]
pub struct AffineLattice {
    pub offset: IOwnedVector,
    pub lattice: Lattice,
}

impl AffineLattice {
    /// Creates an empty lattice.
    pub fn empty() -> Self {
        Self {
            offset: Vector::empty(),
            lattice: Lattice::empty(),
        }
    }

    /// Creates an affine lattice from an offset and a basis.
    pub fn from_offset_basis(offset: IOwnedVector, basis: IOwnedMatrix) -> Self {
        Self {
            offset,
            lattice: Lattice::from_basis(basis),
        }
    }

    /// Is this lattice empty?
    pub fn is_empty(&self) -> bool {
        self.offset.is_empty()
    }

    /// Returns a random point on the lattice mod 2^bits.
    pub fn sample_point(&self, bits: u32) -> IOwnedVector {
        self.lattice.sample_point_impl(bits, self.offset.clone())
    }

    /// Returns the shortest vector on the **affine(!)** lattice.
    /// This solves the CVP on the normal lattice with
    /// target offset.
    pub fn svp<WT: WorkingType>(
        &self, rad: Option<WT::Scalar>, ty: WT
    ) -> Option<IOwnedVector>
        where WT::Scalar: VectorNorm
    {
        self.lattice.cvp_planes(self.offset.view(), rad, ty)
            .map(|v| &self.offset - v)
    }
}

/// Size reduce the basis.
/// This essentially is Gram-Schmidt
/// but rounding the coefficients to integers.
pub fn size_reduce(basis: &mut IOwnedMatrix) {
    for i in 0..basis.nrows() {
        for j in 0..i {
            let b_i = &basis[i];
            let b_j = &basis[j];
            let dot = b_i.dot(b_j);
            let b_j_norm_sqr = b_j.norm_sqr();

            // Rounded div. rug has a function for this.
            let q: BigInt = if dot.sign() != b_j_norm_sqr.sign() {
                (dot + &b_j_norm_sqr / 2) / &b_j_norm_sqr
            } else {
                (dot - &b_j_norm_sqr / 2) / &b_j_norm_sqr
            };

            if q.is_zero() {
                continue;
            }
            let s = b_j * &q;
            basis[i] -= &s;
        }
    }
}

/// Performs LLL basis reduction using rational numbers.
pub fn lll<WT: WorkingType>(
    basis: &mut IOwnedMatrix,
    delta: &WT::Scalar,
    ty: WT
) where WT::Scalar: InnerProduct
{
    let n = basis.nrows();
    let mut swap_condition = true;

    while swap_condition {
        size_reduce(basis);

        // Lovasz condition
        swap_condition = false;
        for i in 0..n-1 {
            let b = &basis[i];
            let c = &basis[i+1];

            let lhs = ty.from_int(&c.norm_sqr());

            let b_norm_sqr = ty.from_int(&b.norm_sqr());
            let q = ty.from_int(&c.dot(b)) / &b_norm_sqr;
            let rhs = b_norm_sqr * &(-WT::square(q) + delta);

            if lhs < rhs {
                basis.swap_rows(i, i + 1);
                swap_condition = true;
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
pub fn cvp_rounding<WT: WorkingType>(
    basis: &IOwnedMatrix, t: &IVectorView, ty: WT
) -> IOwnedVector {
    let mut a = ty.zero_matrix(basis.ncols(), basis.nrows());

    for r in 0..a.nrows() {
        for c in 0..a.ncols() {
            a[(r, c)] = ty.from_int(&basis[(c, r)]);
        }
    }

    let b = OwnedVector::from_iter(
        t.dim(), t.iter().map(|i| ty.from_int(i))
    );

    // If the system has full rank,
    // then we can just use an exact solution.
    let x = if t.dim() == basis.nrows() {
        solve_linear(a, b, ty)
    }

    // Otherwise we can treat it as a Ordinary Least Squares problem,
    // i.e. find the point in the subspace spanned by the vectors
    // that is closest to the given point.
    else {
        let a_t = a.transpose();
        solve_linear(a_t * &a, a_t * &b, ty)
    };

    let x = x.expect("Basis vectors were not independent.");

    let mut res = IOwnedVector::zero(x.dim());
    for (i, f) in res.iter_mut().zip(x.iter()) {
        *i = WT::to_int(f);
    }

    res
}

/// Gram-Schmidt of the rows of the matrix `a`.
/// The rows will span the same subspace as the original rows,
/// even if they are not linearly independent.
/// The rows of the result are not normalized.
/// Call `gram_schmidt_orthonormal` if you want that.
pub fn gram_schmidt<T>(mut a: OwnedMatrix<T>) -> OwnedMatrix<T>
where
    for<'a> T: Clone + InnerProduct + SubAssign<&'a T>
        + Mul<&'a T, Output = T> + Div<&'a T, Output = T>,
{
    for i in 0..a.nrows() {
        for j in 0..i {
            let f = a[j].dot(&a[i]) / &a[j].norm_sqr();
            let p = &a[j] * &f;
            a[i] -= &p;
        }
    }

    a
}

/// Gram-Schmidt of the rows of the matrix `a`.
/// The rows of the result are normalized.
/// Call `gram_schmidt` if you don't want that.
/// This functions requires that the rows are linearly independent
/// (unlike `gram_schmidt`).
pub fn gram_schmidt_orthonormal<T>(mut a: OwnedMatrix<T>) -> OwnedMatrix<T>
where
    for<'a> T: Clone + InnerProduct + VectorNorm
        + Mul<&'a T, Output = T> + Div<&'a T, Output = T>
        + SubAssign<&'a T> + DivAssign<&'a T>,
{
    a = gram_schmidt(a);
    for i in 0..a.nrows() {
        a[i].normalize();
    }
    a
}

/// Computes the RQ-decomposition (row QR-decomposition) of a matrix.
/// Returns a lower triangular matrix `R` and an orthogonal matrix `Q`,
/// such that `R * Q = A`.
pub fn rq_decomposition<T>(a: &OwnedMatrix<T>) -> (OwnedMatrix<T>, OwnedMatrix<T>)
where
    for<'a> T: Clone + InnerProduct + VectorNorm
        + Mul<&'a T, Output = T> + Div<&'a T, Output = T>
        + SubAssign<&'a T> + DivAssign<&'a T>,
{
    let q = gram_schmidt_orthonormal(a.clone());
    let r = a * q.transpose();
    (r, q)
}

/// Approximates the CVP using Babai's nearest plane algorithm.
/// The return value is a pair of the vector and the vector of coefficients
/// of that vector if `compute_coefficients` is true.
/// Otherwise the second element of the pair is an empty vector.
pub fn cvp_nearest_plane<WT: WorkingType>(
    basis: &IOwnedMatrix,
    t: &IVectorView,
    compute_coefficients: bool,
    ty: WT
) -> (IOwnedVector, IOwnedVector) {
    let q = gram_schmidt(basis.transform(|i| ty.from_int(i)));
    let mut b = t.to_owned();
    let mut coeff = if compute_coefficients {
        OwnedVector::zero(basis.nrows())
    } else {
        OwnedVector::empty()
    };
    for i in (0..basis.nrows()).rev() {
        let a = &q[i];
        // c = round(⟨a, b⟩ / ⟨a, a⟩)
        let c = b.iter().zip(a.iter())
            .map(|(i, f)| ty.from_int(i) * f)
            .fold(ty.zero(), WT::Scalar::add);
        let c = c / &a.norm_sqr();
        let c = WT::to_int(&c);
        b -= &(&basis[i] * &c);
        if compute_coefficients {
            coeff[i] = c;
        }
    }

    (t - b, coeff)
}

/// Solves the CVP exactly.
/// - `basis` is the matrix that contains the basis as rows.
/// - `t` is the target vector.
/// - `prec` is the precision of the floating point numbers used.
/// - `rad_sqr` is an optional float that contains the square of the maximum
/// distance to search for.
///
/// The returned vector is the vector of coefficients.
/// This will always return some vector unless no vector is within `r` of the target.
/// This algorithm is a generalization Babai's nearest plane algorithm
/// that searches all planes that could contain the closest vector.
/// It is the simplest one I could think of.
pub fn cvp_planes<WT: WorkingType>(
    basis: &IOwnedMatrix, t: &IVectorView, rad_sqr: Option<WT::Scalar>, ty: WT
) -> Option<IOwnedVector>
where
    WT::Scalar: VectorNorm
{
    assert!(basis.ncols() == t.dim(),
        "Mismatch of basis/target vector dimension.");
    let bf = basis.transform(|i| ty.from_int(i));

    // Q is the Gram-Schmidt orthonormalization and
    // R is the basis matrix with respect to the Gram-Schmidt basis.
    let (r, q) = rq_decomposition(&bf);
    assert!(r.ncols() == r.nrows());

    // Write the target vector in that basis.
    // If the vector not in the span of the basis,
    // then this will project it into the span.
    let qt = &q * &t.transform(|i| ty.from_int(i));

    // We just do this to avoid having to pass around an Option.
    let rad = rad_sqr.unwrap_or_else(|| ty.infinity());

    // Multiply the coefficients of the vector.
    return cvp_impl(r.ncols() - 1, &r, qt.view(), &rad, ty).map(|v| v.0);

    /// This actually finds the closest point.
    /// `rad` is the squared norm.
    /// This returns the coordinates of the closest point in the basis
    /// as the first entry and the actual point as the second entry.
    fn cvp_impl<WT: WorkingType>(
        i: usize,
        r: &OwnedMatrix<WT::Scalar>,
        qt: &VectorView<WT::Scalar>,
        rad: &WT::Scalar,
        ty: WT
    ) -> Option<(IOwnedVector, OwnedVector<WT::Scalar>)> {
        // One dimensional lattice.
        if i == 0 {
            let qt = &qt[0];
            let r = &r[(0, 0)];
            // `m` is the index of the closest point.
            // `d` is the distance to it.
            let m = WT::round(qt.clone() / r);
            let plane = r.clone() * &m;
            let d = WT::square(plane.clone() - qt);
            return (&d <= rad).then(|| (
                IOwnedVector::from_array([WT::to_int(&m)]),
                OwnedVector::from_array([plane])
            ));
        }

        let qtc = &qt[i];
        let rc = &r[(i, i)];

        // Index of the closest plane before rounding.
        let start_index_fl = qtc.clone() / rc;
        let start_index = WT::round(start_index_fl.clone());

        // Suppose the start index was -0.4.
        // We would want to check the planes in the order
        // 0, -1, 1, -2, 2.
        // But if the start index was 0.4, we would want
        // 0, 1, -1, 2, -2
        // So depending on whether we round up or down
        // to the integer start index, we will negate
        // the offset from the start index.
        let negate_offset = (start_index_fl - &start_index) < ty.zero();

        // The current plane's offset.
        let mut offset = 0isize;

        let mut min_dist = rad.clone();
        let mut min = None;

        // Iterate over all possible planes.
        loop {
            // Index of the plane we are considering.
            let index = WT::add_isize(start_index.clone(), offset);

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
            let d = WT::square(index.clone() * rc - qtc);

            // If the plane is not in the radius,
            // then the next one in the loop definitely is not
            // either, because of the way we iterate over the planes.
            if d > min_dist {
                break;
            }

            // We can use a smaller radius inside the plane.
            // The target to its projection to any point
            // in the plane form a right triangle.
            // So by Pythagoras the squared distance in the
            // plane can only be <= rad - d.
            // rad and d are already the square of the distance.
            let plane_dist = min_dist.clone() - d;

            // Compute the point in the plane we need to be close to now.
            let point_in_plane = qt - &r[i] * &index;

            // Recursively find the closest point.
            let Some((mut v, mut w)) = cvp_impl(
                i - 1, r, point_in_plane.view(), &plane_dist, ty
            ) else {
                continue
            };
            assert!(v.dim() == i);

            // v is the new coordinate vector of the point.
            // It is the coordinates of the closest point of
            // the previous call, concat the index of the plane.
            v.append(WT::to_int(&index));

            // w is the closest point.
            // It is the closest point of the previous call,
            // plus the index of the plane times the basis vector
            // of the current iteration.
            w += &(&r[i][0..i] * &index);
            w.append(index.clone() * rc);

            // Compute the distance to the point.
            let d = (&w - &qt[0..i+1]).norm_sqr();

            // If the distance is smaller than the current minimal dist,
            // then we have found the new best point.
            if d <= min_dist {
                min = Some((v, w));
                min_dist = d;
            }
        }

        min
    }
}

/// Solves a square system of linear equations.
fn solve_linear<WT: WorkingType>(
    mut a: OwnedMatrix<WT::Scalar>, mut b: OwnedVector<WT::Scalar>, ty: WT
) -> Option<OwnedVector<WT::Scalar>> {
    assert!(a.nrows() == a.ncols(),
        "This function only supports non-singular square systems.");
    for i in 0..a.ncols() {
        // Choose a pivot in the c-th column.
        let pivot = a.col(i)
            .iter()
            .enumerate()
            .skip(i)
            .filter(|e| !WT::is_zero(e.1))
            .max_by(|e, f| WT::cmp_abs(e.1, f.1))?.0;
        // Swap the pivot row with the current row.
        a.swap_rows(pivot, i);
        let pivot = a[(i, i)].clone();
        for r in i+1..a.nrows() {
            let fac = a[(r, i)].clone() / &pivot;
            for c in i+1..a.ncols() {
                let s = fac.clone() * &a[(i, c)];
                a[(r, c)] -= &s;
            }
            let s = fac * &b[i];
            b[r] -= &s;
        }
    }

    let mut result = ty.zero_vector(a.ncols());
    for i in (0..a.ncols()).rev() {
        let mut sum = b[i].clone();
        for j in i+1..a.ncols() {
            sum -= &(a[(i, j)].clone() * &result[j]);
        }

        sum /= &a[(i, i)];
        result[i] = sum;
    }

    Some(result)
}

/// Type that is used internally by lattice algorithms.
/// This is an ugly trait that is just the collection
/// of things that are needed by the algorithms.
/// Feel free to extend it for your own needs.
pub trait WorkingType: Copy
{
    type Scalar: Clone + Debug
        + Add<Self::Scalar, Output = Self::Scalar>
        + Sub<Self::Scalar, Output = Self::Scalar>
        + for<'a> Add<&'a Self::Scalar, Output = Self::Scalar>
        + for<'a> Sub<&'a Self::Scalar, Output = Self::Scalar>
        + for<'a> Mul<&'a Self::Scalar, Output = Self::Scalar>
        + for<'a> Div<&'a Self::Scalar, Output = Self::Scalar>
        + for<'a> AddAssign<&'a Self::Scalar>
        + for<'a> SubAssign<&'a Self::Scalar>
        + for<'a> MulAssign<&'a Self::Scalar>
        + for<'a> DivAssign<&'a Self::Scalar>
        + Neg<Output = Self::Scalar>
        + PartialOrd
        + InnerProduct;

    fn zero(self) -> Self::Scalar;

    fn eps(self) -> Self::Scalar;

    fn infinity(self) -> Self::Scalar {
        panic!("`infinity` is not supported for this type.");
    }

    fn is_zero(s: &Self::Scalar) -> bool;

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering;

    #[allow(clippy::wrong_self_convention)]
    fn from_int(self, i: &BigInt) -> Self::Scalar;

    fn to_int(s: &Self::Scalar) -> BigInt;

    fn round(s: Self::Scalar) -> Self::Scalar;

    fn zero_vector(self, dim: usize) -> OwnedVector<Self::Scalar> {
        OwnedVector::from_iter(dim, std::iter::repeat_with(|| self.zero()))
    }

    fn zero_matrix(self, r: usize, c: usize) -> OwnedMatrix<Self::Scalar> {
        OwnedMatrix::from_iter(r, c, std::iter::repeat_with(|| self.zero()))
    }

    fn square(s: Self::Scalar) -> Self::Scalar;

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar;
}

#[derive(Copy, Clone)]
pub struct F64;
impl WorkingType for F64 {
    type Scalar = f64;
    fn zero(self) -> Self::Scalar {
        0.0
    }

    fn eps(self) -> Self::Scalar {
        f64::EPSILON
    }

    fn infinity(self) -> Self::Scalar {
        f64::INFINITY
    }

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.abs().partial_cmp(&b.abs()).unwrap()
    }

    fn from_int(self, i: &BigInt) -> Self::Scalar {
        i.to_f64().unwrap()
    }

    fn to_int(s: &Self::Scalar) -> BigInt {
        BigInt::from_f64(s.round()).unwrap()
    }

    fn round(s: Self::Scalar) -> Self::Scalar {
        s.round()
    }

    fn square(s: Self::Scalar) -> Self::Scalar {
        s * s
    }

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar {
        s + i as Self::Scalar
    }
}

#[derive(Copy, Clone)]
pub struct F32;
impl WorkingType for F32 {
    type Scalar = f32;
    fn zero(self) -> Self::Scalar {
        0.0
    }

    fn eps(self) -> Self::Scalar {
        f32::EPSILON
    }

    fn infinity(self) -> Self::Scalar {
        f32::INFINITY
    }

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.abs().partial_cmp(&b.abs()).unwrap()
    }

    fn from_int(self, i: &BigInt) -> Self::Scalar {
        i.to_f32().unwrap()
    }

    fn to_int(s: &Self::Scalar) -> BigInt {
        BigInt::from_f32(s.round()).unwrap()
    }

    fn round(s: Self::Scalar) -> Self::Scalar {
        s.round()
    }

    fn square(s: Self::Scalar) -> Self::Scalar {
        s * s
    }

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar {
        s + i as Self::Scalar
    }
}

#[derive(Copy, Clone)]
pub struct Rat;
impl WorkingType for Rat {
    type Scalar = BigRational;
    fn zero(self) -> Self::Scalar {
        BigRational::zero()
    }

    fn eps(self) -> Self::Scalar {
        BigRational::zero()
    }

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        //use num_traits::Signed;
        //a.abs().cmp(&b.abs())
        (a.numer().magnitude() * b.denom().magnitude()).cmp(
            &(b.numer().magnitude() * a.denom().magnitude()))
    }

    fn from_int(self, i: &BigInt) -> Self::Scalar {
        BigRational::from_integer(i.clone())
    }

    fn to_int(s: &Self::Scalar) -> BigInt {
        s.round().to_integer()
    }

    fn round(s: Self::Scalar) -> Self::Scalar {
        s.round()
    }

    fn square(s: Self::Scalar) -> Self::Scalar {
        &s * &s
    }

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar {
        s + BigRational::from_isize(i).unwrap()
    }
}

#[test]
fn cvp_exact_dim20() {
    let l = Lattice::from_basis(IOwnedMatrix::from_rows(&[
        [  1,  -74,   20,   19,   -5,   21,  -19,   18,  -54,  -56,  -40,
          -38,  -58,   54,  -34,   -1,   -3,    0,  -27,    8],
        [-44,   19,  -20,   14,  -34,  -61,   53,  -31,   42,   42,   27,
           40,  -77,  -59,    2,  -26,    9,    3,  -49,   33],
        [-19,  -12,   18,  -31,   63,    8,   52,   52,  -50,  -19,   22,
           -1,   51,   -8,  -57,   70,  -62,  -45,   48,  -51],
        [-66,   -5,   15,   34,    2,  -50,   82,   28,   16,   18,   17,
           66,   23,  -38,   39,  -36,  -66,   19,  -64,  -62],
        [-15,   10,  -65,  -46,  -55,  -22,  -88,  -49,   46,  -19,   -4,
           -6,   -5,  -20,   23,  -24,  -86,  -29,   -7,   16],
        [ 35,  -57,    9,   41,    9,   -7,   63,   11,  -72,  -23,   -2,
         -103,   62,   -9,   64,    1,   25,   10,   48,  -41],
        [ 18,   11,   31,    2,  -51,   48,   11,  -33,   69,   93,  -12,
          -52,   -3, -100,   14,  -15,   85,  -61,    6,   27],
        [-58,  -54,  -29,  -28,   93,  -88,    3,  -63,   -3,   -9,   26,
           22,   89,   26,   11,   -2,  -21,   -5,   37,  -36],
        [ 24,   72,  -79,  -38,   50,   17,  -54,  -24,   38,   -9,  -64,
          -32,  -10,   70,  -67,  -88,   -4,  -34,   -4,   -3],
        [-24,  -54,   11,   34,  -14,   -5, -132,   46,   51,   67,   24,
           -3,  -10,   -6,   22,   38,  -15,  -23,   29,  -39],
        [ 72,    9,   49,  -19,   66,   -6,  -53,  -77,   40,   -9,  -52,
          -20,   90,  -12,   58, -107,  -47,  -64,  -66,  -10],
        [ 48,  -28,   41,  -76,   17,  -13,  -20,  -16,  -15,   75,  -30,
           51,  -55,   31,  -50,    3,  -60,  -34,   13,   31],
        [ 62,   -9,   13,  -10,   39,   50,   81,   94,  -38,    7,  -62,
          -49,   34,   61,   45,   30,  -51,  -78,  -70,   -4],
        [ 34,   92,  -16,   -3,  113,  -24,    1,   40,  -30,   91,  -57,
           -6,  -28,   -7,   13,   64,   18,   24,  -33,  -10],
        [ -1,  -20,  -45,   44,  -75,   31,   -9,  -47,   74,   -7,   64,
           77,  -41,    9,   52,  -33,  -83,  118,   63,    1],
        [ -9,  -37,   11, -106,  -13,    1,   74,   11,   89,  -10,   61,
           43,  -17,  -45,   -7,    5, -103,   43,  -36,   46],
        [ 60, -101,  -48,  -23,   37,  -45,    1,   23,   52,  -18,  -19,
          -36, -125,  -23,   22,  -32,  -28,  -41,   16,  -34],
        [ 19,  -47,  -85,   17,    2,    1,   12,   19,   27,  -21,   43,
           43,   -9,    8,  -60,   36,   83,   65,  -50,   58],
        [ 59,   -4,   51,   36,  -10,  -12,  -19,  -54,   93,   12,   23,
           31,   77,   18,   45,   57,   46,   52,  -21,  -51],
        [-57,   49,   26,   10,  -22,   32,  -52,  -71,   -9,   31,   45,
           28,  -28,  -30,   24,   44,   88,    2,  -63, -105],
    ]));
    let t = IOwnedVector::from_entries([-48, 69, -76, 36, -72, 31, -53,
        -7, 54, 74, 6, -82, -13, -32, 7, 53, -60, -44, 38, -97]);
    assert_eq!(l.cvp_planes(t.view(), None, F64).unwrap(), [-30i32, 35, -98,
        61, -27, 75, -32, -3, 70, 8, 3, -77, -29, -103, 61, 58, -71, 41, 37, -40]
        .iter().map(|i| BigInt::from(*i)).collect::<Vec<_>>());
}

#[test]
fn gram_schmidt_test() {
    let a = OwnedMatrix::<f64>::from_rows(&[
        [1., 2., 3.],
        [3., 4., 5.],
    ]);
    let (r, q) = rq_decomposition(&a);
    println!("q: {:?}\nr: {:?}", q, r);
    println!("{:?}", &r * &q);
    println!("{:?}", &q * VectorView::from_slice(&[
        5., 6., 7.
    ]));
}

#[test]
fn nearest_plane_example() {
    let l = Lattice::from_basis(IOwnedMatrix::from_rows(&[
        [2, 3, 1],
        [4, 1, -3],
        [2, 2, 2],
    ]));

    let t = IOwnedVector::from_entries([4, 2, 7]);
    let c = l.cvp_nearest_plane(t.view(), F64);
    println!("{:?}", c);
}

#[test]
fn babai_rounding_example() {
    let gen = Matrix::from_rows(&[
        [-97, 75, -97, 75, 22],
        [101, 38, 101, 38, 117],
        [256, 0, 0, 0, 0],
        [0, 256, 0, 0, 0],
        [0, 0, 256, 0, 0],
        [0, 0, 0, 256, 0],
        [0, 0, 0, 0, 256],
    ]);

    let mut lattice = Lattice::from_generating_set(gen);
    println!("{:?}", lattice);
    lattice.lll(&BigRational::new(99.into(), 100.into()), Rat);
    println!("{:?}", lattice);
    let lattice = AffineLattice {
        offset: Vector::from_entries([1, 1, 0, 0, 0]),
        lattice,
    };

    let sample = lattice.sample_point(8);
    println!("{:?}", sample);
    let closest = lattice.lattice.cvp_rounding(sample.view(), F64);
    println!("{:?}", closest);
}

#[test]
fn babai_rounding_identity_dim_2() {
    use rand::random;
    let lattice = Lattice::from_basis(Matrix::from_rows(&[
        [1, 0],
        [0, 1],
    ]));

    for _ in 0..256 {
        let v = Vector::from_array([random::<u64>(), random::<u64>()]);
        let r = lattice.cvp_rounding(v.view(), Rat);
        assert_eq!(r, v);
    }
}

#[test]
fn babai_rounding_identity_dim_2_subspace() {
    use rand::random;
    let lattice = Lattice::from_basis(Matrix::from_rows(&[
        [1, 0, 0],
        [0, 1, 0],
    ]));

    for _ in 0..256 {
        let v = Vector::from_array(
            [random::<u32>(), random::<u32>(), random::<u32>()]
        );
        let r = lattice.cvp_rounding(v.view(), F64);
        assert_eq!(r.as_slice()[..2], v.as_slice()[..2]);
    }
}

#[test]
fn babai_rounding_linear_dim_3() {
    let lattice = Lattice::from_basis(Matrix::from_rows(&[
        [3, 3, 3]
    ]));

    assert_eq!(lattice.cvp_rounding_coeff(
        Vector::from_entries([2, 2, 2]).view(), F64
    ), [BigInt::from(1)]);
    assert_eq!(lattice.cvp_rounding_coeff(
        Vector::from_entries([2, -2, 0]).view(), F64
    ), [BigInt::from(0)]);
}