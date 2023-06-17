use crate::select;
use crate::diophantine::hermite_normal_form;
use crate::{matrix::*, vector::*};
use rug::ops::NegAssign;
use rug::{Integer, Rational, Complete, Float};
use num_traits::{Zero, One, NumAssign};
use std::fmt::Debug;
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
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
    pub fn at<S: VectorStorage<Integer> + ?Sized>(
        &self, coefficients: &Vector<Integer, S>
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

        let mut rand = rug::rand::RandState::new();
        rand.seed(&rand::random::<u64>().into());

        let mut s = initial;
        for b in self.basis.rows() {
            let f = Integer::from(Integer::random_bits(bits, &mut rand));
            s += &(b * &f);
        }

        s.map_mut(|i| i.keep_signed_bits_mut(bits));
        s
    }

    /// Returns a random point on the lattice mod 2^bits.
    pub fn sample_point(&self, bits: u32) -> IOwnedVector {
        self.sample_point_impl(bits, Vector::zero(self.ambient_dim()))
    }

    /// Size reduce the basis.
    /// This essentially is Gram-Schmidt
    /// but rounding the coefficients to integers.
    fn size_reduce(&mut self) {
        for i in 0..self.basis.nrows() {
            for j in 0..i {
                let b_i = &self.basis[i];
                let b_j = &self.basis[j];
                let q = b_i.dot(b_j).div_rem_round(b_j.norm_sqr()).0;
                let s = b_j * &q;
                self.basis[i] -= &s;
            }
        }
    }

    /// Performs LLL basis reduction using rational numbers.
    pub fn lll(&mut self, delta: &Rational) {
        let n = self.basis.nrows();
        let mut swap_condition = true;

        while swap_condition {
            self.size_reduce();

            // Lovasz condition
            swap_condition = false;
            for i in 0..n-1 {
                let b = &self.basis[i];
                let c = &self.basis[i+1];

                let lhs = Rational::from(c.norm_sqr());

                let b_norm_sqr = b.norm_sqr();
                let q = Rational::from((c.dot(b), b_norm_sqr.clone()));
                let rhs = b_norm_sqr * (delta - q.square());

                if lhs < rhs {
                    self.basis.swap_rows(i, i + 1);
                    swap_condition = true;
                    break;
                }
            }
        }
    }

    pub fn cvp_planes_coeff<WT: WorkingType>(
        &self, t: &IVectorView, rad: Option<WT::Scalar>, ty: WT
    ) -> Option<IOwnedVector> {
        cvp_planes(&self.basis, t, rad, ty)
    }

    pub fn cvp_planes<WT: WorkingType>(
        &self, t: &IVectorView, rad: Option<WT::Scalar>, ty: WT
    ) -> Option<IOwnedVector> {
        Some(self.at(&self.cvp_planes_coeff(t, rad, ty)?))
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
    T: for<'a> Div<&'a T, Output = T>,
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
    T: for<'a> Div<&'a T, Output = T>,
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
    T: Clone + for<'a> Div<&'a T, Output = T>,
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
/// - `r` is an optional float that contains the maximum distance to search for.
/// The returned vector is the vector of coefficients.
/// This will always return some vector unless no vector is within `r` of the target.
/// This algorithm is a generalization Babai's nearest plane algorithm
/// that searches all planes that could contain the closest vector.
/// It is the simplest one I could think of.
pub fn cvp_planes<WT: WorkingType>(
    basis: &IOwnedMatrix, t: &IVectorView, rad: Option<WT::Scalar>, ty: WT
) -> Option<IOwnedVector> {
    assert!(basis.ncols() == t.dim(), "Mismatch of basis/target vector dimension.");
    // assert!(rad.as_ref().map_or(true, |rad| rad.prec() == prec),
    //     "rad needs to have the given precision.");
    let bf = basis.transform(|i| ty.from_int(i));

    // Q is the Gram-Schmidt orthonormalization and
    // R is the basis matrix with respect to the Gram-Schmidt basis.
    let (r, q) = rq_decomposition(&bf);
    assert!(r.ncols() == r.nrows());

    // Write the target vector in that basis.
    // If the vector not in the span of the basis,
    // then this will project it into the span.
    let qt = &q * &t.transform(|i| ty.from_int(i));

    // If a radius is given, we square it,
    // otherwise we use the squared norm of the target vector,
    // because the origin is a point on the lattice.
    // The actual value doesn't really matter,
    // it will get replaced with the distance
    // to Babai's nearest plane vector in the first descent.
    // We just do this to avoid having to pass around an Option.
    let rad = rad.map_or_else(|| qt.norm_sqr(), WT::square);

    /// Utility function for comparing a distance with the radius.
    /// If the radius is None, then this always accepts.
    fn in_radius<T: PartialOrd>(d: &T, rad: &T) -> bool {
        d.partial_cmp(rad).unwrap().is_le()
    }

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
            return in_radius(&d, rad).then(|| (
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
            offset.neg_assign();

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
            // either, by the way we iterate over the planes.
            if !in_radius(&d, &min_dist) {
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
            if in_radius(&d, &min_dist) {
                min = Some((v, w));
                min_dist = d;
            }
        }

        min
    }

    // Multiply the coefficients by the basis.
    cvp_impl(r.ncols() - 1, &r, qt.view(), &rad, ty).map(|v| v.0)
}

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
pub trait WorkingType: Copy {
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
        + PartialOrd;

    fn zero(self) -> Self::Scalar;

    fn is_zero(s: &Self::Scalar) -> bool;

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering;

    #[allow(clippy::wrong_self_convention)]
    fn from_int(self, i: &Integer) -> Self::Scalar;

    fn to_int(s: &Self::Scalar) -> Integer;

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

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.abs().partial_cmp(&b.abs()).unwrap()
    }

    fn from_int(self, i: &Integer) -> Self::Scalar {
        i.to_f64()
    }

    fn to_int(s: &Self::Scalar) -> Integer {
        Integer::from_f64(s.round()).unwrap()
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

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.abs().partial_cmp(&b.abs()).unwrap()
    }

    fn from_int(self, i: &Integer) -> Self::Scalar {
        i.to_f32()
    }

    fn to_int(s: &Self::Scalar) -> Integer {
        Integer::from_f32(s.round()).unwrap()
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
pub struct FP(u32);
impl WorkingType for FP {
    type Scalar = Float;
    fn zero(self) -> Self::Scalar {
        Float::new(self.0)
    }

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.cmp_abs(b).unwrap()
    }

    fn from_int(self, i: &Integer) -> Self::Scalar {
        Float::with_val(self.0, i)
    }

    fn to_int(s: &Self::Scalar) -> Integer {
        s.to_integer().unwrap()
    }

    fn round(s: Self::Scalar) -> Self::Scalar {
        s.round()
    }

    fn square(s: Self::Scalar) -> Self::Scalar {
        s.square()
    }

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar {
        s + i
    }
}

#[derive(Copy, Clone)]
pub struct Rat;
impl WorkingType for Rat {
    type Scalar = Rational;
    fn zero(self) -> Self::Scalar {
        Rational::new()
    }

    fn is_zero(s: &Self::Scalar) -> bool {
        s.is_zero()
    }

    fn cmp_abs(a: &Self::Scalar, b: &Self::Scalar) -> Ordering {
        a.cmp_abs(b)
    }

    fn from_int(self, i: &Integer) -> Self::Scalar {
        Rational::from(i)
    }

    fn to_int(s: &Self::Scalar) -> Integer {
        s.round_ref().into()
    }

    fn round(s: Self::Scalar) -> Self::Scalar {
        s.round()
    }

    fn square(s: Self::Scalar) -> Self::Scalar {
        s.square()
    }

    fn add_isize(s: Self::Scalar, i: isize) -> Self::Scalar {
        s + i
    }
}


#[test]
fn cvp_temp_test() {
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
    //let now = std::time::Instant::now();
    //for _ in 0..100 {
    //    let v = cvp_planes(&b, &t, 53, None);
    //}
    //println!("{}", now.elapsed().as_secs_f64());

    fn bench(name: &str, f: impl Fn() -> IOwnedVector, t: &IOwnedVector) {
        let now = std::time::Instant::now();
        let mut v = IOwnedVector::empty();
        for _ in 0..100 {
            v = f();
        }
        let time = now.elapsed().as_secs_f64() * 1000.;
        let d = (&v - t).norm_sqr();
        let d = Float::with_val(d.signed_bits(), d).sqrt();
        println!("{name}: {d:.3} in {time:.3}ms ({v:?})");
    }

    let now = std::time::Instant::now();
    bench("exact solution float", || l.cvp_planes(t.view(), None, FP(53)).unwrap(), &t);
    bench("exact solution f64", || l.cvp_planes(t.view(), None, F64).unwrap(), &t);
    bench("exact solution f32", || l.cvp_planes(t.view(), None, F32).unwrap(), &t);
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
    assert_eq!(l.cvp_planes(t.view(), None, FP(53)).unwrap(), [-30i32, 35, -98,
        61, -27, 75, -32, -3, 70, 8, 3, -77, -29, -103, 61, 58, -71, 41, 37, -40]);
}

#[test]
fn gram_schmidt_test() {
    let a = OwnedMatrix::<i32>::from_rows(&[
        [1, 2, 3],
        [3, 4, 5],
    ]);
    let a = a.to_float(53);
    let (r, q) = rq_decomposition(&a);
    println!("q: {:?}\nr: {:?}", q, r);
    println!("{:?}", &r * &q);
    println!("{:?}", &q * FVectorView::from_slice(&[
        Float::with_val(53, 5), Float::with_val(53, 6), Float::with_val(53, 7)
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
    lattice.lll(&Rational::from((99, 100)));
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
    use rand::random;
    let lattice = Lattice::from_basis(Matrix::from_rows(&[
        [3, 3, 3]
    ]));

    assert_eq!(lattice.cvp_rounding_coeff(
        Vector::from_entries([2, 2, 2]).view(), F64
    ), [1]);
    assert_eq!(lattice.cvp_rounding_coeff(
        Vector::from_entries([2, -2, 0]).view(), F64
    ), [0]);
}