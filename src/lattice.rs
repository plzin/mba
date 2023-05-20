use crate::diophantine::hermite_normal_form;
use crate::{matrix::*, vector::*};
use nalgebra::{DMatrix, DVector, Scalar};
use rug::{Integer, Rational, Complete, Float};
use num_traits::{Zero, One, NumAssign};
use std::ops::{Add, Sub, Mul, Div, AddAssign, SubAssign, MulAssign, DivAssign};
use std::cmp::Ordering;

/// An integer lattice.
#[derive(Debug)]
pub struct Lattice {
    /// The basis matrix of the lattice.
    /// The basis vectors are the rows of the matrix.
    pub basis: IMatrix,
}

impl Lattice {
    /// Creates an empty lattice.
    pub fn empty() -> Self {
        Self {
            basis: Matrix::empty()
        }
    }

    /// Is the lattice empty?
    pub fn is_empty(&self) -> bool {
        self.basis.is_empty()
    }

    /// The lattice basis are the rows of the matrix.
    pub fn from_basis(basis: IMatrix) -> Self {
        Self { basis }
    }

    /// The rows of the matrix generate the lattice
    /// but are potentially linearly dependent.
    /// This function will compute the Hermite normal form
    /// and remove zero rows.
    pub fn from_generating_set(mut generating_set: IMatrix) -> Self {
        hermite_normal_form(&mut generating_set);
        let rank = generating_set.rows - generating_set.rows().rev()
            .take_while(|r| r.iter().all(|i| i.is_zero()))
            .count();
        generating_set.shrink(rank);
        Self {
            basis: generating_set
        }
    }

    /// Returns the rank of the lattice, i.e. the number if basis vectors.
    pub fn rank(&self) -> usize {
        self.basis.rows
    }

    /// Returns the dimension of the ambient space.
    pub fn ambient_dim(&self) -> usize {
        self.basis.cols
    }
    
    /// Returns the vector on the lattice that is the linear
    /// combination of the basis vectors with the given coefficients.
    pub fn at<T: AsRef<[Integer]>>(&self, coefficients: T) -> IVector {
        let coefficients = coefficients.as_ref();
        assert!(coefficients.len() == self.rank());

        self.basis.rows().zip(coefficients)
            .fold(Vector::zero(self.ambient_dim()), |acc, (e, c)| acc + &(c * e))
    }

    /// Samples a point from the lattice that is added to the initial vector.
    pub(self) fn sample_point_impl(
        &self, bits: u32, initial: IVector
    ) -> IVector {
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
    pub fn sample_point(&self, bits: u32) -> IVector {
        self.sample_point_impl(bits, Vector::zero(self.ambient_dim()))
    }

    /// Approximate the CVP using Babai's rounding technique with floats,
    /// but returns the coefficients of the linear combination of the basis vectors.
    pub fn cvp_rounding_coeff(&self, v: &IVector) -> IVector {
        cvp_rounding_float(&self.basis, v)
    }

    /// Approximate the CVP using Babai's rounding technique with floats.
    pub fn cvp_rounding(&self, v: &IVector) -> IVector {
        self.at(self.cvp_rounding_coeff(v))
    }

    /// Approximate the CVP using Babai's rounding technique with rational numbers.
    pub fn cvp_rounding_coeff_exact(&self, v: &IVector) -> IVector {
        cvp_rounding_exact(&self.basis, v)
    }

    /// Approximate the CVP using Babai's rounding technique with rational numbers.
    pub fn cvp_rounding_exact(&self, v: &IVector) -> IVector {
        self.at(self.cvp_rounding_coeff_exact(v))
    }

    /// Size reduce the basis.
    /// This essentially is Gram-Schmidt
    /// but rounding the coefficients to integers.
    fn size_reduce(&mut self) {
        for i in 0..self.basis.rows {
            for j in 0..i {
                let b_i = &self.basis[i];
                let b_j = &self.basis[j];
                let q = b_i.dot(&b_j).div_rem_round(b_j.norm_sqr()).0;
                let s = b_j * &q;
                self.basis[i] -= &s;
            }
        }
    }

    /// Performs LLL basis reduction using rational numbers.
    pub fn lll(&mut self, delta: &Rational) {
        let n = self.basis.rows;
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
}


/// A lattice that is offset from the origin by a vector.
/// Mathematicians would call this a lattice coset.
#[derive(Debug)]
pub struct AffineLattice {
    pub offset: IVector,
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
    pub fn from_offset_basis(offset: IVector, basis: IMatrix) -> Self {
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
    pub fn sample_point(&self, bits: u32) -> IVector {
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
pub fn cvp_rounding<T: Field>(basis: &IMatrix, v: &IVector) -> IVector {
    use nalgebra::{DMatrix, DVector, LU};
    let mut a = DMatrix::<T>::zeros(
        v.dim, basis.rows
    );
    for c in 0..a.ncols() {
        for r in 0..a.nrows() {
            a[(r, c)] = T::from_int(&basis[c][r]);
        }
    }

    let b = DVector::<T>::from_iterator(v.dim, v.iter().map(T::from_int));

    // If the system has full rank,
    // then we can just use an exact solution.
    let x = if v.dim == basis.rows {
        T::solve_linear(a, b)
    }
    
    // Otherwise we can treat it as a Ordinary Least Squares problem,
    // i.e. find the point in the subspace spanned by the vectors
    // that is closest to the given point.
    else {
        let a_t = a.transpose();
        let b_new = &a_t * &b;
        T::solve_linear(a_t * a, b_new)
    };

    let x = x.expect("Basis vectors were not independent.");

    let mut res = Vector::zero(x.nrows());
    for (i, f) in res.iter_mut().zip(x.iter()) {
        *i = f.to_int();
    }

    res
}

/// Approximates the CVP using Babai's rounding technique.
/// The returns value is a vector of the coefficients
/// of the linear combination of the basis vectors,
/// so if you want the actual point,
/// you still have to matrix multiply with the basis matrix.
/// In practice, it is a good idea to reduce the basis (e.g. using LLL)
/// so that the approximation is good.
/// 
/// This uses 64-bit floating point arithmetic instead of
/// exact rational arithmetic.
/// Make sure the entries in the basis fit into 64-bit floats,
/// but even then, this can lead to inaccuracies.
/// Given the standard basis this will fail to return
/// the correct coefficients when the entries in v
/// can't be represented by 64-bit floats.
pub fn cvp_rounding_float(basis: &IMatrix, v: &IVector) -> IVector {
    cvp_rounding::<f64>(basis, v)
}

/// Approximates the CVP using Babai's rounding technique.
/// The returns value is a vector of the coefficients
/// of the linear combination of the basis vectors,
/// so if you want the actual point,
/// you still have to matrix multiply with the basis matrix.
/// In practice, it is a good idea to reduce the basis (e.g. using LLL)
/// so that the approximation is good.
pub fn cvp_rounding_exact(basis: &IMatrix, v: &IVector) -> IVector {
    cvp_rounding::<Rational>(basis, v)
}

#[test]
fn babai_rounding_example() {
    let gen = Matrix::from_rows(&[
        IVector::from_array([-97, 75, -97, 75, 22]),
        IVector::from_array([101, 38, 101, 38, 117]),
        IVector::from_array([256, 0, 0, 0, 0]),
        IVector::from_array([0, 256, 0, 0, 0]),
        IVector::from_array([0, 0, 256, 0, 0]),
        IVector::from_array([0, 0, 0, 256, 0]),
        IVector::from_array([0, 0, 0, 0, 256]),
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
    let closest = lattice.lattice.cvp_rounding(&sample);
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
        let r = lattice.cvp_rounding_exact(&v);
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
        let r = lattice.cvp_rounding(&v);
        assert_eq!(r.as_slice()[..2], v.as_slice()[..2]);
    }
}

#[test]
fn babai_rounding_linear_dim_3() {
    use rand::random;
    let lattice = Lattice::from_basis(Matrix::from_rows(&[
        IVector::from_array([3, 3, 3])
    ]));

    assert_eq!(
        lattice.cvp_rounding_coeff(&Vector::from_array([2, 2, 2])).as_slice(),
        &[1]);
    assert_eq!(
        lattice.cvp_rounding_coeff(&Vector::from_array([2, -2, 0])).as_slice(),
        &[0]);
}

/// Field for linear algebra.
/// This will either be floating point numbers or rug::Rational.
pub trait Field: 'static + Clone + PartialEq + Zero + One
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<Self, Output = Self>
    + Div<Self, Output = Self>
    + AddAssign<Self>
    + SubAssign<Self>
    + MulAssign<Self>
    + DivAssign<Self>
    + std::fmt::Debug
{
    /// Converts an integer to this type.
    fn from_int(i: &Integer) -> Self;

    /// Converts this type to an integer.
    fn to_int(&self) -> Integer;

    /// Solve a linear system of equations with coefficients in the field.
    fn solve_linear(a: DMatrix<Self>, b: DVector<Self>) -> Option<DVector<Self>>;
}

impl Field for Rational {
    fn from_int(i: &Integer) -> Self {
        Rational::from(i)
    }

    fn to_int(&self) -> Integer {
        self.round_ref().into()
    }

    fn solve_linear(a: DMatrix<Self>, b: DVector<Self>) -> Option<DVector<Self>> {
        solve_linear_rat(a, b) 
    }
}

impl Field for f32 {
    fn from_int(i: &Integer) -> Self {
        i.to_f32()
    }

    fn to_int(&self) -> Integer {
        Integer::from_f32(self.round()).unwrap()
    }

    fn solve_linear(a: DMatrix<Self>, b: DVector<Self>) -> Option<DVector<Self>> {
        solve_linear(a, b)
    }
}

impl Field for f64 {
    fn from_int(i: &Integer) -> Self {
        i.to_f64()        
    }

    fn to_int(&self) -> Integer {
        Integer::from_f64(self.round()).unwrap()
    }

    fn solve_linear(a: DMatrix<Self>, b: DVector<Self>) -> Option<DVector<Self>> {
        solve_linear(a, b)
    }
}


fn solve_linear_rat(
    mut a: DMatrix<Rational>, mut b: DVector<Rational>
) -> Option<DVector<Rational>> {
    assert!(a.nrows() == a.ncols(),
        "This function only supports non-singular square systems.");
    for i in 0..a.ncols() {
        // Choose a pivot in the c-th column.
        let pivot = a.column(i)
            .iter()
            .enumerate()
            .skip(i)
            .filter(|e| e.1 != &Rational::new())
            .max_by(|e, f| e.1.cmp_abs(f.1))?.0;
        a.swap_rows(pivot, i);
        let pivot = a[(i, i)].clone();
        for r in i+1..a.nrows() {
            let fac = (&a[(r, i)] / &pivot).complete();
            for c in i+1..a.ncols() {
                let s = (&fac * &a[(i, c)]).complete();
                a[(r, c)] -= s;
            }
            let s = fac * &b[i];
            b[r] -= s;
        }
    }

    let mut result = DVector::<Rational>::zeros(a.ncols());
    for i in (0..a.ncols()).rev() {
        let mut sum = b[i].clone();
        for j in i+1..a.ncols() {
            sum -= (&a[(i, j)] * &result[j]).complete();
        }

        sum /= &a[(i, i)];
        result[i] = sum;
    }

    Some(result)
}

/// PartialOrd should not return None for any of the elements in the matrix.
/// We can't use Ord because of the floating point types.
fn solve_linear<T: NumAssign + PartialOrd + Copy + Scalar>(
    mut a: DMatrix<T>, mut b: DVector<T>
) -> Option<DVector<T>> {
    assert!(a.nrows() == a.ncols(),
        "This function only supports non-singular square systems.");
    let abs = |x: T| if x < T::zero() { T::zero() - x } else { x };
    for i in 0..a.ncols() {
        // Choose a pivot in the c-th column.
        let pivot = a.column(i)
            .iter()
            .enumerate()
            .skip(i)
            .filter(|e| *e.1 != T::zero())
            .max_by(|e, f| abs(*e.1).partial_cmp(&abs(*f.1)).unwrap())?.0;
        a.swap_rows(pivot, i);
        let pivot = a[(i, i)].clone();
        for r in i+1..a.nrows() {
            let fac = a[(r, i)] / pivot;
            for c in i+1..a.ncols() {
                let s = fac * a[(i, c)];
                a[(r, c)] -= s;
            }
            let s = b[i].clone() * fac;
            b[r] -= s;
        }
    }

    let mut result = DVector::zeros(a.ncols());
    for i in (0..a.ncols()).rev() {
        let mut sum = b[i].clone();
        for j in i+1..a.ncols() {
            sum -= a[(i, j)] * result[j];
        }

        sum /= a[(i, i)];
        result[i] = sum;
    }

    Some(result)
}
