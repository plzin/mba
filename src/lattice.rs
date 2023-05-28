use crate::diophantine::hermite_normal_form;
use crate::{matrix::*, vector::*};
use rug::ops::NegAssign;
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
    pub fn cvp_rounding_coeff(&self, v: &IVV) -> IVector {
        cvp_rounding_float(&self.basis, v)
    }

    /// Approximate the CVP using Babai's rounding technique with floats.
    pub fn cvp_rounding(&self, v: &IVV) -> IVector {
        self.at(self.cvp_rounding_coeff(v))
    }

    /// Approximate the CVP using Babai's rounding technique with rational numbers.
    pub fn cvp_rounding_coeff_exact(&self, v: &IVV) -> IVector {
        cvp_rounding_exact(&self.basis, v)
    }

    /// Approximate the CVP using Babai's rounding technique with rational numbers.
    pub fn cvp_rounding_exact(&self, v: &IVV) -> IVector {
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
                let q = b_i.dot(b_j).div_rem_round(b_j.norm_sqr()).0;
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
pub fn cvp_rounding<T: Field>(basis: &'_ IMatrix, v: &'_ IVV) -> IVector
    where for<'a, 'b, 'c> &'a Matrix<T>: Mul<&'b Matrix<T>, Output = Matrix<T>>
        + Mul<&'c VV<T>, Output = Vector<T>>,
{
    let mut a = Matrix::<T>::zero(
        basis.cols, basis.rows
    );
    for r in 0..a.rows {
        for c in 0..a.cols {
            a[(r, c)] = T::from_int(&basis[(c, r)]);
        }
    }

    let b = Vector::<T>::from_iter(v.dim(), v.iter().map(T::from_int));

    // If the system has full rank,
    // then we can just use an exact solution.
    let x = if v.dim() == basis.rows {
        T::solve_linear(a, b)
    }
    
    // Otherwise we can treat it as a Ordinary Least Squares problem,
    // i.e. find the point in the subspace spanned by the vectors
    // that is closest to the given point.
    else {
        let a_t = a.transposed();
        T::solve_linear(&a_t * &a, &a_t * b.view())
    };

    let x = x.expect("Basis vectors were not independent.");

    let mut res = Vector::zero(x.dim);
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
pub fn cvp_rounding_float(basis: &IMatrix, v: &IVV) -> IVector {
    cvp_rounding::<f64>(basis, v)
}

/// Approximates the CVP using Babai's rounding technique.
/// The returns value is a vector of the coefficients
/// of the linear combination of the basis vectors,
/// so if you want the actual point,
/// you still have to matrix multiply with the basis matrix.
/// In practice, it is a good idea to reduce the basis (e.g. using LLL)
/// so that the approximation is good.
pub fn cvp_rounding_exact(basis: &IMatrix, v: &IVV) -> IVector {
    cvp_rounding::<Rational>(basis, v)
}

/// Approximates the CVP using Babai's nearest plane algorithm.
pub fn cvp_nearest_plane_float(basis: &IMatrix, v: &IVV, prec: u32) -> IVector {
    let q = gram_schmidt(basis.to_float(prec), false);
    let mut b = v.to_owned();
    for i in (0..basis.rows).rev() {
        let a = &q[i];
        // c = round(⟨a, b⟩ / ⟨a, a⟩)
        let c = b.iter().zip(a.iter())
            .map(|(i, f)| Float::with_val(prec, i) * f)
            .fold(Float::with_val(prec, 0), |acc, f| acc + f);
        let c = c / a.norm_sqr();
        let c = c.to_integer().unwrap();
        b -= &(&c * &basis[i]);
    }

    v - b
}

/// Solves the CVP exactly.
/// - `basis` is the matrix that contains the basis as rows.
/// - `t` is the target vector.
/// - `prec` is the precision of the floating point numbers used.
/// - `r` is an optional float that contains the maximum distance to search for.
/// This will always return some vector unless no vector is within `r` of the target.
/// This algorithm is a generalization Babai's nearest plane algorithm
/// that searches all planes that could contain the closest vector.
/// It is the simplest one I could think of.
pub fn cvp_planes(basis: &IMatrix, t: &IVV, prec: u32, mut rad: Option<Float>) -> Option<IVector> {
    assert!(basis.cols == t.dim(), "Mismatch of basis/target vector dimension.");
    assert!(rad.as_ref().map_or(true, |rad| rad.prec() == prec),
        "rad needs to have the given precision.");
    let bf = basis.to_float(prec);

    // Q is the Gram-Schmidt orthonormalization and
    // R is the basis matrix with respect to the Gram-Schmidt basis.
    let (r, q) = rq_decomposition(&bf);
    assert!(r.cols == r.rows);

    // Write the target vector in that basis.
    // If the vector not in the span of the basis,
    // then this will project it into the span.
    let qt = &q * &t.to_float(prec);
    //println!("Lattice in the new basis: {:?}", r);
    //println!("Target vector in the new basis: {:?}", qt);

    rad = rad.map(|f| f.square());

    /// Utility function for comparing a distance with the radius.
    /// If the radius is None, then this always accepts.
    fn in_radius(d: &Float, rad: &Option<Float>) -> bool {
        rad.as_ref().map_or(true, |rad| d.total_cmp(rad).is_le())
    }

    /// This actually finds the closest point.
    /// `rad` is the squared norm.
    fn cvp_impl(i: usize, r: &FMatrix, qt: &FVV, prec: u32, rad: &Option<Float>) -> Option<IVector> {
        // One dimensional lattice.
        if i == 0 {
            let qt = &qt[0];
            let r = &r[(0, 0)];
            // Index of the closest point.
            let m = Float::with_val(prec, qt / r).round();
            let d = Float::with_val(prec, &m * r - qt);
            let d = d.square();
            return in_radius(&d, rad)
                .then(|| IVector::from_array([m.to_integer().unwrap()]));
        }

        let qtc = &qt[i];
        let rc = &r[(i, i)];

        // Index of the closest plane before rounding.
        let start_index_fl = Float::with_val(prec, qtc / rc);
        let start_index = Float::with_val(prec, start_index_fl.round_ref());
        //println!("Start index: {start_index} ({start_index_fl})");

        // Suppose the start index was -0.4.
        // We would want to check the planes in the order
        // 0, -1, 1, -2, 2.
        // But if the start index was 0.4, we would want
        // 0, 1, -1, 2, -2
        // So depending on whether we round up or down
        // to the integer start index, we will negate
        // the offset from the start index.
        let negate_offset = (start_index_fl - &start_index).is_sign_negative();

        // The current plane's offset.
        let mut offset = 0isize;

        let mut min_dist = rad.clone();
        let mut min = None;

        // Iterate over all possible planes.
        loop {
            // Index of the plane we are considering.
            let index = Float::with_val(prec, &start_index + offset);

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
            let d = index.clone() * rc - qtc;
            let d = d.square();
            //println!("Index: {index} Distance: {d}");

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
            let plane_dist = min_dist.as_ref().map(|f| f - d);

            // Compute the point in the plane we need to be close to now.
            let point_in_plane = qt - &index * &r[i];

            // Recursively find the closest point.
            let v = cvp_impl(i - 1, r, &point_in_plane, prec, &plane_dist);
            let Some(mut v) = v else {
                continue
            };
            assert!(v.dim == i);

            // v is the new index vector of the point.
            // It is the index of the closest point of the previous call,
            // plus the index of the plane.
            v.append(index.to_integer().unwrap());

            // Compute the actual vector of the index vector,
            // i.e. v * r[:i+1, :i+1]. Since there is currently no
            // way to multiply by a submatrix, just do it manually.
            let iter = (0..v.dim).map(|i| v.iter().zip(r.column(i))
                .map(|(l, r)| Float::with_val(prec, l) * r)
                .fold(Float::with_val(prec, 0), |acc, f| acc + f)
            );
            let vv = FVector::from_iter(v.dim, iter);

            // Compute the distance to the point.
            let d = (&vv - &qt[0..i+1]).norm_sqr();

            // If the distance is smaller than the current minimal dist,
            // then we have found the new best point.
            if in_radius(&d, &min_dist) {
                min = Some(v);
                min_dist = Some(d);
            }
        }

        min
    }

    // Multiply the coefficients by the basis.
    cvp_impl(r.cols - 1, &r, &qt, prec, &rad)
        .map(|v| v.view() * basis)
}

#[test]
fn cvp_temp_test() {
    let b = IMatrix::from_rows(&[
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
    ]);
    let t = IVector::from_entries([-48, 69, -76, 36, -72, 31, -53, -7, 54,
        74, 6, -82, -13, -32, 7, 53, -60, -44, 38, -97]);
    //let now = std::time::Instant::now();
    //for _ in 0..100 {
    //    let v = cvp_planes(&b, &t, 53, None);
    //}
    //println!("{}", now.elapsed().as_secs_f64());

    fn bench(name: &str, f: impl Fn() -> IVector, t: &IVector) {
        let now = std::time::Instant::now();
        let mut v = IVector::empty();
        for _ in 0..100 {
            v = f();
        }
        let time = now.elapsed().as_secs_f64() * 1000.;
        let d = (&v - t).norm();
        println!("{name}: {d:.3} in {time:.3}ms");
    }

    let now = std::time::Instant::now();
    bench("exact solution", || cvp_planes(&b, &t, 53, None).unwrap(), &t);
    bench("nearest plane", || cvp_nearest_plane_float(&b, &t, 53), &t);
    bench("rounding", || cvp_rounding_float(&b, &t), &t);
}

#[test]
fn cvp_exact_dim20() {
    let b = IMatrix::from_rows(&[
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
    ]);
    let t = IVector::from_entries([-48, 69, -76, 36, -72, 31, -53, -7, 54,
        74, 6, -82, -13, -32, 7, 53, -60, -44, 38, -97]);
    assert!(cvp_planes(&b, &t, 53, None).unwrap().view() == &[-30, 35, -98,
        61, -27, 75, -32, -3, 70, 8, 3, -77, -29, -103, 61, 58, -71, 41, 37, -40]);
}

/// Gram-Schmidt of the rows of the matrix `a`.
/// If `orthonormal` is true, the result will be normalized.
fn gram_schmidt(mut a: FMatrix, orthonormal: bool) -> FMatrix {
    for i in 0..a.rows {
        for j in 0..i {
            let f = a[j].dot(&a[i]) / a[j].norm_sqr();
            let p = &f * &a[j];
            a[i] -= &p;
        }
    }

    if orthonormal {
        for i in 0..a.rows {
            let n = a[i].norm();
            a[i] /= &n;
        }
    }

    a
}

/// Computes the RQ-decomposition (row QR-decomposition) of a matrix.
/// Returns a lower triangular matrix `R` and an orthogonal matrix `Q`,
/// such that `R * Q = A`.
pub fn rq_decomposition(a: &FMatrix) -> (FMatrix, FMatrix) {
    let q = gram_schmidt(a.clone(), true);
    let r = a.mul_transpose(&q);
    (r, q)
}

#[test]
fn gram_schmidt_test() {
    let a = Matrix::<i32>::from_rows(&[
        [1, 2, 3],
        [3, 4, 5],
    ]);
    let a = a.to_float(53);
    let (r, q) = rq_decomposition(&a);
    println!("q: {:?}\nr: {:?}", q, r);
    println!("{:?}", &r * &q);
    println!("{:?}", &q * &FVector::from_array(
        [Float::with_val(53, 5), Float::with_val(53, 6), Float::with_val(53, 7)]));
}

#[test]
fn nearest_plane_example() {
    let a = IMatrix::from_rows(&[
        [2, 3, 1],
        [4, 1, -3],
        [2, 2, 2],
    ]);

    let t = IVector::from_entries([4, 2, 7]);
    let c = cvp_nearest_plane_float(&a, &t, 53);
    println!("{:?}", c);
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
    fn solve_linear(a: Matrix<Self>, b: Vector<Self>) -> Option<Vector<Self>>;
}

impl Field for Rational {
    fn from_int(i: &Integer) -> Self {
        Rational::from(i)
    }

    fn to_int(&self) -> Integer {
        self.round_ref().into()
    }

    fn solve_linear(a: Matrix<Self>, b: Vector<Self>) -> Option<Vector<Self>> {
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

    fn solve_linear(a: Matrix<Self>, b: Vector<Self>) -> Option<Vector<Self>> {
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

    fn solve_linear(a: Matrix<Self>, b: Vector<Self>) -> Option<Vector<Self>> {
        solve_linear(a, b)
    }
}


fn solve_linear_rat(
    mut a: Matrix<Rational>, mut b: Vector<Rational>
) -> Option<Vector<Rational>> {
    assert!(a.rows == a.cols,
        "This function only supports non-singular square systems.");
    for i in 0..a.cols {
        // Choose a pivot in the c-th column.
        let pivot = a.column(i)
            .enumerate()
            .skip(i)
            .filter(|e| e.1 != &Rational::new())
            .max_by(|e, f| e.1.cmp_abs(f.1))?.0;
        a.swap_rows(pivot, i);
        let pivot = a[(i, i)].clone();
        for r in i+1..a.rows {
            let fac = (&a[(r, i)] / &pivot).complete();
            for c in i+1..a.cols {
                let s = (&fac * &a[(i, c)]).complete();
                a[(r, c)] -= s;
            }
            let s = fac * &b[i];
            b[r] -= s;
        }
    }

    let mut result = Vector::<Rational>::zero(a.cols);
    for i in (0..a.cols).rev() {
        let mut sum = b[i].clone();
        for j in i+1..a.cols {
            sum -= (&a[(i, j)] * &result[j]).complete();
        }

        sum /= &a[(i, i)];
        result[i] = sum;
    }

    Some(result)
}

/// PartialOrd should not return None for any of the elements in the matrix.
/// We can't use Ord because of the floating point types.
fn solve_linear<T: NumAssign + PartialOrd + Copy>(
    mut a: Matrix<T>, mut b: Vector<T>
) -> Option<Vector<T>> {
    assert!(a.rows == a.cols,
        "This function only supports non-singular square systems.");
    let abs = |x: T| if x < T::zero() { T::zero() - x } else { x };
    for i in 0..a.cols {
        // Choose a pivot in the c-th column.
        let pivot = a.column(i)
            .enumerate()
            .skip(i)
            .filter(|e| *e.1 != T::zero())
            .max_by(|e, f| abs(*e.1).partial_cmp(&abs(*f.1)).unwrap())?.0;
        a.swap_rows(pivot, i);
        let pivot = a[(i, i)];
        for r in i+1..a.rows {
            let fac = a[(r, i)] / pivot;
            for c in i+1..a.cols {
                let s = fac * a[(i, c)];
                a[(r, c)] -= s;
            }
            let s = b[i] * fac;
            b[r] -= s;
        }
    }

    let mut result = Vector::zero(a.cols);
    for i in (0..a.cols).rev() {
        let mut sum = b[i];
        for j in i+1..a.cols {
            sum -= a[(i, j)] * result[j];
        }

        sum /= a[(i, i)];
        result[i] = sum;
    }

    Some(result)
}
