//! Binary permutation polynomials.

use num_traits::{Zero, One, Euclid};
use rand::distributions::{Uniform, Distribution};
use num_bigint::{BigInt, RandBigInt};
use num_integer::Integer;
use crate::{poly::Poly, keep_bits_mut, matrix::IOwnedMatrix, vector::IOwnedVector, diophantine};

/// Returns a pair of permutation polynomials mod 2^n.
/// The functions are inverses of each other.
/// `degree` specifies the degree of the first polynomial.
pub fn perm_pair(zi: &ZeroIdeal, degree: usize) -> (Poly, Poly) {
    assert!(degree >= 1, "Can't create a permutation polynomial \
            of degree 0 as this would be a constant.");
    // Generate a random permutation polynomial.
    let p = random_perm_poly(zi, degree);
    assert!(is_perm_poly(&p), "Generated an invalid polynomial.");

    // Find the inverse.
    let q = compute_inverse(&p, zi);

    (p, q)
}

/// Returns a random permutation polynomial.
pub fn random_perm_poly(zi: &ZeroIdeal, degree: usize) -> Poly {
    let rng = &mut rand::thread_rng();
    let mut p: Vec<_> =  (0..=degree)
        .map(|_| rng.gen_bigint(zi.bits as u64))
        .collect();

    // Make sure that this is a permutation polynomial.
    // The coefficient of degree 1 has to be odd.
    if p[1].is_even() {
        p[1] += 1;
        keep_bits_mut(&mut p[1], zi.bits);
    }

    let mut rng = rand::thread_rng();

    // Make sure the sum of the even coefficients (except 0) is even.
    if p.iter().skip(2).step_by(2).fold(false, parity) {
        let dist = Uniform::from(1..=degree/2);
        let i = dist.sample(&mut rng);
        p[2*i] += 1;
        keep_bits_mut(&mut p[2*i], zi.bits);
    }

    // Make sure the sum of the odd coefficients (except 1) is even.
    if p.iter().skip(3).step_by(2).fold(false, parity) {
        let dist = Uniform::from(1..=(degree-1)/2);
        let i = dist.sample(&mut rng);
        p[2*i+1] += 1;
        keep_bits_mut(&mut p[2*i+1], zi.bits);
    }

    Poly { coeffs: p }.truncated()
}

/// Computes the composition p(q(x)) of two polynomials.
/// It does not simplify the result as much as possible,
/// so it is recommended to call [Poly::simplify] afterwards.
/// It does however reduce the degree of the result.
pub fn compose(p: &Poly, q: &Poly, zi: &ZeroIdeal) -> Poly {
    // Iterate over the coefficients in reverse order.
    let mut iter = p.coeffs.iter().rev();

    let mut r = Poly::constant(
        iter.next().map_or_else(Zero::zero, |c| c.clone())
    );

    // The last coefficient is the initial value.
    for c in iter {
        // Multiply the current value by a and add the next coefficient.
        r *= q;
        r += c;
        r.reduce(zi);
    }

    r
}

/// Computes the inverse of a permutation polynomial using Newton's Method.
pub fn compute_inverse(f: &Poly, zi: &ZeroIdeal) -> Poly {
    assert!(is_perm_poly(f), "Can't invert {f} as it is not a permutation.");

    // Simplify p.
    let p = f.clone().simplified(zi);

    // Initialize q with the initial guess: p(x) = x.
    let mut q = Poly::from_vec(vec![0, 1]);

    let mut it = 0;

    // Do the newton method.
    loop {
        // Not proven to always work so make sure
        // to stop after a certain number of iterations.
        assert!(it <= zi.bits * 2, "Failed to compute the inverse \
                in a reasonable number of iterations.");
        // Compute the composition.
        let mut comp = compose(&p, &q, zi)
            .simplified(zi);

        // Do we already have p(q(x)) = x?
        if comp.is_id() {
            return q;
        }

        // Subtract x.
        // This is the quantity we want to make 0.
        comp.coeffs[1] -= 1;

        // Update the guess.
        let qd = q.derivative();
        q -= &(&qd * &comp);
        q.simplify(zi);

        it += 1;
    }
}

/// Computes the inverse of a permutation polynomial by using f as a generator.
pub fn compute_inverse_generator(f: &Poly, zi: &ZeroIdeal) -> Poly {
    assert!(is_perm_poly(f), "Can't invert {f} as it is not a permutation.");

    // The inverse will contain f^(2^n-1)
    let mut inverse = f.clone();

    for _ in 0..zi.bits {
        // Compose inverse and f.
        let g = compose(&inverse, f, zi).simplified(zi);

        // If it is the identity, then we found the inverse.
        if g.is_id() {
            return inverse;
        }

        // Compose f^(2^n-1) and f^(2^n) to
        // get f^(2^(n+1)-1) for the next iteration.
        inverse = compose(&inverse, &g, zi).simplified(zi);
    }

    panic!("Failed to compute the inverse if {} mod 2^{}", f, zi.bits);
}

/// Compute the inverse using interpolation.
pub fn compute_inverse_interpolation(f: &Poly, zi: &ZeroIdeal) -> Poly {
    // Construct a system of linear congruences.
    let rows = zi.gen.last().unwrap().len();
    let cols = zi.gen.last().unwrap().len();

    // Construct the Vandermonde matrix.
    let mut a = IOwnedMatrix::zero(rows, cols);
    let mut i = BigInt::zero();
    for r in 0..rows {
        let mut j = BigInt::one();
        let x = f.eval_bits(&i, zi.bits);
        for c in 0..cols {
            a[(r, c)] = j.clone();
            j *= &x;
        }

        i += 1;
    }

    // Construct the vector of values of the polynomial.
    let mut b = IOwnedVector::zero(rows);
    for r in 0..rows {
        b[r] = BigInt::from(r);
    }

    let l = diophantine::solve_congruences(a, b, zi.bits);

    for b in l.lattice.basis.rows() {
        let k = Poly::from_vec(b.iter().cloned().collect());
        assert!(k.clone().simplified(zi).is_zero(),
            "Polynomial in kernel is not null: {k}");
    }

    Poly::from_vec(l.offset.iter().cloned().collect()).simplified(zi)
}

/// Computes a set of generators for the "zero ideal" of Z_{2^n}[x],
/// that is all polynomials in Z_{2^n}[x] s.t. p(x) = 0 for all Z_{2^n}.
/// One polynomial p_0(x) = 0 = 2^n is ignored,
/// but it technically should be here.
/// This function could be a lot more efficient:
/// We only reduce the coefficients mod 2^n after the whole polynomial
/// is expanded. Instead we could reduce during the expansion,
/// specifically in mul_linfac.
/// Or calculate with n bit integers the whole time, instead of rug::Integers.
fn zero_ideal(n: u32) -> Vec<Poly> {
    let mut gen = Vec::new();

    // div stores how often 2 divides i!.
    // It is successively updated.
    let mut div = 0;
    for i in (2u32..).step_by(2) {
        div += i.trailing_zeros();

        // Compute the exponent.
        let e = if n <= div { 0 } else { n - div };

        // Let's build the polynomial.
        let mut p = Poly::one();

        for j in 0..i {
            // Multiply the current polynomial by (x-j).
            p.mul_linfac(&j.into());
        }

        p.shl_coeff(e);
        p.mod_coeff(n);
        p.truncate();

        gen.push(p);

        if e == 0 {
            break;
        }
    }

    gen
}

/// This is the zero ideal usually specified in papers,
/// but it is not minimal.
#[allow(dead_code)]
fn zero_ideal_redundant(n: u32) -> Vec<Poly> {
    let mut gen = Vec::new();

    // div stores how often 2 divides i!.
    // It is successively updated.
    let mut div = 0;
    for i in 2u32.. {
        div += i.trailing_zeros();

        // If the exponent would be negative
        // then add the last generator and stop.
        if n <= div {
            let mut p = Poly::one();

            for j in 0..i {
                // Multiply the current polynomial by (x-j).
                p.mul_linfac(&j.into());
            }

            p.mod_coeff(n);
            p.truncate();

            gen.push(p);
            break;
        }

        // Compute the exponent.
        let e = n - div;

        // Let's build the polynomial.
        let mut p = Poly::one();

        for j in 0..i {
            // Multiply the current polynomial by (x-j).
            p.mul_linfac(&j.into());
        }

        p.shl_coeff(e);
        p.mod_coeff(n);
        p.truncate();

        gen.push(p);
    }

    gen
}

/// Used internally as a fold function.
/// Computes the "parity" of a list of integers,
/// that is whether their sum is even or odd,
/// depending on the initial value of the accumulator.
fn parity(acc: bool, i: &BigInt) -> bool {
    i.is_odd() ^ acc
}

/// Is this function a permutation polynomial?
pub fn is_perm_poly(f: &Poly) -> bool {
    return f.coeffs.get(1).map_or(false, |i| i.is_odd())
        && f.coeffs.iter().skip(2).step_by(2).fold(true, parity)
        && f.coeffs.iter().skip(3).step_by(2).fold(true, parity);
}

/// The zero ideal is the ideal of all polynomials
/// that evaluate to zero everywhere in the polynomial ring Z_{2^n}\[X\].
pub struct ZeroIdeal {
    /// The coefficients are mod 2^n.
    bits: u32,

    /// The generators.
    gen: Vec<Poly>,
}

impl ZeroIdeal {
    /// Initializes the zero ideal.
    pub fn init(bits: u32) -> Self {
        assert!(bits > 0, "Not a valid ring.");
        Self { bits, gen: zero_ideal(bits) }
    }

    /// The number of bits.
    pub fn bits(&self) -> u32 {
        self.bits
    }

    /// Get generators for the zero ideal.
    pub fn generators(&self) -> &[Poly] {
        &self.gen
    }
}

impl Poly {
    pub fn simplified(mut self, zi: &ZeroIdeal) -> Self {
        self.simplify(zi);
        self
    }

    /// Simplifies a polynomial by adding a polynomial in the zero ideal
    /// to reduce the degree of the polynomial.
    pub fn simplify(&mut self, zi: &ZeroIdeal) {
        // We reduce all polynomials with degree larger or equal to that of
        // zi.last with zi.last, because it has a leading coefficient of 1,
        // so each coefficient of higher degree
        // can be eliminated in a single step.
        // Once the degree is less than that of zi.last we try to reduce with the
        // other generators as much as possible.

        if self.len() == 0 {
            return;
        }

        let mut coeff = self.len() - 1;

        for gen in zi.gen.iter().rev() {
            let gen_len = gen.len();

            while coeff + 1 >= gen_len {
                let m = self.coeffs[coeff].div_euclid(&gen.coeffs[gen_len-1]);
                if !m.is_zero() {
                    let iter = self.coeffs[coeff+1-gen_len..=coeff]
                        .iter_mut().zip(gen.coeffs.iter());

                    for (p, g) in iter {
                        *p -= &m * g;
                    }
                }
                coeff -= 1;
            }
        }

        // Reduce all coefficients mod 2^n.
        self.mod_coeff(zi.bits);

        // Remove leading coefficients that are 0.
        self.truncate();
    }

    /// Reduces the degree of a polynomial if possible.
    /// This function should be used during computations
    /// where the degree could otherwise explode.
    /// It is a simplified version of the [Poly::simplify] algorithm,
    /// that only uses the generator of the highest degree.
    pub fn reduce(&mut self, zi: &ZeroIdeal) {
        if self.len() == 0 {
            return;
        }

        let gen = zi.gen.last().unwrap();
        let gen_len = gen.len();
        while self.len() >= gen_len {
            let p_len = self.len();
            let (c, rest) = self.coeffs.split_last_mut().unwrap();
            if !c.is_zero() {
                let iter = rest[p_len-gen_len..].iter_mut()
                    .zip(gen.coeffs.iter());

                for (p, g) in iter {
                    *p -= &*c * g;
                }
            }
            self.coeffs.pop();
        }

        // Reduce all coefficients mod 2^n.
        self.mod_coeff(zi.bits);

        // Remove leading coefficients that are 0.
        self.truncate();
    }
}

#[test]
fn check_inverse_8() {
    let zi = ZeroIdeal::init(8);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..10 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&zi, d);
        let comp = compose(&p, &q, &zi).simplified(&zi);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
fn check_inverse_16() {
    let zi = ZeroIdeal::init(8);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..10 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&zi, d);
        let comp = compose(&p, &q, &zi).simplified(&zi);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
fn check_inverse_32() {
    let zi = ZeroIdeal::init(32);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..10 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&zi, d);
        let comp = compose(&p, &q, &zi).simplified(&zi);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
pub fn check_inverse_64() {
    let zi = ZeroIdeal::init(64);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..10 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&zi, d);
        let comp = compose(&p, &q, &zi).simplified(&zi);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
pub fn zero_ideal_test() {
    let zi = ZeroIdeal::init(4);
    for (i, g) in zi.gen.iter().enumerate() {
        println!("{}: {}", i, g);
        for x in 0..16 {
            let y = g.eval_bits(&x.into(), 4);
            if !y.is_zero() {
                panic!("f({}) = {}", x, y);
            }
        }
    }
}

/// Computes the order of a polynomial by brute force.
/// Should only be used for testing and only with tiny n.
pub fn order(p: &Poly, zi: &ZeroIdeal) -> usize {
    assert!(is_perm_poly(p), "The order is only defined for permutation polynomials!");
    assert!(zi.bits <= 8, "'order' should only be used with tiny n!");

    if p.is_id() {
        return 1;
    }

    let mut f = p.clone();
    for i in 1..1usize << zi.bits {
        f = compose(&f, p, zi).simplified(zi);
        if f.is_id() {
            return i + 1;
        }
    }

    panic!("We shouldn't get here. \
        Either p is not a permutation polynomial or the composition is wrong.");
}

#[test]
pub fn test_order() {
    let zi = ZeroIdeal::init(8);
    for d in 0..10usize {
        let p = random_perm_poly(&zi, d % 10 + 1);
        let o = order(&p, &zi);
        assert!(o.is_power_of_two(), "Order of {} is {}!", p, o);
    }
}

//#[test]
//pub fn test_inverse_generator() {
//    let zi = ZeroIdeal::init(64);
//    for d in 0..10usize {
//        let p = random_perm_poly(&zi, d % 10 + 1).simplified(&zi);
//        let q = compute_inverse_generator(&p, &zi);
//        //let q = compute_inverse(&p, &zi);
//        assert!(compose(&p, &q, &zi).simplified(&zi).is_id(),
//            "compute_inverse_generator returned wrong inverse {} of {}", q, p);
//    }
//}