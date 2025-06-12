//! Binary permutation polynomials.

use rand::distr::{Uniform, Distribution};
use rand::Rng;
use crate::solver;
use crate::matrix::OwnedMatrix;
use crate::rings::{Ring, RingElement as _, BinaryRing};
use crate::poly::Poly;
use crate::vector::OwnedVector;

/// Returns a pair of permutation polynomials mod 2^n.
/// The functions are inverses of each other.
/// `degree` specifies the degree of the first polynomial.
pub fn perm_pair<R: BinaryRing, Rand: Rng>(
    rng: &mut Rand,
    zi: &ZeroIdeal<R>,
    degree: usize,
    r: &R,
) -> (Poly<R>, Poly<R>) {
    assert!(degree >= 1, "Can't create a permutation polynomial \
            of degree 0 as this would be a constant.");
    // Generate a random permutation polynomial.
    let p = random_perm_poly(rng, degree, r);
    assert!(is_perm_poly(&p), "Generated an invalid polynomial.");

    // Find the inverse.
    let q = compute_inverse(&p, zi, r);

    (p, q)
}

/// Returns a random permutation polynomial.
pub fn random_perm_poly<R: BinaryRing, Rand: Rng>(
    rng: &mut Rand,
    degree: usize,
    r: &R,
) -> Poly<R> {
    let mut p: Vec<_> =  (0..=degree)
        .map(|_| r.random(rng))
        .collect();

    // Make sure that this is a permutation polynomial.
    // The coefficient of degree 1 has to be odd.
    if R::is_even(&p[1]) {
        r.inc_assign(&mut p[1]);
    }

    // Make sure the sum of the even coefficients (except 0) is even.
    if p.iter().skip(2).step_by(2).fold(false, parity::<R>) {
        let dist = Uniform::new_inclusive(1, degree / 2).unwrap();
        let i = dist.sample(rng);
        r.inc_assign(&mut p[2*i]);
    }

    // Make sure the sum of the odd coefficients (except 1) is even.
    if p.iter().skip(3).step_by(2).fold(false, parity::<R>) {
        let dist = Uniform::new_inclusive(1, (degree - 1) / 2).unwrap();
        let i = dist.sample(rng);
        r.inc_assign(&mut p[2*i+1]);
    }

    Poly { coeffs: p }.truncated()
}

/// Computes the composition p(q(x)) of two polynomials.
/// It does not simplify the result as much as possible,
/// so it is recommended to call [Poly::simplify] afterwards.
/// It does however reduce the degree of the result.
pub fn compose<R: BinaryRing>(
    p: &Poly<R>,
    q: &Poly<R>,
    zi: &ZeroIdeal<R>,
    ring: &R,
) -> Poly<R> {
    // Iterate over the coefficients in reverse order.
    let mut iter = p.coeffs.iter().rev();

    let mut r = Poly::constant(
        iter.next().map_or_else(R::zero, |c| c.clone())
    );

    // The last coefficient is the initial value.
    for c in iter {
        // Multiply the current value by a and add the next coefficient.
        r.mul_assign(q, ring);
        r.add_assign_const(c, ring);
        r.reduce(zi, ring);
    }

    r
}

/// Computes the inverse of a permutation polynomial using Newton's Method.
pub fn compute_inverse<R: BinaryRing>(
    f: &Poly<R>,
    zi: &ZeroIdeal<R>,
    r: &R,
) -> Poly<R> {
    assert!(is_perm_poly(f), "Can't invert {f} as it is not a permutation.");

    // Simplify p.
    let p = f.clone().simplified(zi, r);

    // Initialize q with the initial guess: p(x) = x.
    let mut q = Poly::from_vec(vec![R::zero(), R::one()]);

    let mut it = 0;

    // Do the newton method.
    loop {
        // Not proven to always work so make sure
        // to stop after a certain number of iterations.
        assert!(it <= r.bits() * 2, "Failed to compute the inverse \
                in a reasonable number of iterations.");
        // Compute the composition.
        let mut comp = compose(&p, &q, zi, r)
            .simplified(zi, r);

        // Do we already have p(q(x)) = x?
        if comp.is_id() {
            return q;
        }

        // Subtract x.
        // This is the quantity we want to make 0.
        r.dec_assign(&mut comp.coeffs[1]);

        // Update the guess.
        let qd = q.derivative(r);
        q.sub_assign(&(qd.mul(&comp, r)), r);
        q.simplify(zi, r);

        it += 1;
    }
}

/// Computes the inverse of a permutation polynomial by using f as a generator.
pub fn compute_inverse_generator<R: BinaryRing>(
    f: &Poly<R>,
    zi: &ZeroIdeal<R>,
    r: &R,
) -> Poly<R> {
    assert!(is_perm_poly(f), "Can't invert {f} as it is not a permutation.");

    // The inverse will contain f^(2^n-1)
    let mut inverse = f.clone();

    for _ in 0..r.bits() {
        // Compose inverse and f.
        let g = compose(&inverse, f, zi, r).simplified(zi, r);

        // If it is the identity, then we found the inverse.
        if g.is_id() {
            return inverse;
        }

        // Compose f^(2^n-1) and f^(2^n) to
        // get f^(2^(n+1)-1) for the next iteration.
        inverse = compose(&inverse, &g, zi, r).simplified(zi, r);
    }

    panic!("Failed to compute the inverse if {} mod 2^{}", f, r.bits());
}

/// Compute the inverse using interpolation.
pub fn compute_inverse_interpolation<R: BinaryRing>(
    f: &Poly<R>,
    zi: &ZeroIdeal<R>,
    ring: &R,
) -> Poly<R> {
    // Construct a system of linear congruences.
    let rows = zi.generators.last().unwrap().len();
    let cols = zi.generators.last().unwrap().len();

    // Construct the Vandermonde matrix.
    let mut a = OwnedMatrix::<R>::zero(rows, cols);
    let mut i = R::zero();
    for r in 0..rows {
        let mut j = R::one();
        let x = f.eval(&i, ring);
        for c in 0..cols {
            a[(r, c)] = j.clone();
            ring.mul_assign(&mut j, &x);
        }

        ring.inc_assign(&mut i);
    }

    // Construct the vector of values of the polynomial.
    let mut b = OwnedVector::zero(rows);
    for r in 0..rows {
        b[r] = ring.element_from_usize(r);
    }

    let l = solver::solve_via_integer_diagonalize(a, b, ring);

    for b in l.lattice.basis.rows() {
        let k = Poly::from_vec(b.iter().cloned().collect());
        assert!(k.clone().simplified(zi, ring).is_zero(),
            "Polynomial in kernel is not null: {k}");
    }

    Poly::from_vec(l.offset.iter().cloned().collect()).simplified(zi, ring)
}

/// Computes a set of generators for the "zero ideal" of Z_{2^n}[x],
/// that is all polynomials in Z_{2^n}[x] s.t. p(x) = 0 for all Z_{2^n}.
/// One polynomial p_0(x) = 0 = 2^n is ignored,
/// but it technically should be here.
fn zero_ideal<R: BinaryRing>(r: &R) -> Vec<Poly<R>> {
    let mut generators = Vec::new();

    // div stores how often 2 divides i!.
    // It is successively updated.
    let mut div = 0;
    for i in (2u32..).step_by(2) {
        div += i.trailing_zeros();

        // Compute the exponent.
        let e = r.bits().saturating_sub(div);

        // Let's build the polynomial.
        let mut p = Poly::one();

        for j in 0..i {
            // Multiply the current polynomial by (x-j).
            p.mul_linfac(&r.element_from_usize(j as usize), r);
        }

        p.shl_coeff(e, r);
        p.truncate();

        generators.push(p);

        if e == 0 {
            break;
        }
    }

    generators
}

/// This is the zero ideal usually specified in papers,
/// but it is not minimal.
#[allow(dead_code)]
fn zero_ideal_redundant<R: BinaryRing>(r: &R) -> Vec<Poly<R>> {
    let mut generators = Vec::new();

    // div stores how often 2 divides i!.
    // It is successively updated.
    let mut div = 0;
    for i in 2u32.. {
        div += i.trailing_zeros();

        // If the exponent would be negative
        // then add the last generator and stop.
        if r.bits() <= div {
            let mut p = Poly::one();

            for j in 0..i {
                // Multiply the current polynomial by (x-j).
                p.mul_linfac(&r.element_from_usize(j as usize), r);
            }

            p.truncate();

            generators.push(p);
            break;
        }

        // Compute the exponent.
        let e = r.bits() - div;

        // Let's build the polynomial.
        let mut p = Poly::one();

        for j in 0..i {
            // Multiply the current polynomial by (x-j).
            p.mul_linfac(&r.element_from_usize(j as usize), r);
        }

        p.shl_coeff(e, r);
        p.truncate();

        generators.push(p);
    }

    generators
}

/// Used internally as a fold function.
/// Computes the "parity" of a list of integers,
/// that is whether their sum is even or odd,
/// depending on the initial value of the accumulator.
fn parity<R: BinaryRing>(acc: bool, i: &R::Element) -> bool {
    R::is_odd(i) ^ acc
}

/// Is this function a permutation polynomial?
pub fn is_perm_poly<R: BinaryRing>(f: &Poly<R>) -> bool {
    f.coeffs.get(1).is_some_and(|i| R::is_odd(i))
        && f.coeffs.iter().skip(2).step_by(2).fold(true, parity::<R>)
        && f.coeffs.iter().skip(3).step_by(2).fold(true, parity::<R>)
}

/// The zero ideal is the ideal of all polynomials
/// that evaluate to zero everywhere in the polynomial ring Z_{2^n}\[X\].
pub struct ZeroIdeal<R: Ring> {
    /// The generators.
    generators: Vec<Poly<R>>,
}

impl<R: Ring> ZeroIdeal<R> {
    /// Initializes the zero ideal.
    pub fn init(r: &R) -> Self
    where
        R: BinaryRing
    {
        Self { generators: zero_ideal(r) }
    }

    /// Get generators for the zero ideal.
    pub fn generators(&self) -> &[Poly<R>] {
        &self.generators
    }
}

impl<R: BinaryRing> Poly<R> {
    pub fn simplified(mut self, zi: &ZeroIdeal<R>, r: &R) -> Self {
        self.simplify(zi, r);
        self
    }

    /// Simplifies a polynomial by adding a polynomial in the zero ideal
    /// to reduce the degree of the polynomial.
    pub fn simplify(&mut self, zi: &ZeroIdeal<R>, r: &R) {
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

        for g in zi.generators.iter().rev() {
            let gen_len = g.len();

            while coeff + 1 >= gen_len {
                let m = R::euclidean_div(
                    &self.coeffs[coeff],
                    &g.coeffs[gen_len-1]
                );
                if !m.is_zero() {
                    let iter = self.coeffs[coeff+1-gen_len..=coeff]
                        .iter_mut().zip(g.coeffs.iter());

                    for (p, g) in iter {
                        r.mul_sub_assign(p, &m, g);
                    }
                }
                coeff -= 1;
            }
        }

        // Remove leading coefficients that are 0.
        self.truncate();
    }

    /// Reduces the degree of a polynomial if possible.
    /// This function should be used during computations
    /// where the degree could otherwise explode.
    /// It is a simplified version of the [Poly::simplify] algorithm,
    /// that only uses the generator of the highest degree.
    pub fn reduce(&mut self, zi: &ZeroIdeal<R>, r: &R) {
        if self.len() == 0 {
            return;
        }

        let g = zi.generators.last().unwrap();
        let gen_len = g.len();
        while self.len() >= gen_len {
            let p_len = self.len();
            let (c, rest) = self.coeffs.split_last_mut().unwrap();
            if !c.is_zero() {
                let iter = rest[p_len-gen_len..].iter_mut()
                    .zip(g.coeffs.iter());

                for (p, g) in iter {
                    r.mul_sub_assign(p, c, g);
                }
            }
            self.coeffs.pop();
        }

        // Remove leading coefficients that are 0.
        self.truncate();
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::rings::{U8, U16, U32, U64, U128};
    use rand::{rngs::StdRng, SeedableRng as _};

    fn check_composition<R: BinaryRing>(r: &R) {
        let rng = &mut StdRng::seed_from_u64(0);
        let zi = ZeroIdeal::init(r);
        let deg = Uniform::new(1, 16).unwrap();
        for _ in 0..100 {
            let deg_p = deg.sample(rng);
            let deg_q = deg.sample(rng);
            let p = Poly::random(rng, deg_p, r);
            let q = Poly::random(rng, deg_q, r);

            let comp = compose(&p, &q, &zi, r);

            for _ in 0..5 {
                let x = r.random(rng);
                assert_eq!(
                    comp.eval(&x, r),
                    p.eval(&q.eval(&x, r), r)
                );
            }
        }
    }

    #[test]
    fn check_composition_8() {
        check_composition(&U8);
    }

    #[test]
    fn check_composition_16() {
        check_composition(&U16);
    }

    #[test]
    fn check_composition_32() {
        check_composition(&U32);
    }

    #[test]
    fn check_composition_64() {
        check_composition(&U64);
    }

    #[test]
    fn check_composition_128() {
        check_composition(&U128);
    }

    fn check_inverse<R: BinaryRing>(r: &R) {
        let mut rng = StdRng::seed_from_u64(0);
        let zi = ZeroIdeal::init(r);
        let deg = Uniform::new(1, 16).unwrap();
        for _ in 0..10 {
            let d = deg.sample(&mut rng);
            let (p, q) = perm_pair(&mut rng, &zi, d, r);
            let comp = compose(&p, &q, &zi, r).simplified(&zi, r);
            assert!(comp.is_id(), "Incorrect inverse.");
        }
    }

    #[test]
    fn check_inverse_8() {
        check_inverse(&U8);
    }

    #[test]
    fn check_inverse_16() {
        check_inverse(&U16);
    }

    #[test]
    fn check_inverse_32() {
        check_inverse(&U32);
    }

    #[test]
    pub fn check_inverse_64() {
        check_inverse(&U64);
    }

    #[test]
    pub fn check_inverse_128() {
        check_inverse(&U128);
    }

    #[test]
    pub fn zero_ideal_test() {
        //let r = UniformBigInt::new(4);
        let r = &U8;
        let zi = ZeroIdeal::init(r);
        for (i, g) in zi.generators.iter().enumerate() {
            println!("{i}: {g}");
            for x in 0usize..16 {
                let y = g.eval(&r.element_from_usize(x), r);
                if !y.is_zero() {
                    panic!("f({x}) = {y}");
                }
            }
        }
    }

    /// Computes the order of a polynomial by brute force.
    /// Should only be used for testing and only with tiny moduli.
    pub fn order<R: BinaryRing>(
        p: &Poly<R>,
        zi: &ZeroIdeal<R>,
        r: &R,
    ) -> usize {
        assert!(is_perm_poly(p), "The order is only defined for permutation polynomials!");
        assert!(r.bits() <= 8, "'order' should only be used with tiny n!");

        if p.is_id() {
            return 1;
        }

        let mut f = p.clone();
        for i in 1..1usize << r.bits() {
            f = compose(&f, p, zi, r).simplified(zi, r);
            if f.is_id() {
                return i + 1;
            }
        }

        panic!("We shouldn't get here. \
            Either p is not a permutation polynomial or the composition is wrong.");
    }

    #[test]
    pub fn test_order() {
        let r = &U8;
        let mut rng = StdRng::seed_from_u64(0);
        let zi = ZeroIdeal::init(r);
        for d in 0..10usize {
            let p = random_perm_poly(&mut rng, d % 10 + 1, r);
            let o = order(&p, &zi, r);
            assert!(o.is_power_of_two(), "Order of {p} is {o}!");
        }
    }

    #[test]
    #[ignore = "very slow and `compute_inverse_generator` shouldn't be used"]
    pub fn test_inverse_generator() {
        let r = &U64;
        let mut rng = StdRng::seed_from_u64(0);
        let zi = ZeroIdeal::init(r);
        for d in 0..10usize {
            let p = random_perm_poly(&mut rng, d % 10 + 1, r).simplified(&zi, r);
            let q = compute_inverse_generator(&p, &zi, r);
            assert!(compose(&p, &q, &zi, r).simplified(&zi, r).is_id(),
                "compute_inverse_generator returned wrong inverse {q} of {p}");
        }
    }
}