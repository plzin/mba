#![allow(dead_code)]

use rug::{Integer, Complete, rand::{ThreadRandState, ThreadRandGen}};
use rand::{distributions::{Distribution, Uniform}, prelude::ThreadRng, RngCore};

use crate::poly::Poly;

struct Generator(ThreadRng);
impl Generator {
    pub fn new(rng: ThreadRng) -> Self {
        Self(rng)
    }
}

impl ThreadRandGen for Generator {
    fn gen(&mut self) -> u32 {
        self.0.next_u32()
    }
}


/// Returns a pair of permutation polynomials mod 2^n.
/// The functions are inverses of each other.
/// The degree of one polynomial can be specified.
pub fn perm_pair(qr: &QuotientRing, degree: usize) -> (Poly, Poly) {
    assert!(degree >= 1, "Can't create a permutation polynomial \
            of degree 0 as this would be a constant.");
    // Generate a random permutation polynomial.
    let p = random_perm_poly(qr, degree);
    assert!(is_perm_poly(&p), "Generated an invalid polynomial.");

    // Find the inverse.
    let q = compute_inverse(&p, qr);

    (p, q)

    // let p = Poly::from_vec(vec![3, 3, 5, 1, 7, 9]);
    // let q = self.compute_inverse(&p);
    // println!("p(x) = {}", p);
    // println!("q(x) = {}", q);
    // let mut r = self.compose(&p, &q);
    // self.simp_poly(&mut r);
    // println!("p(q(x)) = {}", r);

    // let p = Poly::from_vec(vec![0, 1, 2]);
    // let q = Poly::from_vec(vec![0, 1, 0, 1]);

    // let r = self.compose(&p, &q);

    // println!("p(x) = {}", p);
    // println!("q(x) = {}", q);
    // println!("p(q(x)) = {}", r);

    // let mut p = Poly::from_vec(vec![3, 0, 8, 2, 15, 5, 3, 1]);
    // println!("{}", p);

    // for i in 0u32..(1<<n) {
    //     println!("p({}) = {}", i, p.eval(&i.into()).keep_bits(n));
    // }

    // simp_poly(&mut p, &zero_ideal, n);

    // println!("{}", p);
    // for i in 0u32..(1<<n) {
    //     println!("p({}) = {}", i, p.eval(&i.into()).keep_bits(n));
    // }


    // Verify that all generators do indeed compute the zero polynomial.
    // for (i, p) in zero_ideal.iter().enumerate() {
    //     for j in 0u32..1 << n {
    //         let r = p.eval(&j.into()).keep_bits(n);
    //         println!("p_{}({}) = {}", i, j, r);
    //     }
    // }
}

/// Returns a random permutation polynomial.
fn random_perm_poly(qr: &QuotientRing, degree: usize) -> Poly {
    let mut gen = Generator::new(rand::thread_rng());
    let mut rng = ThreadRandState::new_custom(&mut gen);
    let mut p: Vec<_> =  (0..=degree)
        .map(|_| Integer::from(Integer::random_bits(qr.n, &mut rng)))
        .collect();

    // Make sure that this is a permutation polynomial.
    // The coefficient of degree 1 has to be odd.
    if p[1].is_even() {
        p[1] += 1;
        p[1].keep_bits_mut(qr.n);
    }

    let mut rng = rand::thread_rng();

    // Make sure the sum of the even coefficients (except 0) is even.
    if p.iter().skip(2).step_by(2).fold(false, parity) {
        let dist = Uniform::from(1..=degree/2);
        let i = dist.sample(&mut rng);
        p[2*i] += 1;
        p[2*i].keep_bits_mut(qr.n);
    }

    // Make sure the sum of the odd coefficients (except 1) is even.
    if p.iter().skip(3).step_by(2).fold(false, parity) {
        let dist = Uniform::from(1..=(degree-1)/2);
        let i = dist.sample(&mut rng);
        p[2*i+1] += 1;
        p[2*i+1].keep_bits_mut(qr.n);
    }

    Poly::from_vec(p)
}

/// Computes the composition of two polynomials.
fn compose(p: &Poly, q: &Poly, qr: &QuotientRing) -> Poly {
    // Iterate over the coefficients in reverse order.
    let mut iter = p.coeffs.iter().rev();

    let mut r = Poly::constant(
        iter.next().map_or_else(|| Integer::new(), |c| c.clone())
    );

    // The last coefficient is the initial value.
    for c in iter {
        // Multiply the current value by a and add the next coefficient.
        r *= q;
        r += c;
        r.reduce(qr);
    }

    r
}

/// Computes the inverse of a permutation polynomial.
fn compute_inverse(f: &Poly, qr: &QuotientRing) -> Poly {
    assert!(is_perm_poly(&f),
        "Can't invert the function as it is not a permutation polynomial");
    // Simplify p.
    let p = f.clone().simplified(qr);

    // Initialize q with the initial guess: p(x) = x.
    let mut q = Poly::from_vec(vec![0, 1]);

    let mut it = 0;

    // Do the newton method.
    loop {
        // Not proven to always work so make sure
        // to stop after a certain number of iterations.
        assert!(it <= qr.n * 2, "Failed to compute the inverse \
                in a reasonable number of iterations.");
        //println!("compute_inverse: it={}: {}", it, q);
        // Compute the composition.
        let mut comp = compose(&p, &q, qr)
            .simplified(qr);

        // Do we already have p(q(x)) = x?
        if &comp.coeffs == &[0, 1] {
            //println!("compute_inverse: Found the inverse after {} iterations.", it);
            return q;
        }

        // Subtract x.
        // This is the quantity we want to make 0.
        comp.coeffs[1] -= 1;

        // Update the guess.
        let qd = q.derivative();
        q -= &(&qd * &comp);
        q.simplify(qr);

        it += 1;
    }
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
    for i in 1u32.. {
        div += i.trailing_zeros();

        // If the exponent would be negative
        // then add the last generator and stop.
        if n <= i + div {
            let mut p = Poly::one();

            for j in 0..2*i {
                // Multiply the current polynomial by (x-j).
                p.mul_linfac(&j.into());
            }

            p.mod_coeff(n);

            gen.push(p);
            break;
        }

        // Compute the exponent.
        let e = n - i - div;

        // Let's build the polynomial.
        let mut p = Poly::one();

        for j in 0..2*i {
            // Multiply the current polynomial by (x-j).
            p.mul_linfac(&j.into());
        }

        p.shl_coeff(e);
        p.mod_coeff(n);

        gen.push(p);
    }

    gen
}

/// Used internally as a fold function.
/// Computes the "parity" of a list of integers,
/// that is whether they are even or odd,
/// depending on the intial value of the accumulator.
fn parity(acc: bool, i: &Integer) -> bool {
    match i.is_odd() {
        true => !acc,
        false => acc,
    }
}

/// Is this function a permutation polynomial?
fn is_perm_poly(f: &Poly) -> bool {
    return f.coeffs.get(1).map_or(false, |i| i.is_odd())
        && f.coeffs.iter().skip(2).step_by(2).fold(true, parity)
        && f.coeffs.iter().skip(3).step_by(2).fold(true, parity);
}

/// Stores some information about the quotiont ring.
/// Which are annoying to pass to all functions.
pub struct QuotientRing {
    /// The coefficients are mod 2^n.
    n: u32,

    /// The generators of the "zero ideal".
    /// The zero ideal contains all polynomials p s.t. p(x) = 0 for all x.
    zi: Vec<Poly>,
}

impl QuotientRing {
    /// Initializes the quotient ring.
    /// This ring polynomials from Z/2^nZ[x] mod some simplifications.
    pub fn init(n: u32) -> Self {
        assert!(n > 0, "Not a valid ring.");
        let zi = zero_ideal(n);
        //for (i, g) in zi.iter().enumerate() {
        //    println!("{}: {}", i, g);
        //}

        Self { n, zi }
    }

}

impl Poly {
    fn simplified(mut self, qr: &QuotientRing) -> Self {
        self.simplify(qr);
        self
    }

    /// Simplifies a polynomial by adding a polynomial in the zero ideal
    /// to reduce the degree of the polynomial.
    fn simplify(&mut self, qr: &QuotientRing) {
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

        for gen in qr.zi.iter().rev() {
            let gen_len = gen.len();

            while coeff + 1 >= gen_len {
                let m = (&self.coeffs[coeff] / &gen.coeffs[gen_len-1])
                    .complete();
                if m != 0 {
                    let iter = (&mut self.coeffs[coeff+1-gen_len..=coeff])
                        .iter_mut().zip(gen.coeffs.iter());

                    for (p, g) in iter {
                        *p -= &m * g;
                    }
                }
                coeff -= 1;
            }
        }

        // Reduce all coefficients mod 2^n.
        self.mod_coeff(qr.n);

        // Remove leading coefficients that are 0.
        self.truncate();
    }

    /// Reduces the degree of a polynomial if possible.
    /// This function should be used during computations
    /// where the degree could otherwise explode.
    /// It is a simplified version of the simplify algorithm,
    /// that only uses the generator of the highest degree.
    fn reduce(&mut self, qr: &QuotientRing) {
        if self.len() == 0 {
            return;
        }

        let gen = qr.zi.last().unwrap();
        let gen_len = gen.len();
        while self.len() >= gen_len {
            let p_len = self.len();
            let (c, rest) = self.coeffs.split_last_mut().unwrap();
            if *c != 0 {
                let iter = (&mut rest[p_len-gen_len..]).iter_mut()
                    .zip(gen.coeffs.iter());

                for (p, g) in iter {
                    *p -= &*c * g;
                }
            }
            self.coeffs.pop();
        }

        // Reduce all coefficients mod 2^n.
        self.mod_coeff(qr.n);

        // Remove leading coefficients that are 0.
        self.truncate();
    }
}


#[test]
fn check_inverse_8() {
    let qr = QuotientRing::init(8);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&qr, d);
        let comp = compose(&p, &q, &qr).simplified(&qr);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
fn check_inverse_16() {
    let qr = QuotientRing::init(8);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&qr, d);
        let comp = compose(&p, &q, &qr).simplified(&qr);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
fn check_inverse_32() {
    let qr = QuotientRing::init(32);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&qr, d);
        let comp = compose(&p, &q, &qr).simplified(&qr);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}

#[test]
pub fn check_inverse_64() {
    let qr = QuotientRing::init(64);
    let deg = Uniform::from(1..16);
    let mut rng = rand::thread_rng();
    for _ in 0..100 {
        let d = deg.sample(&mut rng);
        let (p, q) = perm_pair(&qr, d);
        let comp = compose(&p, &q, &qr).simplified(&qr);
        assert!(comp.is_id(), "Incorrect inverse.");
    }
}
