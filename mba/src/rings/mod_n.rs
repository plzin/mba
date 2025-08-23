use super::*;
use crate::formatter::Formatter;

/// The ring of integers mod n. Each element is stored as a [`BigUint`].
#[derive(Clone, PartialEq, Eq)]
pub struct BigIntModN {
    n: BigUint,
}

impl BigIntModN {
    /// Create a new ring.
    pub fn new(n: BigUint) -> Self {
        assert!(n > One::one(), "modulus must be greater than 1");
        Self { n }
    }

    /// Returns the modulus.
    pub fn modulus(&self) -> &BigUint {
        &self.n
    }
}

impl Ring for BigIntModN {
    type Element = BigUint;

    fn is_domain(&self) -> bool {
        // TODO: If the modulus is prime then the ring is a domain.
        // For now we just assume that that the ring isn't a domain.
        false
    }

    fn negative_one(&self) -> Self::Element {
        self.n.clone() - 1u8
    }

    fn neg_assign(&self, e: &mut Self::Element) {
        neg_assign_mod(e, &self.n);
    }

    fn add_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        *l += r;
        reduce_simple(l, &self.n);
    }

    fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        if r > l {
            *l += &self.n;
        }
        *l -= r
    }

    fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        *l *= r;
        *l %= &self.n;
    }

    fn inc_assign(&self, e: &mut Self::Element) {
        e.inc();
        if e == &self.n {
            e.truncate(0);
        }
    }

    fn dec_assign(&self, e: &mut Self::Element) {
        if RingElement::is_zero(e) {
            *e = self.modulus() - 1u32;
        } else {
            e.dec();
        }
    }

    fn is_unit(&self, e: &Self::Element) -> bool {
        RingElement::is_one(&e.gcd(&self.n))
    }

    fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
        e.modinv(&self.n)
    }

    fn is_zero_divisor(&self, e: &Self::Element) -> bool {
        !self.is_unit(e)
    }

    fn random<R: rand::Rng>(&self, rng: &mut R) -> Self::Element {
        rng.random_biguint_below(&self.n)
    }

    fn element_from_usize(&self, n: usize) -> Self::Element {
        BigUint::from(n) % &self.n
    }

    fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
        n.clone() % &self.n
    }

    fn data_type_name(&self, _formatter: Formatter) -> impl std::fmt::Display {
        BigIntModNDataType { modulus: self.n.clone() }
    }
}

impl_ordered_ring_uint!(BigIntModN);

impl IntDivRing for BigIntModN {
    fn rounded_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
        (l + (r >> 1)) / r
    }

    fn euclidean_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
        l / r
    }

    fn euclidean_rem(l: &Self::Element, r: &Self::Element) -> Self::Element {
        l % r
    }
}

impl std::fmt::Debug for BigIntModN {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BigIntModN({})", self.n)
    }
}

pub struct BigIntModNDataType {
    modulus: BigUint,
}

impl std::fmt::Display for BigIntModNDataType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "mod{}_t", self.modulus)
    }
}

#[cfg(test)]
mod test_mod_n {
    use rand::{SeedableRng, rngs::StdRng};

    use super::{test::*, *};

    #[test]
    fn test_inverse_mod_n() {
        let rng = &mut StdRng::seed_from_u64(0);

        for _ in 0..10 {
            let n = loop {
                let n = rng.random_biguint(64);
                if n > One::one() {
                    break n;
                }
            };
            let r = BigIntModN::new(n.clone());
            test_inverse(&r);
        }
    }
}
