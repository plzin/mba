use super::*;
use crate::formatter::Formatter;

/// The integers mod 2^n. So esentially `n` bit integers with wrapping
/// semantics. The elements are stored as [`BigUint`]s. The arithmetic could
/// definitely be more efficient if [`num_bigint`] supported this use case
/// better. We usually calculate the full result and then truncate to the number
/// of bits, instead of ignoring all bits above the limit during the
/// computation.
#[derive(Clone, PartialEq, Eq)]
pub struct BinaryBigInt {
    /// The number of bits of an element or log2 of the modulus if you will.
    bits: u32,

    /// We are storing the modulus here because it lets us avoid allocations
    /// for a lot of operations, e.g. `negative_one`, `neg_assign`.
    /// This wouldn't be necessary if `num_bigint` supported this use case
    /// better.
    modulus: BigUint,
}

impl BinaryBigInt {
    pub fn new(bits: u32) -> Self {
        assert!(bits > 0, "bits must be greater than 0");
        Self { bits, modulus: biguint_pow2(bits) }
    }

    /// Returns the number of bits of a ring element.
    pub fn bits(&self) -> u32 {
        self.bits
    }

    /// Returns the modulus.
    pub fn modulus(&self) -> &BigUint {
        &self.modulus
    }
}

impl Ring for BinaryBigInt {
    type Element = BigUint;

    fn is_domain(&self) -> bool {
        self.bits == 1
    }

    fn negative_one(&self) -> Self::Element {
        self.modulus() - 1u8
    }

    fn neg_assign(&self, e: &mut Self::Element) {
        neg_assign_mod(e, self.modulus());
    }

    fn add_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        *l += r;
        l.truncate(self.bits as u64);
    }

    fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        if r > l {
            l.set_bit(self.bits as u64, true);
        }
        *l -= r;
    }

    fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        *l *= r;
        l.truncate(self.bits as u64);
    }

    fn inc_assign(&self, e: &mut Self::Element) {
        e.inc();
        if e == &self.modulus {
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
        e.is_odd()
    }

    fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
        // TODO: This can be done way more efficiently. Computing the
        // multiplicative inverse mod a power of two can be done with Newton's
        // method.

        // This check is a fast fail.
        if !self.is_unit(e) {
            return None;
        }

        e.modinv(&self.modulus)
    }

    fn is_zero_divisor(&self, e: &Self::Element) -> bool {
        e.is_even()
    }

    fn random<R: rand::Rng>(&self, rng: &mut R) -> Self::Element {
        rng.random_biguint(self.bits as u64)
    }

    fn element_from_usize(&self, n: usize) -> Self::Element {
        BigUint::from(n).truncated(self.bits as u64)
    }

    fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
        n.clone().truncated(self.bits as u64)
    }

    fn data_type_name(&self, formatter: Formatter) -> impl std::fmt::Display {
        BinaryBigIntDataType { bits: self.bits, formatter }
    }
}

impl BinaryRing for BinaryBigInt {
    fn bits(&self) -> u32 {
        self.bits
    }

    fn bit(e: &Self::Element, i: u32) -> bool {
        e.bit(i as u64)
    }

    fn not_assign(&self, e: &mut Self::Element) {
        not_assign_mod(e, self.modulus());
    }

    fn and_assign(l: &mut Self::Element, r: &Self::Element) {
        *l &= r;
    }

    fn or_assign(l: &mut Self::Element, r: &Self::Element) {
        *l |= r;
    }

    fn xor_assign(l: &mut Self::Element, r: &Self::Element) {
        *l ^= r;
    }

    fn shl_assign(&self, l: &mut Self::Element, r: u32) {
        if r >= self.bits {
            l.truncate(0);
        } else {
            *l <<= r;
            l.truncate(self.bits as u64);
        }
    }

    fn count_ones(e: &Self::Element) -> u32 {
        e.count_ones() as u32
    }

    fn min_bits(e: &Self::Element) -> u32 {
        e.bits() as u32
    }

    fn to_representative(e: &Self::Element) -> BigUint {
        e.clone()
    }

    #[cfg(target_pointer_width = "32")]
    fn to_usize(e: &Self::Element) -> usize {
        e.iter_u32_digits().next().unwrap_or(0) as usize
    }

    #[cfg(target_pointer_width = "64")]
    fn to_usize(e: &Self::Element) -> usize {
        e.iter_u64_digits().next().unwrap_or(0) as usize
    }
}

impl_ordered_ring_uint!(BinaryBigInt);

impl IntDivRing for BinaryBigInt {
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

impl std::fmt::Debug for BinaryBigInt {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "BitwiseBigInt({})", self.bits)
    }
}

pub struct BinaryBigIntDataType {
    bits: u32,
    formatter: Formatter,
}

impl std::fmt::Display for BinaryBigIntDataType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.formatter {
            Formatter::C => write!(f, "uint{}_t", self.bits),
            Formatter::Rust => write!(f, "Wrapping<u{}>", self.bits),
            Formatter::Tex => write!(f, "uint{}", self.bits),
            Formatter::LLVM => write!(f, "i{}", self.bits),
        }
    }
}

#[cfg(test)]
mod test_bitwise_bigint {
    use rand::{Rng, SeedableRng, rngs::StdRng};

    use super::{test::*, *};
    #[test]
    fn test_inverse_bitwise_bigint() {
        let rng = &mut StdRng::seed_from_u64(0);

        for _ in 0..10 {
            let bits: u8 = loop {
                let bits = rng.random();
                if bits > 0 {
                    break bits;
                }
            };
            let r = BinaryBigInt::new(bits as u32);
            test_inverse(&r);
        }
    }
}
