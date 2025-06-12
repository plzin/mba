use super::*;
use super::mod_inv_pow_two::ModInv as _;
use crate::formatter::Formatter;

// Implements an order on a ring whose elements are represented using
// `BigUint`s, so currently `BigIntModN` and `BitwiseBigInt`.
// TODO: We could also implement this as if we are using the representatives
// `-floor(n/2)..=ceil(n/2)`. In that case we could rename `BigIntModN` to
// `BigUintModN` with the current order an a `BigIntModN` with the other order.
macro_rules! impl_ordered_ring_uint {
    ($ring:ident) => {
        impl OrderedRing for $ring {
            fn cmp(
                &self,
                l: &Self::Element,
                r: &Self::Element,
            ) -> std::cmp::Ordering {
                l.cmp(r)
            }

            fn is_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l < r
            }

            fn is_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l <= r
            }

            fn is_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l >= r
            }

            fn is_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l > r
            }

            fn is_positive(&self, e: &Self::Element) -> bool {
                !RingElement::is_zero(e)
            }

            fn is_negative(&self, _: &Self::Element) -> bool {
                false
            }

            fn abs_assign(&self, _: &mut Self::Element) {}

            fn abs(&self, e: Self::Element) -> Self::Element {
                e
            }

            fn cmp_abs(
                &self,
                l: &Self::Element,
                r: &Self::Element,
            ) -> std::cmp::Ordering {
                l.cmp(r)
            }

            fn is_abs_eq(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l == r
            }

            fn is_abs_ne(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l != r
            }

            fn is_abs_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l < r
            }

            fn is_abs_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l <= r
            }

            fn is_abs_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l >= r
            }

            fn is_abs_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l > r
            }
        }
    }
}

pub(crate) use impl_ordered_ring_uint;

macro_rules! uint_ring {
    ($ring:ident, $uint:ident, $c_type:literal) => {
        #[derive(Clone, PartialEq, Eq, Debug)]
        pub struct $ring;

        impl_ring_element!($uint);

        impl Ring for $ring {
            type Element = $uint;

            fn is_domain(&self) -> bool {
                false
            }

            fn negative_one(&self) -> Self::Element {
                $uint::MAX
            }

            fn neg_assign(&self, e: &mut Self::Element) {
                *e = e.wrapping_neg();
            }

            fn add_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_add(*r);
            }

            fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_sub(*r);
            }

            fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_mul(*r);
            }

            fn is_unit(&self, e: &Self::Element) -> bool {
                e.is_odd()
            }

            fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
                self.is_unit(e).then(|| self.mod_inv(e))
            }

            fn is_zero_divisor(&self, e: &Self::Element) -> bool {
                e.is_even()
            }

            fn random<R: rand::Rng>(&self, rand: &mut R) -> Self::Element {
                rand.random()
            }

            fn element_from_usize(&self, n: usize) -> Self::Element {
                n as $uint
            }

            #[cfg(target_pointer_width = "32")]
            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                let mut e = 0;
                for (i, digit) in n.iter_u32_digits()
                    .take(($uint::BITS.div_ceil(32)) as usize)
                    .enumerate()
                {
                    e |= (digit as $uint) << (i * 32);
                }
                e
            }

            #[cfg(target_pointer_width = "64")]
            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                let mut e = 0;
                for (i, digit) in n.iter_u64_digits()
                    .take(($uint::BITS.div_ceil(64)) as usize)
                    .enumerate()
                {
                    e |= (digit as $uint) << (i * 64);
                }
                e
            }

            fn data_type_name(&self, formatter: Formatter) -> impl std::fmt::Display {
                match formatter {
                    Formatter::C => $c_type,
                    Formatter::Rust => concat!("Wrapping<", stringify!($uint), ">"),
                    Formatter::Tex => $c_type,
                }
            }
        }

        impl BinaryRing for $ring {
            fn bits(&self) -> u32 {
                $uint::BITS
            }

            fn bit(e: &Self::Element, i: u32) -> bool {
                *e & (1 << i) != 0
            }

            fn not_assign(&self, e: &mut Self::Element) {
                *e = !*e;
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
                *l <<= r;
            }

            fn count_ones(e: &Self::Element) -> u32 {
                e.count_ones()
            }

            fn min_bits(e: &Self::Element) -> u32 {
                $uint::BITS - e.leading_zeros()
            }

            fn to_representative(&e: &Self::Element) -> BigUint {
                e.into()
            }

            fn to_usize(&e: &Self::Element) -> usize {
                e as usize
            }
        }

        impl_ordered_ring_uint!($ring);

        impl IntDivRing for $ring {
            fn rounded_div(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                // This isn't my favorite way to do this. But there are problems
                // with the naive implementation (`(l + (r >> 1)) / r`), i.e.
                // it can easily overflow. This can be avoided by casting to
                // bigger types, but this is also not really ideal as there is
                // no bigger type than `u128`, but might be faster, for the
                // other types. The same is true about casting to floats.

                // One x86 at least the quotient and remainder can be computed
                // in a single instruction (`div`).
                let q = l / r;
                let rem = l % r;
                q + if rem >= r.div_ceil(2) {
                    1
                } else {
                    0
                }
            }

            fn euclidean_div(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                l.overflowing_div(r).0
            }

            fn euclidean_rem(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                l.overflowing_rem(r).0
            }
        }
    };
}

uint_ring!(U8, u8, "uint8_t");
uint_ring!(U16, u16, "uint16_t");
uint_ring!(U32, u32, "uint32_t");
uint_ring!(U64, u64, "uint64_t");
uint_ring!(U128, u128, "uint128_t");


#[cfg(test)]
mod test_primitive_uint {
    use super::*;
    use super::test::*;

    #[test]
    fn test_inverse_u8() {
        test_inverse(&U8);
    }

    #[test]
    fn test_inverse_u16() {
        test_inverse(&U16);
    }

    #[test]
    fn test_inverse_u32() {
        test_inverse(&U32);
    }

    #[test]
    fn test_inverse_u64() {
        test_inverse(&U64);
    }

    #[test]
    fn test_inverse_u128() {
        test_inverse(&U128);
    }

    #[test]
    fn test_element_from_biguint_u8() {
        test_element_from_biguint(&U8);
    }

    #[test]
    fn test_element_from_biguint_u16() {
        test_element_from_biguint(&U16);
    }

    #[test]
    fn test_element_from_biguint_u32() {
        test_element_from_biguint(&U32);
    }

    #[test]
    fn test_element_from_biguint_u64() {
        test_element_from_biguint(&U64);
    }

    #[test]
    fn test_element_from_biguint_u128() {
        test_element_from_biguint(&U128);
    }

    #[test]
    fn test_rounded_div_u8() {
        test_rounded_div(&U8);
    }

    #[test]
    fn test_rounded_div_u16() {
        test_rounded_div(&U16);
    }

    #[test]
    fn test_rounded_div_u32() {
        test_rounded_div(&U32);
    }

    #[test]
    fn test_rounded_div_u64() {
        test_rounded_div(&U64);
    }

    #[test]
    fn test_rounded_div_u128() {
        test_rounded_div(&U128);
    }

    fn test_rounded_div_special_cases<R>(r: &R)
    where
        R: BinaryRing,
        R::Element: Copy + std::ops::Div<Output = R::Element>,
    {
        let max = r.sub(R::zero(), &R::one());
        let half_floor = max / r.element_from_usize(2);
        let half_ceil = r.inc(half_floor);
        let pairs = [
            (R::zero(), R::one()),
            (R::zero(), half_floor),
            (R::zero(), half_ceil),
            (R::zero(), max),
            (R::one(), R::one()),
            (R::one(), half_floor),
            (R::one(), half_ceil),
            (R::one(), max),
            (half_floor, R::one()),
            (half_floor, half_floor),
            (half_floor, half_ceil),
            (half_floor, max),
            (half_ceil, R::one()),
            (half_ceil, half_floor),
            (half_ceil, half_ceil),
            (half_ceil, max),
            (max, R::one()),
            (max, half_floor),
            (max, half_ceil),
            (max, max),
        ];

        for (a, b) in pairs {
            let q = R::rounded_div(&a, &b);
            let rational = BigRational::new(
                R::to_representative(&a).into(),
                R::to_representative(&b).into(),
            );
            let check = r.element_from_biguint(rational.round().numer().magnitude());
            assert_eq!(q, check,
                "rounded_div({a}, {b}) = {check} but got {q}");
        }
    }

    #[test]
    fn test_rounded_div_special_cases_u8() {
        test_rounded_div_special_cases(&U8);
    }

    #[test]
    fn test_rounded_div_special_cases_u16() {
        test_rounded_div_special_cases(&U16);
    }

    #[test]
    fn test_rounded_div_special_cases_u32() {
        test_rounded_div_special_cases(&U32);
    }

    #[test]
    fn test_rounded_div_special_cases_u64() {
        test_rounded_div_special_cases(&U64);
    }

    #[test]
    fn test_rounded_div_special_cases_u128() {
        test_rounded_div_special_cases(&U128);
    }
}