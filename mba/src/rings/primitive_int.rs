use super::*;
use super::mod_inv_pow_two::ModInv as _;
use crate::formatter::Formatter;

macro_rules! int_ring {
    ($ring:ident, $int:ident, $uring:ident, $c_type:literal) => {
        #[derive(Clone, PartialEq, Eq, Debug)]
        pub struct $ring;

        impl_ring_element!($int);

        impl Ring for $ring {
            type Element = $int;

            fn is_domain(&self) -> bool {
                false
            }

            fn negative_one(&self) -> Self::Element {
                -1
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
                self.is_unit(e).then(|| $uring.mod_inv(&e.cast_unsigned()).cast_signed())
            }

            fn is_zero_divisor(&self, e: &Self::Element) -> bool {
                e.is_even()
            }

            fn random<R: rand::Rng>(&self, rand: &mut R) -> Self::Element {
                rand.random()
            }

            fn element_from_usize(&self, n: usize) -> Self::Element {
                n as $int
            }

            #[cfg(target_pointer_width = "32")]
            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                let mut e = 0;
                for (i, digit) in n.iter_u32_digits()
                    .take(($int::BITS.div_ceil(32)) as usize)
                    .enumerate()
                {
                    e |= (digit as $int) << (i * 32);
                }
                e
            }

            #[cfg(target_pointer_width = "64")]
            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                let mut e = 0;
                for (i, digit) in n.iter_u64_digits()
                    .take(($int::BITS.div_ceil(64)) as usize)
                    .enumerate()
                {
                    e |= (digit as $int) << (i * 64);
                }
                e
            }

            fn data_type_name(&self, formatter: Formatter) -> impl std::fmt::Display {
                match formatter {
                    Formatter::C => $c_type,
                    Formatter::Rust => concat!("Wrapping<", stringify!($int), ">"),
                    Formatter::Tex => $c_type,
                }
            }
        }

        impl BinaryRing for $ring {
            fn bits(&self) -> u32 {
                $int::BITS
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
                $int::BITS - e.leading_zeros()
            }

            fn to_representative(e: &Self::Element) -> BigUint {
                e.cast_unsigned().into()
            }

            fn to_usize(&e: &Self::Element) -> usize {
                e as usize
            }
        }

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
                *e > 0
            }

            fn is_negative(&self, e: &Self::Element) -> bool {
                *e < 0
            }

            fn abs_assign(&self, e: &mut Self::Element) {
                *e = e.abs();
            }

            fn abs(&self, e: Self::Element) -> Self::Element {
                e.abs()
            }

            fn cmp_abs(
                &self,
                l: &Self::Element,
                r: &Self::Element,
            ) -> std::cmp::Ordering {
                l.abs().cmp(&r.abs())
            }

            fn is_abs_eq(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() == r.abs()
            }

            fn is_abs_ne(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() != r.abs()
            }

            fn is_abs_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() < r.abs()
            }

            fn is_abs_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() <= r.abs()
            }

            fn is_abs_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() >= r.abs()
            }

            fn is_abs_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
                l.abs() > r.abs()
            }
        }

        impl IntDivRing for $ring {
            fn rounded_div(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                // See `primitive_uint.rs` for why this seems more complicated
                // than it should be. Note also that in this case `i32::MIN`
                // divided by `-1` is `i32::MIN` instead of panicking.

                let q = l.overflowing_div(r).0;
                let rem = l.overflowing_rem(r).0;

                // Let's hope the compiler implements the `div_ceil` as just
                // `(abs(r) + 1) >> 1` because it doesn't have to worry about
                // `abs(r) + 1` overflowing.
                q + if rem.unsigned_abs() >= r.unsigned_abs().div_ceil(2) {
                    // This is basically `l.signum() * r.signum()`, but more
                    // efficient in the case that `l` and `r` are non-zero.
                    ((l ^ r) >> ($int::BITS - 1)) | 1
                } else {
                    0
                }
            }

            fn euclidean_div(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                l.overflowing_div_euclid(r).0
            }

            fn euclidean_rem(
                &l: &Self::Element,
                &r: &Self::Element,
            ) -> Self::Element {
                l.overflowing_rem_euclid(r).0
            }
        }
    };
}

int_ring!(I8, i8, U8, "int8_t");
int_ring!(I16, i16, U16, "int16_t");
int_ring!(I32, i32, U32, "int32_t");
int_ring!(I64, i64, U64, "int64_t");
int_ring!(I128, i128, U128, "int128_t");

#[cfg(test)]
mod test_primitive_int {
    use super::*;
    use super::test::*;

    #[test]
    fn test_inverse_i8() {
        test_inverse(&I8);
    }

    #[test]
    fn test_inverse_i16() {
        test_inverse(&I16);
    }

    #[test]
    fn test_inverse_i32() {
        test_inverse(&I32);
    }

    #[test]
    fn test_inverse_i64() {
        test_inverse(&I64);
    }

    #[test]
    fn test_inverse_i128() {
        test_inverse(&I128);
    }

    #[test]
    fn test_element_from_biguint_i8() {
        test_element_from_biguint(&I8);
    }

    #[test]
    fn test_element_from_biguint_i16() {
        test_element_from_biguint(&I16);
    }

    #[test]
    fn test_element_from_biguint_i32() {
        test_element_from_biguint(&I32);
    }

    #[test]
    fn test_element_from_biguint_i64() {
        test_element_from_biguint(&I64);
    }

    #[test]
    fn test_element_from_biguint_i128() {
        test_element_from_biguint(&I128);
    }

    #[test]
    fn test_rounded_div_i8() {
        test_rounded_div(&I8);
    }

    #[test]
    fn test_rounded_div_i16() {
        test_rounded_div(&I16);
    }

    #[test]
    fn test_rounded_div_i32() {
        test_rounded_div(&I32);
    }

    #[test]
    fn test_rounded_div_i64() {
        test_rounded_div(&I64);
    }

    #[test]
    fn test_rounded_div_i128() {
        test_rounded_div(&I128);
    }

    fn test_rounded_div_special_cases<R>(
        r: &R,
        min: R::Element,
        max: R::Element,
    )
    where
        R: BinaryRing,
        R::Element: Copy + std::ops::Div<Output = R::Element>,
        BigInt: From<R::Element>,
    {
        let pairs = [
            (R::zero(), R::one()),
            (R::zero(), min),
            (R::zero(), max),
            (R::one(), R::one()),
            (R::one(), min),
            (R::one(), max),
            (min, R::one()),
            (min, min),
            (min, max),
            (max, R::one()),
            (max, min),
            (max, max),
        ];

        for (a, b) in pairs {
            for i in 0..4 {
                // Check all combinations of positive and negative `a` and `b`.
                let a = if i & 1 == 0 { a } else { r.neg(a) };
                let b = if i & 2 == 0 { b } else { r.neg(b) };

                let q = R::rounded_div(&a, &b);
                let check = r.element_from_bigint(BigRational::new(
                    BigInt::from(a),
                    BigInt::from(b),
                ).round().numer());
                assert_eq!(q, check,
                    "rounded_div({a}, {b}) = {check} but got {q}");
            }
        }
    }

    #[test]
    fn test_rounded_div_special_cases_i8() {
        test_rounded_div_special_cases(&I8, i8::MIN, i8::MAX);
    }

    #[test]
    fn test_rounded_div_special_cases_i16() {
        test_rounded_div_special_cases(&I16, i16::MIN, i16::MAX);
    }

    #[test]
    fn test_rounded_div_special_cases_i32() {
        test_rounded_div_special_cases(&I32, i32::MIN, i32::MAX);
    }

    #[test]
    fn test_rounded_div_special_cases_i64() {
        test_rounded_div_special_cases(&I64, i64::MIN, i64::MAX);
    }

    #[test]
    fn test_rounded_div_special_cases_i128() {
        test_rounded_div_special_cases(&I128, i128::MIN, i128::MAX);
    }
}