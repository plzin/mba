use super::*;
use crate::formatter::Formatter;

macro_rules! float_field {
    (
        $field:ident,
        $float:ident,
        $from_big:ident,
        $c_type:literal,
        $rust_type:literal
    ) => {
        /// Floating point numbers aren't actually a ring (or a field) for
        /// multiple reasons:
        /// - Addition and multiplication are not associative
        /// - Infinities, NaNs, -0
        ///
        /// But we act like it so we can store them in vectors, matrices.
        ///
        /// Many functions assume that the floats are not NaN or infinity. Your
        /// computer won't blow up when that happens, but some methods may
        /// return unexpected results. E.g. `is_unit` only checks for zero.
        #[derive(Clone, PartialEq, Eq, Debug)]
        pub struct $field;

        impl_ring_element!($float);

        impl Ring for $field {
            type Element = $float;

            fn is_domain(&self) -> bool {
                true
            }

            fn negative_one(&self) -> Self::Element {
                -1.
            }

            fn neg_assign(&self, e: &mut Self::Element) {
                *e = -*e;
            }

            fn add_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l += r;
            }

            fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l -= r;
            }

            fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l *= r;
            }

            fn is_unit(&self, e: &Self::Element) -> bool {
                !RingElement::is_zero(e)
            }

            fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
                self.is_unit(e).then(|| e.recip())
            }

            fn is_zero_divisor(&self, e: &Self::Element) -> bool {
                *e == 0.
            }

            fn random<R: rand::Rng>(&self, _: &mut R) -> Self::Element {
                panic!("Can't sample float uniformly at random.");
            }

            fn element_from_usize(&self, n: usize) -> Self::Element {
                n as $float
            }

            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                use num_traits::cast::ToPrimitive;
                n.$from_big().unwrap()
            }

            fn data_type_name(
                &self,
                formatter: Formatter,
            ) -> impl std::fmt::Display {
                match formatter {
                    Formatter::C => $c_type,
                    Formatter::Rust => $rust_type,
                    Formatter::Tex => $c_type,
                }
            }
        }

        impl SqrtRing for $field {
            fn sqrt(e: &Self::Element) -> Self::Element {
                e.sqrt()
            }
        }

        impl Field for $field {
            fn div_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l /= r;
            }
        }

        impl_int_div_for_field!($field);

        impl OrderedRing for $field {
            fn cmp(
                &self,
                l: &Self::Element,
                r: &Self::Element,
            ) -> std::cmp::Ordering {
                l.partial_cmp(r).unwrap()
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
                *e >= 0.
            }

            fn is_negative(&self, e: &Self::Element) -> bool {
                *e < 0.
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
                l.abs().partial_cmp(&r.abs()).unwrap()
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
    };
}

float_field!(F32, f32, to_f32, "float", "f32");
float_field!(F64, f64, to_f64, "double", "f64");
