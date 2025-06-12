use super::*;
use crate::formatter::Formatter;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Q;

impl Ring for Q {
    type Element = BigRational;

    fn is_domain(&self) -> bool {
        true
    }

    fn negative_one(&self) -> Self::Element {
        (-<BigInt as One>::one()).into()
    }

    fn neg_assign(&self, e: &mut Self::Element) {
        neg_assign(e)
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
        !RingElement::is_zero(e.numer())
    }

    fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
        self.is_unit(e).then(|| e.recip())
    }

    fn is_zero_divisor(&self, e: &Self::Element) -> bool {
        RingElement::is_zero(e.numer())
    }

    fn random<R: rand::Rng>(&self, _: &mut R) -> Self::Element {
        panic!("Can't sample rational number uniformly at random.");
    }

    fn element_from_usize(&self, n: usize) -> Self::Element {
        BigInt::from(n).into()
    }

    fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
        BigInt::from(n.clone()).into()
    }

    fn data_type_name(&self, _formatter: Formatter) -> impl std::fmt::Display {
        "BigRational"
    }
}

impl SqrtRing for Q {
    fn sqrt(e: &Self::Element) -> Self::Element {
        BigRational::new(e.numer().sqrt(), e.denom().sqrt())
    }
}

impl Field for Q {
    fn div_assign(&self, l: &mut Self::Element, r: &Self::Element) {
        *l /= r;
    }
}

impl_int_div_for_field!(Q);

impl OrderedRing for Q {
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
        e.is_positive()
    }

    fn is_negative(&self, e: &Self::Element) -> bool {
        e.is_negative()
    }

    fn abs_assign(&self, e: &mut Self::Element) {
        if e.is_negative() {
            self.neg_assign(e);
        }
    }

    fn cmp_abs(
        &self,
        l: &Self::Element,
        r: &Self::Element,
    ) -> std::cmp::Ordering {
        // Both the `Ratio::cmp` method as well as the code below seem
        // unnecessarily inefficient. We should also do the checks the
        // `Ratio::cmp` does and additionaly compute the number of bits of the
        // result of the multiplications (or just the number of BigInt digits)
        // before doing the whole computation. I think this would completely
        // eleminate the need to multiply and allocate in many practical cases.

        //(l.numer().magnitude() * r.denom().magnitude()).cmp(
        //    &(r.numer().magnitude() * l.denom().magnitude()))

        l.abs().cmp(&r.abs())
    }

    fn is_abs_eq(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.denom() == r.denom() && l.numer().magnitude() == r.numer().magnitude()
    }

    fn is_abs_ne(&self, l: &Self::Element, r: &Self::Element) -> bool {
        !self.is_abs_eq(l, r)
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