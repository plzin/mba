use super::*;
use crate::formatter::Formatter;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Z;

impl Ring for Z {
    type Element = BigInt;

    fn is_domain(&self) -> bool {
        true
    }

    fn negative_one(&self) -> Self::Element {
        -<BigInt as One>::one()
    }

    fn neg_assign(&self, e: &mut Self::Element) {
        neg_assign(e);
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

    fn inc_assign(&self, e: &mut Self::Element) {
        e.inc();
    }

    fn dec_assign(&self, e: &mut Self::Element) {
        e.dec();
    }

    fn is_unit(&self, e: &Self::Element) -> bool {
        RingElement::is_one(e.magnitude())
    }

    fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
        if !self.is_unit(e) {
            None
        } else if e.is_negative() {
            Some(self.negative_one())
        } else {
            Some(Self::one())
        }
    }

    fn is_zero_divisor(&self, e: &Self::Element) -> bool {
        RingElement::is_zero(e)
    }

    fn random<R: rand::Rng>(&self, _: &mut R) -> Self::Element {
        panic!("Can't sample integer uniformly at random.");
    }

    fn element_from_usize(&self, n: usize) -> Self::Element {
        BigInt::from(n)
    }

    fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
        n.clone().into()
    }

    fn data_type_name(&self, _formatter: Formatter) -> impl std::fmt::Display {
        "BigInt"
    }
}

impl SqrtRing for Z {
    fn sqrt(e: &Self::Element) -> Self::Element {
        e.sqrt()
    }
}

impl IntDivRing for Z {
    fn rounded_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
        // rug has a function for this.
        if l.sign() == r.sign() {
            (l + (r / 2)) / r
        } else {
            (l - (r / 2)) / r
        }
    }

    fn euclidean_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
        l.div_euclid(r)
    }

    fn euclidean_rem(l: &Self::Element, r: &Self::Element) -> Self::Element {
        l.rem_euclid(r)
    }
}

impl OrderedRing for Z {
    fn cmp(&self, l: &Self::Element, r: &Self::Element) -> std::cmp::Ordering {
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

    fn cmp_abs(&self, l: &Self::Element, r: &Self::Element) -> std::cmp::Ordering {
        l.magnitude().cmp(r.magnitude())
    }

    fn is_abs_eq(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.magnitude() == r.magnitude()
    }

    fn is_abs_ne(&self, l: &Self::Element, r: &Self::Element) -> bool {
        !self.is_abs_eq(l, r)
    }

    fn is_abs_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.magnitude() < r.magnitude()
    }

    fn is_abs_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.magnitude() <= r.magnitude()
    }

    fn is_abs_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.magnitude() >= r.magnitude()
    }

    fn is_abs_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        l.magnitude() > r.magnitude()
    }
}

#[test]
fn rounded_div_test() {
    use rand::SeedableRng as _;
    use rand::distr::{Distribution as _, Uniform};
    let mut rng = rand::rngs::StdRng::seed_from_u64(0);

    let dividend_max = 1000;
    let dividend_dist = Uniform::new(-dividend_max, dividend_max + 1).unwrap();
    let divisor_max = 100;
    let divisor_dist = Uniform::new(-divisor_max, divisor_max + 1).unwrap();

    let iterations = 1000usize;
    for _ in 0..iterations {
        let dividend = dividend_dist.sample(&mut rng);
        let divisor = loop {
            let divisor = divisor_dist.sample(&mut rng);
            if divisor != 0 {
                break divisor;
            }
        };

        let quotient = Z::rounded_div(&dividend.into(), &divisor.into());
        let check = (((dividend as f64) / (divisor as f64)).round() as i32).into();
        assert_eq!(
            quotient, check,
            "rounded_div({dividend}, {divisor}) = {check} but got {quotient}"
        );
    }
}
