use num_bigint::{BigInt, BigUint, Sign};

/// An element of a ring.
/// This exists mostly for convenience, so we can call `e.is_zero` on ring
/// elements, but this should really all be implemented on the `Ring` itself,
/// so we'd have to call `r.is_zero(e)`.
pub trait RingElement: 'static + Clone + PartialEq + std::fmt::Debug + std::fmt::Display {
    /// Returns the "zero" element of the ring.
    fn zero() -> Self;

    /// Is the given element "zero"?
    fn is_zero(&self) -> bool;

    /// Returns the "one" element of the ring.
    /// This is already a bit scary to implement without referencing the ring
    /// because mod 1 this is 0, but we will just disallow this.
    /// (Technically [`RingElement::zero`] is as well because this is just the
    /// type we represent ring elements with e.g. [`BigUint`] and we could have
    /// `BigUint::zero` is the zero element in one ring but `BigUint::from(42)`
    /// is in another ring.)
    fn one() -> Self;

    /// Is the given element "one"?
    fn is_one(&self) -> bool;
}

macro_rules! impl_ring_element {
    ($t:ident) => {
        impl RingElement for $t {
            fn zero() -> Self {
                Zero::zero()
            }

            fn is_zero(&self) -> bool {
                Zero::is_zero(self)
            }

            fn one() -> Self {
                One::one()
            }

            fn is_one(&self) -> bool {
                One::is_one(self)
            }
        }
    };
}

pub(crate) use impl_ring_element;

/// A ring.
///
/// `Clone` is here for convenience. All current implementations have reasonably
/// efficient clones and deriving `Clone` for structs that have a `R: Ring`
/// generic will not work because the `R` is not `Clone`.
pub trait Ring: 'static + Clone + PartialEq + Eq + std::fmt::Debug {
    /// The type of the elements of the ring.
    type Element: RingElement;

    /// Is the ring a domain, i.e. is 0 the only zero divisor?
    fn is_domain(&self) -> bool;

    /// Returns the "zero" element of the ring.
    fn zero() -> Self::Element {
        Self::Element::zero()
    }

    /// Returns the "one" element of the ring.
    fn one() -> Self::Element {
        Self::Element::one()
    }

    /// Returns the additive inverse of 1.
    fn negative_one(&self) -> Self::Element {
        self.sub(Self::zero(), &Self::one())
    }

    /// Negates the element in place.
    fn neg_assign(&self, e: &mut Self::Element);

    /// Negates the element, i.e. computes `0 - e`.
    fn neg(&self, mut e: Self::Element) -> Self::Element {
        self.neg_assign(&mut e);
        e
    }

    /// Add an element to another element.
    fn add_assign(&self, l: &mut Self::Element, r: &Self::Element);

    /// Add two elements.
    fn add(&self, mut l: Self::Element, r: &Self::Element) -> Self::Element {
        self.add_assign(&mut l, r);
        l
    }

    /// Subtract one element from another.
    fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element);

    /// Subtract one element from another.
    fn sub(&self, mut l: Self::Element, r: &Self::Element) -> Self::Element {
        self.sub_assign(&mut l, r);
        l
    }

    /// Subtract an element from another element and store the result in the
    /// right-hand side.
    ///
    /// This exists as an optimization. The naive way to do this is to use
    /// [`Ring::sub`] but that requires cloning the left-hand side. The slightly
    /// less naive way (which is the default implementation) is to rewrite the
    /// subtraction as `a - b = -b + a` and to then use [`Ring::neg_assign`]
    /// followed by [`Ring::add_assign`] which avoids the allocation. But not
    /// only does this function exist as a convenience function for that, but it
    /// is also possible that rewriting it like that is less optimial than a
    /// direct way, so this function can then be overriden.
    ///
    /// Note that this isn't required for addition, because addition is
    /// commutative.
    fn sub_assign_rhs(&self, l: &Self::Element, r: &mut Self::Element) {
        self.neg_assign(r);
        self.add_assign(r, l);
    }

    /// Subtract an element given by value from an element given by reference.
    /// See [`Ring::sub_assign_rhs`].
    fn sub_rhs(
        &self, l: &Self::Element, mut r: Self::Element
    ) -> Self::Element {
        self.sub_assign_rhs(l, &mut r);
        r
    }

    /// Multiply two elements.
    fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element);

    /// Multiply two elements.
    fn mul(&self, mut l: Self::Element, r: &Self::Element) -> Self::Element {
        self.mul_assign(&mut l, r);
        l
    }

    /// Multiply two elements and add the result to another element.
    ///
    /// This exists as an optimization and has a default implementation that
    /// just allocates a new element for the product and adds that, but this can
    /// sometimes be avoided. E.g. [`num_bigint`] has a function that does this
    /// but it is not currently exposed.
    /// <https://github.com/rust-num/num-bigint/blob/575cea47d21f969e541a7668751d4a82825d02bd/src/biguint/multiplication.rs#L67C37-L67C76>
    ///
    /// TODO: Expose this function and use it.
    fn mul_add_assign(
        &self,
        acc: &mut Self::Element,
        a: &Self::Element,
        b: &Self::Element,
    ) {
        self.add_assign(acc, &self.mul(a.clone(), b))
    }

    /// Multiply two elements and add the result to another element.
    /// See [`Ring::mul_add_assign`].
    fn mul_add(
        &self,
        mut acc: Self::Element,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Self::Element {
        self.mul_add_assign(&mut acc, a, b);
        acc
    }

    /// [`Ring::mul_add_assign`] but with [`Ring::sub`].
    fn mul_sub_assign(
        &self,
        acc: &mut Self::Element,
        a: &Self::Element,
        b: &Self::Element,
    ) {
        self.sub_assign(acc, &self.mul(a.clone(), b));
    }

    /// [`Ring::mul_add`] but with [`Ring::sub`].
    fn mul_sub(
        &self,
        mut acc: Self::Element,
        a: &Self::Element,
        b: &Self::Element,
    ) -> Self::Element {
        self.mul_sub_assign(&mut acc, a, b);
        acc
    }

    /// Increment the element.
    fn inc_assign(&self, e: &mut Self::Element) {
        self.add_assign(e, &Self::one());
    }

    /// Increment the element.
    fn inc(&self, mut e: Self::Element) -> Self::Element {
        self.inc_assign(&mut e);
        e
    }

    /// Decrement the element.
    fn dec_assign(&self, e: &mut Self::Element) {
        self.sub_assign(e, &Self::one());
    }

    /// Decrement the element.
    fn dec(&self, mut e: Self::Element) -> Self::Element {
        self.dec_assign(&mut e);
        e
    }

    /// Square a number. This exists as an optimization opportunity, although
    /// no current ring takes advantage of it.
    fn square(&self, a: Self::Element) -> Self::Element {
        self.mul(a.clone(), &a)
    }

    /// Check if an element is a unit.
    /// This can potentially be faster than trying to compute the inverse, so if
    /// you don't need the inverse, you should call this.
    fn is_unit(&self, e: &Self::Element) -> bool;

    /// Compute the multiplicative inverse of an element if it is a unit.
    fn inverse(&self, e: &Self::Element) -> Option<Self::Element>;

    /// Check if an element is a zero divisor.
    fn is_zero_divisor(&self, e: &Self::Element) -> bool;

    /// Returns a random element. This is used when generating values for
    /// variables that don't have one in a valuation. The distribution is
    /// whatever you want but usually uniformly if that makes sense for the
    /// ring.
    /// This should maybe be moved into its own trait because rings whose
    /// whose elements can't be sensibly sampled uniformly at random, e.g. Z, Q,
    /// will just panic.
    fn random<R: rand::Rng>(&self, rng: &mut R) -> Self::Element;

    /// Converts the `usize` `n` into an element.
    /// The element is the result of adding `1` `n`-times to itself.
    fn element_from_usize(&self, n: usize) -> Self::Element;

    /// Converts the [`BigUint`] `n` into an element.
    /// The element is the result of adding `1` `n`-times to itself.
    fn element_from_biguint(&self, n: &BigUint) -> Self::Element;

    /// Converts the [`BigInt`] `n` into an element.
    ///
    /// - If `n` is non-negative, this is the result of adding `1` `n`-times to
    ///   itself.
    ///
    /// - If `n` is negative, this is the result of adding `-1` `-n`-times to
    ///   itself.
    fn element_from_bigint(&self, n: &BigInt) -> Self::Element {
        match n.sign() {
            Sign::NoSign => Self::zero(),
            Sign::Plus => self.element_from_biguint(n.magnitude()),
            Sign::Minus => self.neg(self.element_from_biguint(n.magnitude())),
        }
    }

    /// Parse an element from the given iterator over characters.
    fn parse_element(
        &self,
        it: &mut std::iter::Peekable<impl Iterator<Item = char>>
    ) -> Option<Self::Element> {
        // TODO: allow hexadecimal literals.

        // Is there any digit at all?
        // If not for this check, this function would return 0 in this case.
        if !it.peek().is_some_and(char::is_ascii_digit) {
            return None;
        }

        // Store the 10 so that we don't have to recreate it in each iteration.
        let ten = self.element_from_usize(10);

        // The accumulator.
        let mut num = Self::zero();

        // While there are more characters..
        while let Some(c) = it.peek() {
            // ..check if the character is a digit and get the value.
            let Some(d) = c.to_digit(10) else {
                // If there are no more digits, we are done.
                break;
            };

            // Multiply by 10.
            self.mul_assign(&mut num, &ten);

            // Add the digit.
            self.add_assign(&mut num, &self.element_from_usize(d as usize));

            // Consume the character.
            it.next();
        }

        Some(num)
    }

    /// Parse an element from the given string.
    fn element_from_string(&self, s: &str) -> Option<Self::Element> {
        let mut it = s.chars().peekable();
        self.parse_element(&mut it)
    }

    /// The name of the data type for the elements of this ring. This is used
    /// when formatting expressions. It doesn't have to make sense.
    fn data_type_name(&self, formatter: Formatter) -> impl std::fmt::Display;
}

/// Basically the integers mod a power of two.
pub trait BinaryRing: IntDivRing {
    /// The number of bits.
    fn bits(&self) -> u32;

    /// Returns the bit at the given position.
    fn bit(e: &Self::Element, i: u32) -> bool;

    /// Bitwise not.
    fn not_assign(&self, e: &mut Self::Element);

    /// Bitwise not.
    fn not(&self, mut e: Self::Element) -> Self::Element {
        self.not_assign(&mut e);
        e
    }

    /// Bitwise and.
    fn and_assign(l: &mut Self::Element, r: &Self::Element);

    /// Bitwise and.
    fn and(mut l: Self::Element, r: &Self::Element) -> Self::Element {
        Self::and_assign(&mut l, r);
        l
    }

    /// Bitwise or.
    fn or_assign(l: &mut Self::Element, r: &Self::Element);

    /// Bitwise or.
    fn or(mut l: Self::Element, r: &Self::Element) -> Self::Element {
        Self::or_assign(&mut l, r);
        l
    }

    /// Bitwise xor.
    fn xor_assign(l: &mut Self::Element, r: &Self::Element);

    /// Bitwise xor.
    fn xor(mut l: Self::Element, r: &Self::Element) -> Self::Element {
        Self::xor_assign(&mut l, r);
        l
    }

    /// Shift the bits to the left.
    fn shl_assign(&self, l: &mut Self::Element, r: u32);

    /// Shift the bits to the left.
    fn shl(&self, mut l: Self::Element, r: u32) -> Self::Element {
        self.shl_assign(&mut l, r);
        l
    }

    /// Is the number even (is bit-0 0)?
    fn is_even(e: &Self::Element) -> bool {
        !Self::is_odd(e)
    }

    /// Is the number odd (is bit-0 1)?
    fn is_odd(e: &Self::Element) -> bool {
        Self::bit(e, 0)
    }

    /// Count the number of 1 bits.
    fn count_ones(e: &Self::Element) -> u32;

    /// Compute the minimum number of bits needed for the element.
    fn min_bits(e: &Self::Element) -> u32;

    /// Convert the element to a [`BigUint`].
    fn to_representative(e: &Self::Element) -> BigUint;

    /// Convert the element to a [`usize`].
    fn to_usize(e: &Self::Element) -> usize;
}

/// A field.
pub trait Field: Ring {
    /// Divide an element by another element.
    fn div_assign(&self, l: &mut Self::Element, r: &Self::Element);

    /// Divide an element by another element.
    fn div(&self, mut l: Self::Element, r: &Self::Element) -> Self::Element {
        self.div_assign(&mut l, r);
        l
    }
}

/// A ring where you can take approximate square roots of positive (whatever
/// that means) elements. This is used for vectors that should have an (L2)
/// norm. This is not really well defined but makes enough sense for all the
/// rings it is implemented for to be useful.
pub trait SqrtRing: Ring {
    /// Compute an approximate square root.
    /// This takes the element by reference and always allocates a new one
    /// because I don't think there is any type that can compute the square
    /// root in-place.
    /// It also doesn't take a `self` reference, because all the rings we care
    /// about don't need it.
    fn sqrt(e: &Self::Element) -> Self::Element;
}

/// A ring where you can do integer division.
///
/// This is needed by the [`crate::solver`] algorithms.
pub trait IntDivRing: OrderedRing {
    /// Divide two elements and round it to the nearest ring element.
    fn rounded_div(l: &Self::Element, r: &Self::Element) -> Self::Element;

    /// Divide two elements such that the remainder is non-negative.
    fn euclidean_div(l: &Self::Element, r: &Self::Element) -> Self::Element;

    /// Compute the remainder of the euclidean division.
    fn euclidean_rem(l: &Self::Element, r: &Self::Element) -> Self::Element;
}

// For fields this is trivial because we can always divide by non-zero elements
// with zero remainder.
macro_rules! impl_int_div_for_field {
    ($field:ident) => {
        impl IntDivRing for $field {
            fn rounded_div(
                l: &Self::Element,
                r: &Self::Element,
            ) -> Self::Element {
                l / r
            }

            fn euclidean_div(
                l: &Self::Element,
                r: &Self::Element,
            ) -> Self::Element {
                l / r
            }

            fn euclidean_rem(
                _l: &Self::Element,
                _r: &Self::Element,
            ) -> Self::Element {
                Self::zero()
            }
        }
    };
}

pub(crate) use impl_int_div_for_field;

use crate::formatter::Formatter;

/// A ring where you can compare two elements.
///
/// This plays two roles:
/// - In [`crate::solver`] many of the functions need [`OrderedRing::cmp_abs`]
///   and some need [`OrderedRing::is_negative`].
///
/// - In many lattice algorithms we need to compare two elements of a
///   [`crate::lattice::WorkingType`], e.g. in [`crate::lattice::lll`].
///   And we also need [`OrderedRing::cmp_abs`] for solving systems of linear
///   equations over a field there. For this reason this is also implemented
///   for [`super::F32`] and [`super::F64`] which rust doesn't consider to be
///   ordered ([`Ord`]) because of NaNs. The implementation for them just
///   `unwrap`s the [`PartialOrd::partial_cmp`].
pub trait OrderedRing: Ring {
    /// Compare two elements.
    fn cmp(&self, l: &Self::Element, r: &Self::Element) -> std::cmp::Ordering;

    /// Is `l` less than `r`?
    fn is_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp(l, r).is_lt()
    }

    /// Is `l` less than or equal to `r`?
    fn is_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp(l, r).is_le()
    }

    /// Is `l` greater than or equal to `r`?
    fn is_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp(l, r).is_ge()
    }

    /// Is `l` greater than `r`?
    fn is_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp(l, r).is_gt()
    }

    /// Is the given element greater than 0?
    fn is_positive(&self, e: &Self::Element) -> bool {
        self.is_gt(e, &Self::zero())
    }

    /// Is the given element less than 0?
    fn is_negative(&self, e: &Self::Element) -> bool {
        self.is_lt(e, &Self::zero())
    }

    /// Compute the absolute value of the element, which is defined here as
    /// `is_negative(e) ? neg(e) : e`.
    fn abs_assign(&self, e: &mut Self::Element) {
        if self.is_negative(e) {
            self.neg_assign(e);
        }
    }

    /// Compute the absolute value of the element.
    /// See [`OrderedRing::abs_assign`].
    fn abs(&self, mut e: Self::Element) -> Self::Element {
        self.abs_assign(&mut e);
        e
    }

    /// Compare the absolute value of two elements.
    fn cmp_abs(
        &self,
        l: &Self::Element,
        r: &Self::Element,
    ) -> std::cmp::Ordering {
        self.cmp(&self.abs(l.clone()), &self.abs(r.clone()))
    }

    /// Is the absolute value of `l` equal to the absolute value of `r`?
    fn is_abs_eq(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_eq()
    }

    /// Is the absolute value of `l` not equal to the absolute value of `r`?
    fn is_abs_ne(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_ne()
    }

    /// Is the absolute value of `l` less than the absolute value of `r`?
    fn is_abs_lt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_lt()
    }

    /// Is the absolute value of `l` less than or equal to the absolute value of
    /// `r`?
    fn is_abs_le(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_le()
    }

    /// Is the absolute value of `l` greater than or equal to the absolute value
    /// of `r`?
    fn is_abs_ge(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_ge()
    }

    /// Is the absolute value of `l` greater than the absolute value of `r`?
    fn is_abs_gt(&self, l: &Self::Element, r: &Self::Element) -> bool {
        self.cmp_abs(l, r).is_gt()
    }
}
