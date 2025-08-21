//! A trait for the kind of number
//! ([ring](https://en.wikipedia.org/wiki/Ring_(mathematics))) we are working
//! with. These traits exist and look the way they do for multiple reasons:
//!
//! - So we can have extremely fast implementations for common bit widths,
//!   i.e. 8, 16, 32, 64, that use single instructions for wrapping operations
//!   and don't have to go through an insane amount of `BigInt` overhead.
//! - We can have slower but more general implementations using `BigInt`s
//!   that allow the modulus to be set at runtime rather than compile time.
//!
//! The most important trait is [`Ring`] which stores information about the
//! ring we are working in and all operations on the ring are implemented in
//! this trait. E.g. or [`U8`], [`U16`], it is empty, and for [`BigIntModN`]
//! it contains the modulus.
//!
//! An instance of [`Ring`] needs to be passed to any function that uses the
//! ring which seems unnecessary for [`U8`] but for [`BigIntModN`] we don't
//! know how to (e.g.) add elements otherwise. The other option would be to
//! store the modulus with every element of the ring which would waste memory
//! and we'd have to check that the moduli are equal for arithmetic operations.
//!
//! The [`BinaryRing`] trait is a trait that essentially just means
//! that we are in a ring modulo a power of two.
//! This allows us in [`BinaryBigInt`] to just store the number of bits instead
//! of the whole modulus.
//! It has the bitwise operations on it.

mod binary_bigint;
mod float;
mod integers;
mod mod_inv_pow_two;
mod mod_n;
mod primitive_int;
mod primitive_uint;
mod rationals;
mod traits;
mod var_bits_primitive;

pub use binary_bigint::*;
pub use float::*;
pub use integers::*;
pub use mod_n::*;
pub use primitive_int::*;
pub use primitive_uint::*;
pub use rationals::*;
pub use traits::*;
pub use var_bits_primitive::*;

use num_bigint::{BigInt, BigUint, RandBigInt};
use num_integer::Integer as _;
use num_rational::BigRational;
use num_traits::{Euclid, One, Signed, Zero};

impl_ring_element!(BigUint);
impl_ring_element!(BigInt);
impl_ring_element!(BigRational);

/// Choose the most efficient [`BinaryRing`] for the given number of bits.
#[macro_export]
macro_rules! choose_binary_ring {
    ($e:expr, $r:ident = $bits:expr) => {{
        let bits = $bits;
        if bits < 8 {
            let $r = $crate::rings::VarU8::new(bits);
            $e
        } else if bits == 8 {
            let $r = $crate::rings::U8;
            $e
        } else if bits < 16 {
            let $r = $crate::rings::VarU16::new(bits);
            $e
        } else if bits == 16 {
            let $r = $crate::rings::U16;
            $e
        } else if bits < 32 {
            let $r = $crate::rings::VarU32::new(bits);
            $e
        } else if bits == 32 {
            let $r = $crate::rings::U32;
            $e
        } else if bits < 64 {
            let $r = $crate::rings::VarU64::new(bits);
            $e
        } else if bits == 64 {
            let $r = $crate::rings::U64;
            $e
        } else if bits < 128 {
            let $r = $crate::rings::VarU128::new(bits);
            $e
        } else if bits == 128 {
            let $r = $crate::rings::U128;
            $e
        } else {
            let $r = $crate::rings::BinaryBigInt::new(bits);
            $e
        }
    }};
}

/// This should eventually be the fastest way to initialize a [`BigUint`] with a
/// power of two. But I haven't actually measured it so I am just guessing that
/// this is faster than `(BigUint::one() << bit)`.
pub(crate) fn biguint_pow2(bit: u32) -> BigUint {
    let mut i = <BigUint as Zero>::zero();
    i.set_bit(bit as u64, true);
    i
}

/// When the `m` is guaranteed to be in `[0, 2n)` where `n` is the modulus,
/// reduce mod n. I am guessing checking if `m >= n` and subtracting `n`
/// will be a lot faster than `%`.
pub(crate) fn reduce_simple(m: &mut BigUint, n: &BigUint) {
    if &*m >= n {
        *m -= n;
    }
}

/// Negates an element without allocating. There should really be a `NegAssign`
/// trait for this in `num_traits`. The `rug` crate had this. Or at least a
/// function on `BigInt`/`BigUint`/`BigRational`.
pub(crate) fn neg_assign<T: std::ops::Neg<Output = T> + Default>(e: &mut T) {
    *e = -std::mem::take(e);

    // This is unsound if `neg` can panic.
    //unsafe {
    //    std::ptr::write(e, -std::ptr::read(e));
    //}
}

/// Negates the element mod m without allocating.
/// See [`neg_assign`].
pub(crate) fn neg_assign_mod(e: &mut BigUint, m: &BigUint) {
    if Zero::is_zero(e) {
        return;
    }
    *e = m - std::mem::take(e);

    // This is unsound because `sub` can panic if you use the
    // `Ring` functions with invalid representatives.
    //unsafe {
    //    std::ptr::write(e, m - std::ptr::read(e));
    //}
}

/// Bitwise negates the element mod m without allocating.
/// See [`neg_assign`].
///
/// The way this works is by expressing the negation using arithmetic.
/// Usually to subtract e from a, we compute the two's complement of e and add
/// it to a. The two's complement means negating the bits of e and adding one.
/// So, `0 - e = -e = !e + 1` => `!e = -e - 1 = -(e + 1)`.
/// Since we are doing this arithmetic on `BigUint`s, which don't do wrapping
/// arithmetic, we need to be careful about two things:
/// - Since `0 <= e < m`, `-m <= -(e + 1) < 0`, so we can add `m` (subtract this
///   from `m`) to get a value `0 <= m - (e + 1) < m` which is the correct
///   representative.
/// - We also need to be careful to never use negative numbers during the
///   computation, since `BigUint` can only store non-negative integers.
///   If we wrote this as `-e - 1 + m` we would be using negative numbers, but
///   by subtracting the `e + 1` from `m` it doesn't.
pub(crate) fn not_assign_mod(e: &mut BigUint, m: &BigUint) {
    *e = m - (std::mem::take(e) + 1u32);

    // This is unsound because `sub` can panic when you use the
    // `Ring` functions with invalid representatives.
    //unsafe {
    //    std::ptr::write(e, m - (std::ptr::read(e) + 1u32))
    //}
}

#[cfg(test)]
mod test {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    pub fn test_inverse<R: Ring>(r: &R) {
        let rng = &mut StdRng::seed_from_u64(0);

        for _ in 0..1000 {
            let e = r.random(rng);
            let inv = r.inverse(&e);
            let is_unit = r.is_unit(&e);
            if let Some(inv) = inv {
                let prod = r.mul(e.clone(), &inv);
                assert_eq!(
                    prod,
                    R::one(),
                    "`inverse` returned {inv} but {e} * {inv} = {prod}"
                );
                assert!(
                    is_unit,
                    "`inverse` returned {inv} but `is_unit` returned false"
                );
            } else {
                assert!(
                    !is_unit,
                    "`inverse` returned None but `is_unit` returned true"
                );
            }
        }
    }

    pub fn test_element_from_biguint<R>(r: &R)
    where
        R: BinaryRing,
    {
        let mut rng = StdRng::seed_from_u64(0);
        for _ in 0..1000 {
            let n = r.random(&mut rng);
            let biguint = R::to_representative(&n);
            let check = r.element_from_biguint(&biguint);
            assert_eq!(
                n, check,
                "element_from_biguint({biguint}) = {n} but got {check}"
            );
        }
    }

    pub fn test_rounded_div<R>(r: &R)
    where
        R: BinaryRing,
        BigInt: From<R::Element>,
    {
        let mut rng = StdRng::seed_from_u64(0);
        for _ in 0..1000 {
            let a = r.random(&mut rng);
            let b = loop {
                let b = r.random(&mut rng);
                if !b.is_zero() {
                    break b;
                }
            };
            let q = R::rounded_div(&a, &b);
            let check = r.element_from_bigint(
                &BigRational::new(BigInt::from(a.clone()), BigInt::from(b.clone()))
                    .round()
                    .into_raw()
                    .0,
            );
            assert_eq!(q, check, "rounded_div({a}, {b}) = {check} but got {q}");
        }
    }
}
