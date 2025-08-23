use super::{BinaryBigInt, BinaryRing, U8, U16, U32, U64, U128};
use crate::rings::RingElement as _;

/// This is here to allow us to specialize the `mod_inv` method for each type.
/// For `u8` we compute a lookup table at compile time and can just look up the
/// answer. The other algorithms use Newton's method/Hensel's lemma starting
/// from the inverse mod 2^8 that we store in the lookup table.
pub trait ModInv: BinaryRing {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element;
}

impl ModInv for U8 {
    fn mod_inv(&self, &e: &Self::Element) -> Self::Element {
        MOD_INV_TABLE[e as usize]
    }
}

impl ModInv for U16 {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element {
        mod_inv(e, self)
    }
}

impl ModInv for U32 {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element {
        mod_inv(e, self)
    }
}

impl ModInv for U64 {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element {
        mod_inv(e, self)
    }
}

impl ModInv for U128 {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element {
        mod_inv(e, self)
    }
}

impl ModInv for BinaryBigInt {
    fn mod_inv(&self, e: &Self::Element) -> Self::Element {
        mod_inv(e, self)
    }
}

/// Generic implementation of Newton's method for computing the inverse of an
/// element mod a power of two.
fn mod_inv<R: BinaryRing>(e: &R::Element, r: &R) -> R::Element {
    // We initialize the inverse with the inverse of e mod 2^8.
    let mut x =
        r.element_from_usize(MOD_INV_TABLE[R::to_usize(e) & 0xff] as usize);

    // We need the number two in each iteration.
    let two = r.element_from_usize(2);

    let mut i = 0;

    // We need to do 3 less iterations because we start with 8 correct bits.
    while i < ceil_ilog2(r.bits()).saturating_sub(3) {
        let rhs = r.sub_rhs(&two, r.mul(x.clone(), e));
        r.mul_assign(&mut x, &rhs);
        i += 1;

        // At this point `x` is the inverse of `e` mod `2^{2^i}` or in other
        // words `x` has `2^i` correct bits.
    }

    // Here, we have `i = ceil_ilog2(bits)`, so `x` has
    // `2^ceil_ilog2(u8::BITS) >= u8::BITS` correct bits. So it is the inverse.

    debug_assert!(r.mul(x.clone(), e).is_one());

    x
}

/// Computes `ceil(log2(e))`, which is the number of iterations you need for
/// Newton's method.
const fn ceil_ilog2(e: u32) -> u32 {
    u32::BITS - (e - 1).leading_zeros()
}

/// Computes the inverse of `e` mod 2^8. Assumes that `e` is odd.
const fn compute_mod_inv_table_entry(e: u8) -> u8 {
    // We also Newton's method here, although it is a lot less important.

    // The initial value is 1, because 1 is the inverse of 1 mod 2.
    let mut x = 1u8;

    let mut i = 0;
    while i < ceil_ilog2(u8::BITS) {
        x = x.wrapping_mul(2u8.wrapping_sub(x.wrapping_mul(e)));
        i += 1;
    }

    assert!(x.wrapping_mul(e) == 1);

    x
}

const fn mod_inv_table() -> [u8; 256] {
    let mut table = [0; 256];
    let mut i = 1u8;
    loop {
        table[i as usize] = compute_mod_inv_table_entry(i);
        i = i.wrapping_add(2);
        if i == 1 {
            break;
        }
    }
    table
}

const MOD_INV_TABLE: [u8; 256] = mod_inv_table();
