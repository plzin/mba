use crate::formatter::Formatter;

use super::*;

macro_rules! var_primitive_int {
    ($ring:ident, $uint:ident, $uint_ring:ident) => {
        #[derive(Clone, PartialEq, Eq, Debug)]
        pub struct $ring {
            pub bits: $uint,
            pub mask: $uint,
        }

        impl $ring {
            pub fn new(bits: u32) -> Self {
                assert!((1..=$uint::BITS).contains(&bits));
                Self {
                    bits: bits as $uint,
                    mask: ((1 as $uint) << bits).wrapping_sub(1),
                }
            }
        }

        impl Ring for $ring {
            type Element = $uint;

            fn is_domain(&self) -> bool {
                self.bits == 1
            }

            fn neg_assign(&self, e: &mut Self::Element) {
                *e = e.wrapping_neg() & self.mask;
            }

            fn add_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_add(*r) & self.mask;
            }

            fn sub_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_sub(*r) & self.mask;
            }

            fn mul_assign(&self, l: &mut Self::Element, r: &Self::Element) {
                *l = l.wrapping_mul(*r) & self.mask;
            }

            fn is_unit(&self, e: &Self::Element) -> bool {
                *e != 0
            }

            fn inverse(&self, e: &Self::Element) -> Option<Self::Element> {
                $uint_ring.inverse(e).map(|e| e & self.mask)
            }

            fn is_zero_divisor(&self, e: &Self::Element) -> bool {
                e.is_even()
            }

            fn random<R: rand::Rng>(&self, rng: &mut R) -> Self::Element {
                rng.random::<$uint>() & self.mask
            }

            fn element_from_usize(&self, n: usize) -> Self::Element {
                n as $uint & self.mask
            }

            fn element_from_biguint(&self, n: &BigUint) -> Self::Element {
                $uint_ring.element_from_biguint(n) & self.mask
            }

            fn data_type_name(&self, formatter: Formatter) -> impl std::fmt::Display {
                VarBitsPrimitiveIntDataType {
                    bits: self.bits as u32,
                    formatter,
                }
            }
        }

        impl BinaryRing for $ring {
            fn bits(&self) -> u32 {
                self.bits as u32
            }

            fn bit(e: &Self::Element, i: u32) -> bool {
                e & (1 << i) != 0
            }

            fn not_assign(&self, e: &mut Self::Element) {
                *e = !*e & self.mask;
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
                *l &= self.mask;
            }

            fn count_ones(e: &Self::Element) -> u32 {
                $uint_ring::count_ones(e)
            }

            fn min_bits(e: &Self::Element) -> u32 {
                $uint_ring::min_bits(e)
            }

            fn to_representative(e: &Self::Element) -> BigUint {
                $uint_ring::to_representative(e)
            }

            fn to_usize(e: &Self::Element) -> usize {
                $uint_ring::to_usize(e)
            }
        }

        impl IntDivRing for $ring {
            fn rounded_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
                $uint_ring::rounded_div(l, r)
            }

            fn euclidean_div(l: &Self::Element, r: &Self::Element) -> Self::Element {
                $uint_ring::euclidean_div(l, r)
            }

            fn euclidean_rem(l: &Self::Element, r: &Self::Element) -> Self::Element {
                $uint_ring::euclidean_rem(l, r)
            }
        }

        impl_ordered_ring_uint!($ring);
    };
}

pub struct VarBitsPrimitiveIntDataType {
    bits: u32,
    formatter: Formatter,
}

impl std::fmt::Display for VarBitsPrimitiveIntDataType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.formatter {
            Formatter::C => write!(f, "uint{}_t", self.bits),
            Formatter::Rust => write!(f, "Wrapping<u{}>", self.bits),
            Formatter::Tex => write!(f, "uint{}", self.bits),
            Formatter::LLVM => write!(f, "i{}", self.bits),
        }
    }
}

var_primitive_int!(VarU8, u8, U8);
var_primitive_int!(VarU16, u16, U16);
var_primitive_int!(VarU32, u32, U32);
var_primitive_int!(VarU64, u64, U64);
var_primitive_int!(VarU128, u128, U128);
