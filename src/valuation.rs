//! Basically a key-value store for variable names and their values,
//! but you can specify what to do when a variable is not found.

use num_bigint::{BigInt, RandBigInt};
use num_traits::Zero;
use crate::Symbol;

/// Stores values that should be substituted into variables.
#[derive(Debug)]
pub struct Valuation {
    /// The key value pairs are stored as a Vector
    /// because I doubt a hashmap/tree would be faster
    /// when there are so few variables.
    vals: Vec<(Symbol, BigInt)>,

    /// A function that is called when the value of a variable
    /// is requested but not found in the valuation.
    missing: MissingValue,
}

impl Valuation {
    /// An empty valuation that will panic when any variable is requested.
    pub fn empty() -> Self {
        Self { vals: Vec::new(), missing: MissingValue::panic() }
    }

    /// A valuation that returns zero for any variable.
    pub fn zero() -> Self {
        Self { vals: Vec::new(), missing: MissingValue::zero() }
    }

    /// A valuation that returns a random value for any variable.
    /// The value will be consistent across multiple uses of the same variable.
    /// It will be stored in the valuation.
    pub fn random(bits: u32) -> Self {
        Self { vals: Vec::new(), missing: MissingValue::random(bits) }
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list, it will panic.
    pub fn from_vec_panic(vals: Vec<(Symbol, BigInt)>) -> Self {
        Self { vals, missing: MissingValue::panic() }
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list, it will return zero.
    pub fn from_vec_zero(vals: Vec<(Symbol, BigInt)>) -> Self {
        Self { vals, missing: MissingValue::zero() }
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list,
    /// it will return a random value.
    /// The value will be consistent across multiple uses of the same variable.
    /// It will be stored in the valuation.
    pub fn from_vec_random(vals: Vec<(Symbol, BigInt)>, bits: u32) -> Self {
        Self { vals, missing: MissingValue::random(bits) }
    }

    /// Returns the value of a variable.
    pub fn value(&mut self, name: Symbol) -> &mut BigInt {
        // Feels like this is a borrow checker limitation,
        // rather than a me problem.
        let vals = unsafe {
            std::mem::transmute::<_, &'static mut Vec<(Symbol, BigInt)>>(
                &mut self.vals
            )
        };

        for (n, v) in vals {
            if *n == name {
                return v;
            }
        }

        // If not, use the missing valuation.
        let new_val = match &mut self.missing {
            MissingValue::Panic => panic!("Variable {} not found in valuation.", name),
            MissingValue::Zero => {
                BigInt::zero()
            },
            MissingValue::Random(bits, state) => {
                state.gen_bigint(*bits as u64)
            },
        };

        self.vals.push((name.to_owned(), new_val));
        &mut self.vals.last_mut().unwrap().1
    }

    /// Sets the value of a variable.
    pub fn set_value(&mut self, name: Symbol, value: BigInt) {
        for (n, v) in &mut self.vals {
            if *n == name {
                *v = value;
                return;
            }
        }

        self.vals.push((name.to_owned(), value));
    }
}

/// What should be done for a variable that is not found in the valuation.
enum MissingValue {
    /// Panic.
    Panic,

    /// Return zero.
    Zero,

    /// Return a random value with the given number of bits.
    /// We use a box here, because the random number generator
    /// can be quite large, and we want the [Valuation] to be small.
    Random(u32, Box<dyn rand::RngCore>),
}

impl MissingValue {
    /// Panic on a missing valuation.
    pub fn panic() -> Self {
        Self::Panic
    }

    /// Return zero on a missing valuation.
    pub fn zero() -> Self {
        Self::Zero
    }

    /// Return a random value with the given number of bits.
    /// The value will be stored in the valuation so
    /// subsequent uses of the variable will have the same value.
    /// This uses the [rand::rngs::StdRng] internally.
    pub fn random(bits: u32) -> Self {
        use rand::{thread_rng, rngs::StdRng, SeedableRng};
        // Create a new random state and seed it.
        let rng = Box::new(StdRng::from_rng(thread_rng()).unwrap());
        Self::Random(bits, rng)
    }
}

impl std::fmt::Debug for MissingValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MissingValue::Panic => write!(f, "Panic"),
            MissingValue::Zero => write!(f, "Zero"),
            MissingValue::Random(bits, _) => write!(f, "Random({})", bits),
        }
    }
}