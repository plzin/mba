//! Basically a key-value store for variable names and their values,
//! but you can specify what to do when a variable is not found.

use rand::{Rng, SeedableRng, rngs::StdRng};

use crate::{
    Symbol,
    rings::{Ring, RingElement as _},
};

/// Stores values that should be substituted into variables.
#[derive(Debug)]
pub struct Valuation<R: Ring> {
    /// The key value pairs are stored as a Vector
    /// because I doubt a hashmap/tree would be faster
    /// when there are so few variables.
    vals: Vec<(Symbol, R::Element)>,

    /// A function that is called when the value of a variable
    /// is requested but not found in the valuation.
    missing: MissingValue,
}

impl<R: Ring> Valuation<R> {
    /// An empty valuation that will panic when any variable is requested.
    pub fn empty() -> Self {
        Self::from_vec_panic(Vec::new())
    }

    /// A valuation that returns zero for any variable.
    pub fn zero() -> Self {
        Self::from_vec_zero(Vec::new())
    }

    /// A valuation that returns a random value for any variable.
    /// The value will be consistent across multiple uses of the same variable.
    /// It will be stored in the valuation.
    pub fn random_seeded(seed: u64) -> Self {
        Self::from_vec_random_seeded(Vec::new(), seed)
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list, it will panic.
    pub fn from_vec_panic(vals: Vec<(Symbol, R::Element)>) -> Self {
        Self { vals, missing: MissingValue::panic() }
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list, it will return zero.
    pub fn from_vec_zero(vals: Vec<(Symbol, R::Element)>) -> Self {
        Self { vals, missing: MissingValue::zero() }
    }

    /// Initializes a valuation from a list of pairs of variables and values.
    /// If a variable is requested that is not in the list,
    /// it will return a random value.
    /// The value will be consistent across multiple uses of the same variable.
    /// It will be stored in the valuation.
    pub fn from_vec_random_seeded(
        vals: Vec<(Symbol, R::Element)>,
        seed: u64,
    ) -> Self {
        let rng = Box::new(StdRng::seed_from_u64(seed));
        Self { vals, missing: MissingValue::random(rng) }
    }

    /// Returns a valuation that is zero for all the given variables and panics
    /// for all other variables.
    ///
    /// It is equivalent to calling [`Valuation::from_vec_panic`] with a list
    /// of all the variables with the value zero.
    pub fn vars_zero_or_panic(vars: &[Symbol]) -> Self {
        Self::from_vec_panic(vars.iter().map(|&v| (v, R::zero())).collect())
    }

    /// Returns the value of a variable.
    pub fn value(&mut self, name: Symbol, r: &R) -> &mut R::Element {
        // Feels like this is a borrow checker limitation,
        // rather than a me problem.
        let vals = unsafe {
            std::mem::transmute::<
                &mut Vec<(Symbol, R::Element)>,
                &'static mut Vec<(Symbol, R::Element)>,
            >(&mut self.vals)
        };

        for (n, v) in vals {
            if *n == name {
                return v;
            }
        }

        // If not, use the missing valuation.
        let new_val = match &mut self.missing {
            MissingValue::Panic => {
                panic!("Variable {name} not found in valuation.")
            },
            MissingValue::Zero => R::zero(),
            MissingValue::Random(rng) => r.random(&mut *rng),
        };

        self.vals.push((name.to_owned(), new_val));
        &mut self.vals.last_mut().unwrap().1
    }

    /// Sets the value of a variable.
    pub fn set_value(&mut self, name: Symbol, value: R::Element) {
        for (n, v) in &mut self.vals {
            if *n == name {
                *v = value;
                return;
            }
        }

        self.vals.push((name.to_owned(), value));
    }

    /// Returns the values of all the seen variables.
    pub fn values(&self) -> &[(Symbol, R::Element)] {
        &self.vals
    }

    /// Increments the valuation in a binary way.
    ///
    /// It is analogous to incrementing `n` bit integer where `n` is the number
    /// of variables.
    ///
    /// This is useful for iterating over all possible values of the variables.
    ///
    /// Be aware of three things:
    /// 1. It sets the variables to 0/-1, not 0/1.
    /// 2. The first variable is the analog of the least significant bit, which
    ///    is the reverse of how we usually write numbers.
    /// 3. It does this only for the variables that have a set value.
    ///
    /// Returns `true` if the valuation has wrapped around, `false` otherwise.
    ///
    /// # Example
    /// ```
    /// use mba::valuation::Valuation;
    /// use mba::rings::U8;
    /// use mba::Symbol;
    ///
    /// let x = Symbol::new("x");
    /// let y = Symbol::new("y");
    ///
    /// let mut v = Valuation::vars_zero_or_panic(&[x, y]);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 0));
    /// assert_eq!(v.inc_bin(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (255, 0));
    /// assert_eq!(v.inc_bin(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 255));
    /// assert_eq!(v.inc_bin(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (255, 255));
    /// assert_eq!(v.inc_bin(&U8), true);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 0));
    /// ```
    pub fn inc_bin(&mut self, r: &R) -> bool {
        for (_, v) in &mut self.vals {
            if v.is_zero() {
                *v = r.negative_one();
                return false;
            } else {
                *v = R::zero();
            }
        }

        true
    }

    /// Similar to [`Valuation::inc_bin`], but the variables iterate over their
    /// whole range.
    ///
    /// # Example
    /// ```
    /// use mba::valuation::Valuation;
    /// use mba::rings::U8;
    /// use mba::Symbol;
    ///
    /// let x = Symbol::new("x");
    /// let y = Symbol::new("y");
    ///
    /// let mut v = Valuation::vars_zero_or_panic(&[x, y]);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 0));
    /// assert_eq!(v.inc_full(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (1, 0));
    /// assert_eq!(v.inc_full(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (2, 0));
    /// for _ in 0..253 {
    ///     assert_eq!(v.inc_full(&U8), false);
    /// }
    /// assert_eq!((*v.value(x), *v.value(y)), (255, 0));
    /// assert_eq!(v.inc_full(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 1));
    /// // Simulating many more increments.
    /// v.set_value(x, 254);
    /// v.set_value(y, 255);
    /// assert_eq!(v.inc_full(&U8), false);
    /// assert_eq!((*v.value(x), *v.value(y)), (255, 255));
    /// assert_eq!(v.inc_full(&U8), true);
    /// assert_eq!((*v.value(x), *v.value(y)), (0, 0));
    /// ```
    pub fn inc_full(&mut self, r: &R) -> bool {
        for (_, v) in &mut self.vals {
            r.inc_assign(v);
            if !v.is_zero() {
                return false;
            }
        }

        true
    }

    /// Updates the values of the variables to random values.
    pub fn update_random<Rand: Rng>(&mut self, rng: &mut Rand, r: &R) {
        for (_, v) in &mut self.vals {
            *v = r.random(rng);
        }
    }
}

/// What should be done for a variable that is not found in the valuation.
pub(crate) enum MissingValue {
    /// Panic.
    Panic,

    /// Return zero.
    Zero,

    /// Return a random value with the given number of bits.
    /// We use a box here, because the random number generator
    /// can be quite large, and we want the [`Valuation`] to be small.
    Random(Box<StdRng>),
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

    /// Return a random value.
    /// The value will be stored in the valuation so subsequent uses of the
    /// variable will have the same value.
    pub fn random(rng: Box<StdRng>) -> Self {
        Self::Random(rng)
    }
}

impl std::fmt::Debug for MissingValue {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            MissingValue::Panic => write!(f, "Panic"),
            MissingValue::Zero => write!(f, "Zero"),
            MissingValue::Random(_) => write!(f, "Random"),
        }
    }
}
