//! Polynomials.

use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign};
use num_bigint::BigInt;
use num_traits::{Zero, One};
use crate::{expr::{ExprOp, Expr}, keep_bits, keep_bits_mut};

/// Represents a polynomial with integer coefficients.
///
/// `coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + ...`
///
/// The coefficients should maybe be stored in reverse order.
#[derive(Clone, Debug)]
pub struct Poly {
    pub coeffs: Vec<BigInt>,
}

impl Poly {
    /// Returns the constant zero polynomial.
    pub const fn zero() -> Self {
        Self {
            coeffs: Vec::new()
        }
    }

    /// Returns the constant one polynomial.
    pub fn one() -> Self {
        Self::constant(1.into())
    }

    /// Returns the constant polynomial.
    pub fn constant(a: BigInt) -> Self {
        Self {
            coeffs: vec![a]
        }
    }

    /// Returns a polynomial from a list of coefficients.
    /// This function is not very intuitive because
    /// the coefficients are passed in the reverse order of what
    /// you are used to.
    /// Another reason to change the representation sometime.
    pub fn from_vec<T: Into<BigInt>>(v: Vec<T>) -> Self {
        Self {
            coeffs: v.into_iter().map(|e| e.into()).collect()
        }
    }

    /// Returns the degree of the polynomial.
    /// The degree of the zero polynomial is zero here!
    pub fn degree(&self) -> usize {
        match self.coeffs.len() {
            0 => 0,
            n => n - 1,
        }
    }

    /// Returns the number of coefficient.
    #[allow(clippy::len_without_is_empty)]
    pub fn len(&self) -> usize {
        self.coeffs.len()
    }

    /// Checks whether a truncated polynomial is zero.
    pub fn is_zero(&self) -> bool {
        self.len() == 0
    }

    /// Checks whether a truncated polynomial is the identity polynomial.
    pub fn is_id(&self) -> bool {
        debug_assert!(self.coeffs.last().is_none_or(|i| !i.is_zero()),
            "Truncate the polynomial before checking if it is the identity");
        matches!(&self.coeffs[..], [z, o] if z.is_zero() && o.is_one())
    }

    /// Evaluate the polynomial at a using Horner's method.
    pub fn eval(&self, a: &BigInt) -> BigInt {
        // Iterate over the coefficients in reverse order.
        let mut iter = self.coeffs.iter().rev();

        // The last coefficient is the initial value.
        let mut v = iter.next().map_or_else(BigInt::zero, |c| c.clone());
        for c in iter {
            // Multiply the current value by a and add the next coefficient.
            v *= a;
            v += c;
        }

        v
    }

    /// Evaluate the polynomial mod 2^n.
    pub fn eval_bits(&self, a: &BigInt, n: u32) -> BigInt {
        // Iterate over the coefficients in reverse order.
        let mut iter = self.coeffs.iter().rev();

        // The last coefficient is the initial value.
        let mut v = iter.next()
            .map_or_else(BigInt::zero, |c| keep_bits(c, n));
        for c in iter {
            // Multiply the current value by a and add the next coefficient.
            v *= a;
            v += c;
            keep_bits_mut(&mut v, n);
        }

        v
    }

    /// Truncates leading zero coefficients.
    pub fn truncate(&mut self) {
        for i in (0..self.len()).rev() {
            if !self.coeffs[i].is_zero() {
                break;
            }
            self.coeffs.pop();
        }
    }

    /// Returns the truncated polynomial.
    pub fn truncated(mut self) -> Self {
        self.truncate();
        self
    }

    /// Multiplies the polynomial by a linear factor (x-a).
    pub fn mul_linfac(&mut self, a: &BigInt) {
        // p(x) * (x-a) = p(x) * x - p(x) * a

        // Shift every coefficient to the left
        // which corresponds to a multiplication by x.
        self.coeffs.insert(0, BigInt::zero());

        // Now subtract a times the original polynomial.
        for i in 0..self.coeffs.len()-1 {
            let m = a * &self.coeffs[i+1];
            self.coeffs[i] -= m;
        }
    }

    /// Computes the derivative of the polynomial.
    pub fn derivative(&self) -> Self {
        if self.len() == 0 {
            return Poly::zero();
        }

        let mut coeffs = Vec::with_capacity(self.len()-1);

        for (e, c) in self.coeffs[1..].iter().enumerate() {
            coeffs.push(c * (e + 1) as u64);
        }

        Self { coeffs }
    }

    /// Reduces all coefficients mod 2^n.
    pub fn mod_coeff(&mut self, n: u32) {
        for c in &mut self.coeffs {
            keep_bits_mut(c, n);
        }
    }

    /// Shift the coefficients to the left by m.
    pub fn shl_coeff(&mut self, m: u32) {
        for c in &mut self.coeffs {
            *c <<= m;
        }
    }

    /// Returns an expression that uses
    /// Horner's method to evaluate the polynomial in x.
    pub fn to_expr(&self) -> ExprOp {
        let mut it = self.coeffs.iter().rev();
        let mut e = match it.next() {
            None => ExprOp::zero(),
            Some(c) => ExprOp::Const(c.clone()),
        };

        let x = Expr::new(ExprOp::Var("x".into()));

        for c in it {
            e = ExprOp::Mul(x.clone(), e.into());
            if !c.is_zero() {
                e = ExprOp::Add(ExprOp::Const(c.clone()).into(), e.into());
            }
        }

        e
    }

    /// Parse an expression from a string.
    /// This can either be a "normalized" expression:
    /// ```
    /// use mba::poly::Poly;
    /// // The variable has to be called x or X.
    /// // There may only be one term per monomial, e.g. not `3x + 4x`.
    /// // Spaces are optional.
    /// Poly::parse("3x^2 + 4x + 5").unwrap();
    /// ```
    /// Or a space-separated list of coefficients a_d ... a_0:
    /// ```
    /// use mba::poly::Poly;
    /// Poly::parse("3 4 5").unwrap();
    /// ```
    pub fn parse<T: Into<String>>(str: T) -> Result<Self, String> {
        let mut str = str.into();

        // It should all be ascii.
        if !str.is_ascii() {
            return Err("The string contains non-ascii characters.".to_owned());
        }

        // Make everything lowercase so you can use x or X.
        str.make_ascii_lowercase();

        let mut coeffs = Vec::new();
        let mut str = str.into_bytes();

        if str.contains(&b'x') {
            str.retain(|c| *c != b' ');

            let mut i = 0;
            let mut last_i = usize::MAX;
            while i < str.len() {
                if i == last_i {
                    return Err(format!(
                        "Unexpected input at {i}: {}.", str[i] as char
                    ));
                }
                last_i = i;

                // Parse the sign.
                let sign = match str[i] {
                    b'+' => { i += 1; false },
                    b'-' => { i += 1; true },
                    _ => false,
                };

                // Parse the coefficient.
                let mut c = BigInt::one();

                // Is there a coefficient?
                let is_digit = str.get(i)
                    .ok_or("Unexpected end of input.")?
                    .is_ascii_digit();

                if is_digit {
                    c = Zero::zero();

                    // Parse the number.
                    while str.get(i).is_some_and(u8::is_ascii_digit) {
                        c *= 10;
                        c += str[i] - b'0';
                        i += 1;
                    }

                    // Skip the `*` if it exists.
                    if str.get(i).is_some_and(|c| *c == b'*') {
                        i += 1;
                    }
                }

                if sign {
                    c = -c;
                }

                // Parse the exponent.
                let mut e = 0;

                if str.get(i).is_some_and(|c| *c == b'x') {
                    // Skip past the `x`.
                    i += 1;
                    e = 1;

                    // If there is an exponent, parse it.
                    if str.get(i).is_some_and(|c| *c == b'^') {
                        i += 1;
                        // Is there a number?
                        let is_digit = str.get(i)
                            .ok_or("Unexpected end of input.")?
                            .is_ascii_digit();
                        if !is_digit {
                            return Err(format!(
                                "Failed to parse exponent at {i}: {}.",
                                str[i] as char
                            ));
                        }

                        // Parse the number.
                        e = 0;
                        while str.get(i).is_some_and(u8::is_ascii_digit) {
                            e *= 10;
                            e += (str[i] - b'0') as usize;
                            i += 1;
                        }
                    }
                }

                if e >= coeffs.len() {
                    coeffs.resize(e + 1, Zero::zero());
                }

                coeffs[e] = c;
            }
        } else {
            for c in str.split(|c| *c == b' ') {
                let c = BigInt::parse_bytes(c, 10)
                    .ok_or("Failed to parse coefficient.")?;
                coeffs.push(c);
            }
            coeffs.reverse();
        }

        Ok(Poly { coeffs }.truncated())
    }
}

impl Add for &Poly {
    type Output = Poly;
    fn add(self, rhs: Self) -> Self::Output {
        // Which polynomial is of larger degree?
        let (min, max) = if self.len() >= rhs.len() {
            (rhs, self)
        } else {
            (self, rhs)
        };

        let mut coeffs = Vec::with_capacity(max.len());

        // Add up all coefficients that exist in both.
        self.coeffs.iter()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| coeffs.push(l + r));

        // Push the remaining coefficients.
        for c in &max.coeffs[min.len()..] {
            coeffs.push(c.clone());
        }

        Poly { coeffs }
    }
}

impl AddAssign<&Poly> for Poly {
    fn add_assign(&mut self, rhs: &Poly) {
        // Add the coefficients that exist in both.
        self.coeffs.iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| *l += r);

        // Push the remaining coefficients should rhs have more.
        for c in &rhs.coeffs[self.len()..] {
            self.coeffs.push(c.clone());
        }
    }
}

impl AddAssign<&BigInt> for Poly {
    fn add_assign(&mut self, rhs: &BigInt) {
        match self.coeffs.is_empty() {
            true => self.coeffs.push(rhs.clone()),
            false => self.coeffs[0] += rhs,
        }
    }
}

impl Sub for &Poly {
    type Output = Poly;
    fn sub(self, rhs: Self) -> Self::Output {
        // Subtract the rhs for the coefficients that exist in both.
        let mut coeffs = Vec::with_capacity(self.len());
        self.coeffs.iter()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| coeffs.push(l - r));

        // Push the remaining coefficients or their additive inverses
        // depending on what polynomial has the more coefficients.
        if self.len() >= rhs.len() {
            for c in &self.coeffs[rhs.len()..] {
                coeffs.push(c.clone());
            }
        } else {
            for c in &rhs.coeffs[self.len()..] {
                coeffs.push(-c.clone());
            }
        }

        Poly { coeffs }
    }
}

impl SubAssign<&Poly> for Poly {
    fn sub_assign(&mut self, rhs: &Poly) {
        // Subtract the rhs for the coefficients that exist in both.
        self.coeffs.iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| *l -= r);

        // Push the remaining coefficients of the rhs
        // if it has more coefficients.
        for c in &rhs.coeffs[self.len()..] {
            self.coeffs.push(-c.clone());
        }
    }
}

impl SubAssign<&BigInt> for Poly {
    fn sub_assign(&mut self, rhs: &BigInt) {
        match self.coeffs.is_empty() {
            true => self.coeffs.push(-rhs.clone()),
            false => self.coeffs[0] -= rhs,
        }
    }
}

impl Mul for &Poly {
    type Output = Poly;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut coeffs = vec![BigInt::zero(); self.len() + rhs.len() - 1];
        for (i, c) in rhs.coeffs.iter().enumerate() {
            for (j, d)  in self.coeffs.iter().enumerate() {
                coeffs[i + j] += c * d;
            }
        }

        Poly { coeffs }
    }
}

impl MulAssign<&Poly> for Poly {
    fn mul_assign(&mut self, rhs: &Poly) {
        let r = &*self * rhs;
        *self = r;
    }
}

impl Mul<&BigInt> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &BigInt) -> Self::Output {
        let mut r = self.clone();
        r *= rhs;
        r
    }
}

impl MulAssign<&BigInt> for Poly {
    fn mul_assign(&mut self, rhs: &BigInt) {
        for c in &mut self.coeffs {
            *c *= rhs;
        }
    }
}

impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.coeffs
            .iter()
            .enumerate()
            .rev()
            .filter(|(_, c)| !c.is_zero());

        fn write_term(
            f: &mut std::fmt::Formatter<'_>,
            e: usize,
            c: &BigInt
        ) -> std::fmt::Result {
            if e == 0 {
                return write!(f, "{c}");
            }

            if !c.is_one() {
                write!(f, "{c}")?;
            }

            write!(f, "x")?;

            if e > 1 {
                write!(f, "^{e}")?;
            }

            Ok(())
        }

        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => write_term(f, e, c)?,
        };

        for (e, c) in iter {
            f.write_str(" + ")?;
            write_term(f, e, c)?;
        }

        Ok(())
    }
}