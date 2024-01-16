//! Polynomials.

use std::ops::{Add, AddAssign, Sub, SubAssign, Mul, MulAssign};

use rug::{Integer, Complete};

use crate::expr::{ExprOp, Expr};

/// Represents a polynomial with integer coefficients.
///
/// `coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + ...`
///
/// The coefficients should maybe be stored in reverse order.
#[derive(Clone, Debug)]
pub struct Poly {
    pub coeffs: Vec<Integer>,
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
    pub fn constant(a: Integer) -> Self {
        Self {
            coeffs: vec![a]
        }
    }

    /// Returns a polynomial from a list of coefficients.
    /// This function is not very intuitive because
    /// the coefficients are passed in the reverse order of what
    /// you are used to.
    /// Another reason to change the representation sometime.
    pub fn from_vec<T: Into<Integer>>(v: Vec<T>) -> Self {
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
        debug_assert!(self.coeffs.last().map_or(true, |i| *i != 0),
            "Truncate the polynomial before checking if it is the identity");
        self.coeffs == [0, 1]
    }

    /// Evaluate the polynomial at a using Horner's method.
    pub fn eval(&self, a: &Integer) -> Integer {
        // Iterate over the coefficients in reverse order.
        let mut iter = self.coeffs.iter().rev();

        // The last coefficient is the initial value.
        let mut v = iter.next().map_or_else(Integer::new, |c| c.clone());
        for c in iter {
            // Multiply the current value by a and add the next coefficient.
            v *= a;
            v += c;
        }

        v
    }

    /// Evaluate the polynomial mod 2^n.
    pub fn eval_bits(&self, a: &Integer, n: u32) -> Integer {
        // Iterate over the coefficients in reverse order.
        let mut iter = self.coeffs.iter().rev();

        // The last coefficient is the initial value.
        let mut v = iter.next()
            .map_or_else(Integer::new, |c| c.clone().keep_bits(n));
        for c in iter {
            // Multiply the current value by a and add the next coefficient.
            v *= a;
            v.keep_bits_mut(n);
            v += c;
            v.keep_bits_mut(n);
        }

        v
    }

    /// Truncates leading zero coefficients.
    pub fn truncate(&mut self) {
        for i in (0..self.len()).rev() {
            if self.coeffs[i] == 0 {
                self.coeffs.pop();
            } else {
                break;
            }
        }
    }

    /// Returns the truncated polynomial.
    pub fn truncated(mut self) -> Self {
        self.truncate();
        self
    }

    /// Multiplies the polynomial by a linear factor (x-a).
    pub fn mul_linfac(&mut self, a: &Integer) {
        // p(x) * (x-a) = p(x) * x - p(x) * a

        // Shift every coefficient to the left
        // which corresponds to a multiplication by x.
        self.coeffs.insert(0, Integer::new());

        // Now subtract a times the original polynomial.
        for i in 0..self.coeffs.len()-1 {
            let m = (a * &self.coeffs[i+1]).complete();
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
            coeffs.push((c * (e + 1) as u64).complete());
        }

        Self { coeffs }
    }

    /// Reduces all coefficients mod 2^n.
    pub fn mod_coeff(&mut self, n: u32) {
        for c in &mut self.coeffs {
            c.keep_bits_mut(n);
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
            if c != &Integer::ZERO {
                e = ExprOp::Add(ExprOp::Const(c.clone()).into(), e.into());
            }
        }

        e
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
            .for_each(|(l, r)| coeffs.push((l + r).complete()));

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

impl AddAssign<&Integer> for Poly {
    fn add_assign(&mut self, rhs: &Integer) {
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
            .for_each(|(l, r)| coeffs.push((l - r).complete()));

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

impl SubAssign<&Integer> for Poly {
    fn sub_assign(&mut self, rhs: &Integer) {
        match self.coeffs.is_empty() {
            true => self.coeffs.push(-rhs.clone()),
            false => self.coeffs[0] -= rhs,
        }
    }
}

impl Mul for &Poly {
    type Output = Poly;
    fn mul(self, rhs: Self) -> Self::Output {
        let mut coeffs = vec![Integer::new(); self.len() + rhs.len() - 1];
        for (i, c) in rhs.coeffs.iter().enumerate() {
            for (j, d)  in self.coeffs.iter().enumerate() {
                coeffs[i + j] += (c * d).complete();
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

impl Mul<&Integer> for &Poly {
    type Output = Poly;
    fn mul(self, rhs: &Integer) -> Self::Output {
        let mut r = self.clone();
        r *= rhs;
        r
    }
}

impl MulAssign<&Integer> for Poly {
    fn mul_assign(&mut self, rhs: &Integer) {
        for c in &mut self.coeffs {
            *c *= rhs;
        }
    }
}

impl std::fmt::Display for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.coeffs.iter().enumerate().rev();

        let mut has_terms = false;
        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => if e == 0 {
                write!(f, "{}", c)?
            } else if *c != 0 {
                write!(f, "{}x^{}", c, e)?;
                has_terms = true;
            },
        };

        for (e, c) in iter {
            if e == 0 {
                if *c != 0 || !has_terms {
                    write!(f, " + {}", c)?;
                }
            } else if *c != 0 {
                write!(f, " + {}x^{}", c, e)?;
                has_terms = true;
            }
        }

        Ok(())
    }
}

impl std::fmt::LowerHex for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.coeffs.iter().enumerate().rev();

        let mut has_terms = false;
        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => if e == 0 {
                write!(f, "{:x}", c)?
            } else if *c != 0 {
                write!(f, "{:x}x^{}", c, e)?;
                has_terms = true;
            },
        };

        for (e, c) in iter {
            if e == 0 {
                if *c != 0 || !has_terms {
                    write!(f, " + {:x}", c)?;
                }
            } else if *c != 0 {
                write!(f, " + {:x}x^{}", c, e)?;
                has_terms = true;
            }
        }

        Ok(())
    }
}

impl std::fmt::UpperHex for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.coeffs.iter().enumerate().rev();

        let mut has_terms = false;
        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => if e == 0 {
                write!(f, "{:X}", c)?
            } else if *c != 0 {
                write!(f, "{:X}x^{}", c, e)?;
                has_terms = true;
            },
        };

        for (e, c) in iter {
            if e == 0 {
                if *c != 0 || !has_terms {
                    write!(f, " + {:X}", c)?;
                }
            } else if *c != 0 {
                write!(f, " + {:X}x^{}", c, e)?;
                has_terms = true;
            }
        }

        Ok(())
    }
}

impl std::fmt::Binary for Poly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.coeffs.iter().enumerate().rev();

        let mut has_terms = false;
        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => if e == 0 {
                write!(f, "{:b}", c)?
            } else if *c != 0 {
                write!(f, "{:b}x^{}", c, e)?;
                has_terms = true;
            },
        };

        for (e, c) in iter {
            if e == 0 {
                if *c != 0 || !has_terms {
                    write!(f, " + {:b}", c)?;
                }
            } else if *c != 0 {
                write!(f, " + {:b}x^{}", c, e)?;
                has_terms = true;
            }
        }

        Ok(())
    }
}
