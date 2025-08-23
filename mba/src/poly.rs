//! Polynomials.

use rand::Rng;

use crate::{
    expr::{Expr, ExprOp},
    rings::{BinaryRing, Ring, RingElement},
};

/// Represents a polynomial with integer coefficients.
///
/// `coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + ...`
///
/// The coefficients should maybe be stored in reverse order.
#[derive(Clone, Debug)]
pub struct Poly<R: Ring> {
    pub coeffs: Vec<R::Element>,
}

impl<R: Ring> Poly<R> {
    /// Returns the constant zero polynomial.
    pub const fn zero() -> Self {
        Self { coeffs: Vec::new() }
    }

    /// Returns the constant one polynomial.
    pub fn one() -> Self {
        Self::constant(R::one())
    }

    /// Returns the constant polynomial.
    pub fn constant(c: R::Element) -> Self {
        Self { coeffs: vec![c] }
    }

    /// Returns a polynomial from a list of coefficients.
    /// This function is not very intuitive because
    /// the coefficients are passed in the reverse order of what
    /// you are used to.
    /// Another reason to change the representation sometime.
    pub fn from_vec(v: Vec<R::Element>) -> Self {
        Self { coeffs: v }
    }

    pub fn random<Rand: Rng>(rng: &mut Rand, degree: usize, r: &R) -> Self {
        let p: Vec<_> = (0..=degree).map(|_| r.random(rng)).collect();
        Self { coeffs: p }.truncated()
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
        debug_assert!(
            self.coeffs.last().is_none_or(|i| !i.is_zero()),
            "Truncate the polynomial before checking if it is the identity"
        );
        matches!(&self.coeffs[..], [z, o] if z.is_zero() && o.is_one())
    }

    /// Evaluate the polynomial at a using Horner's method.
    pub fn eval(&self, a: &R::Element, r: &R) -> R::Element {
        // Iterate over the coefficients in reverse order.
        let mut iter = self.coeffs.iter().rev();

        // The last coefficient is the initial value.
        let mut v = iter.next().map_or_else(R::zero, |c| c.clone());
        for c in iter {
            // Multiply the current value by a and add the next coefficient.
            v = r.mul(v, a);
            v = r.add(v, c);
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

    /// Returns a struct that can be used to display the polynomial.
    pub fn display<'a>(&'a self, var: &'a str) -> DisplayPoly<'a, R> {
        DisplayPoly { poly: self, var }
    }

    /// Add two polynomials.
    pub fn add_assign(&mut self, rhs: &Self, ring: &R) {
        // Add the coefficients that exist in both.
        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| ring.add_assign(l, r));

        // Push the remaining coefficients should rhs have more.
        for c in &rhs.coeffs[self.len()..] {
            self.coeffs.push(c.clone());
        }
    }

    /// Add a constant to this polynomial.
    pub fn add_assign_const(&mut self, rhs: &R::Element, ring: &R) {
        match self.coeffs.first_mut() {
            None => self.coeffs.push(rhs.clone()),
            Some(c) => ring.add_assign(c, rhs),
        }
    }

    /// Add two polynomials.
    pub fn add(&self, rhs: &Self, ring: &R) -> Self {
        // Which polynomial is of larger degree?
        let (min, max) = if self.len() >= rhs.len() {
            (rhs, self)
        } else {
            (self, rhs)
        };

        let mut coeffs = Vec::with_capacity(max.len());

        // Add up all coefficients that exist in both.
        self.coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| coeffs.push(ring.add(l.clone(), r)));

        // Push the remaining coefficients.
        for c in &max.coeffs[min.len()..] {
            coeffs.push(c.clone());
        }

        Poly { coeffs }
    }

    /// Subtract one polynomial from another.
    pub fn sub_assign(&mut self, rhs: &Self, ring: &R) {
        // Subtract the rhs for the coefficients that exist in both.
        self.coeffs
            .iter_mut()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| ring.sub_assign(l, r));

        // Push the remaining coefficients of the rhs
        // if it has more coefficients.
        for c in &rhs.coeffs[self.len()..] {
            self.coeffs.push(ring.neg(c.clone()));
        }
    }

    /// Subtract a constant from the polynomial.
    pub fn sub_assign_const(&mut self, rhs: &R::Element, ring: &R) {
        match self.coeffs.first_mut() {
            None => self.coeffs.push(ring.neg(rhs.clone())),
            Some(c) => ring.sub_assign(c, rhs),
        }
    }

    /// Subtract one polynomial from another.
    pub fn sub(&self, rhs: &Self, ring: &R) -> Self {
        // Subtract the rhs for the coefficients that exist in both.
        let mut coeffs = Vec::with_capacity(self.len());
        self.coeffs
            .iter()
            .zip(rhs.coeffs.iter())
            .for_each(|(l, r)| coeffs.push(ring.sub(l.clone(), r)));

        // Push the remaining coefficients or their additive inverses
        // depending on what polynomial has the more coefficients.
        if self.len() >= rhs.len() {
            for c in &self.coeffs[rhs.len()..] {
                coeffs.push(c.clone());
            }
        } else {
            for c in &rhs.coeffs[self.len()..] {
                coeffs.push(ring.neg(c.clone()));
            }
        }

        Poly { coeffs }
    }

    /// Multiply two polynomials.
    pub fn mul(&self, rhs: &Self, ring: &R) -> Self {
        let mut coeffs = vec![R::zero(); self.len() + rhs.len() - 1];
        for (i, c) in rhs.coeffs.iter().enumerate() {
            for (j, d) in self.coeffs.iter().enumerate() {
                ring.mul_add_assign(&mut coeffs[i + j], c, d);
            }
        }

        Poly { coeffs }
    }

    /// Multiply two polynomials.
    /// Note that this is **not** more efficient than [`Self::mul`].
    /// We still allocate space for the result and do not re-use the current
    /// vector.
    pub fn mul_assign(&mut self, rhs: &Self, ring: &R) {
        let r = self.mul(rhs, ring);
        *self = r;
    }

    /// Multiply the polynomial by a constant.
    pub fn mul_assign_const(&mut self, rhs: &R::Element, ring: &R) {
        for c in &mut self.coeffs {
            ring.mul_assign(c, rhs);
        }
    }

    /// Multiply the polynomial by a constant.
    pub fn mul_const(mut self, rhs: &R::Element, ring: &R) -> Self {
        self.mul_assign_const(rhs, ring);
        self
    }

    /// Multiplies the polynomial by a linear factor (x-a).
    pub fn mul_linfac(&mut self, a: &R::Element, r: &R) {
        // p(x) * (x-a) = p(x) * x - p(x) * a

        // Shift every coefficient to the left
        // which corresponds to a multiplication by x.
        self.coeffs.insert(0, R::zero());

        // Now subtract a times the original polynomial.
        for i in 0..self.coeffs.len() - 1 {
            let m = r.mul(a.clone(), &self.coeffs[i + 1]);
            r.sub_assign(&mut self.coeffs[i], &m);
        }
    }

    /// Computes the derivative of the polynomial.
    pub fn derivative(&self, r: &R) -> Self {
        if self.len() == 0 {
            return Poly::zero();
        }

        let mut coeffs = Vec::with_capacity(self.len() - 1);

        for (e, c) in self.coeffs[1..].iter().enumerate() {
            coeffs.push(r.mul(r.element_from_usize(e + 1), c));
        }

        Self { coeffs }
    }

    /// Parse an expression from a string.
    /// This can either be a "normalized" expression:
    /// ```
    /// use mba::{rings::U32, poly::Poly};
    /// // The variable has to be called x or X.
    /// // There may only be one term per monomial, e.g. not `3x + 4x`.
    /// // Spaces are optional.
    /// Poly::parse("3x^2 + 4x + 5".to_owned(), &U32).unwrap();
    /// ```
    /// Or a space-separated list of coefficients a_d ... a_0:
    /// ```
    /// use mba::{rings::U32, poly::Poly};
    /// Poly::parse("3 4 5".to_owned(), &U32).unwrap();
    /// ```
    pub fn parse(mut str: String, r: &R) -> Result<Self, String> {
        // It should all be ascii.
        if !str.is_ascii() {
            return Err("The string contains non-ascii characters.".to_owned());
        }

        // Make everything lowercase so you can use x or X.
        str.make_ascii_lowercase();

        let mut coeffs = Vec::new();
        let mut str = str.into_bytes();

        // Cache the 10.
        let ten = r.element_from_usize(10);

        if str.contains(&b'x') {
            str.retain(|c| *c != b' ');

            let mut i = 0;
            let mut last_i = usize::MAX;
            while i < str.len() {
                if i == last_i {
                    return Err(format!(
                        "Unexpected input at {i}: {}.",
                        str[i] as char
                    ));
                }
                last_i = i;

                // Parse the sign.
                let sign = match str[i] {
                    b'+' => {
                        i += 1;
                        false
                    },
                    b'-' => {
                        i += 1;
                        true
                    },
                    _ => false,
                };

                // Is there a coefficient?
                let is_digit = str
                    .get(i)
                    .ok_or("Unexpected end of input.")?
                    .is_ascii_digit();

                // Parse the coefficient.
                let mut c = if is_digit {
                    let mut c = R::zero();

                    // Parse the number.
                    while str.get(i).is_some_and(u8::is_ascii_digit) {
                        r.mul_assign(&mut c, &ten);
                        r.add_assign(
                            &mut c,
                            &r.element_from_usize((str[i] - b'0') as usize),
                        );
                        i += 1;
                    }

                    // Skip the `*` if it exists.
                    if str.get(i).is_some_and(|c| *c == b'*') {
                        i += 1;
                    }

                    c
                } else {
                    R::one()
                };

                if sign {
                    r.neg_assign(&mut c);
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
                        let is_digit = str
                            .get(i)
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
                    coeffs.resize(e + 1, R::zero());
                }

                coeffs[e] = c;
            }
        } else {
            for c in str.split(|c| *c == b' ') {
                let c = r
                    .parse_element(&mut c.iter().map(|&c| c as char).peekable())
                    .ok_or("Failed to parse coefficient.")?;
                coeffs.push(c);
            }
            coeffs.reverse();
        }

        Ok(Poly { coeffs }.truncated())
    }
}

impl<R: BinaryRing> Poly<R> {
    /// Shift the coefficients to the left by m.
    pub fn shl_coeff(&mut self, m: u32, r: &R) {
        for c in &mut self.coeffs {
            r.shl_assign(c, m);
        }
    }

    /// Returns an expression that uses
    /// Horner's method to evaluate the polynomial in x.
    pub fn to_expr(&self) -> ExprOp<R> {
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
}

pub struct DisplayPoly<'a, R: Ring> {
    poly: &'a Poly<R>,
    var: &'a str,
}

impl<'a, R: Ring> DisplayPoly<'a, R> {
    pub fn new(poly: &'a Poly<R>, var: &'a str) -> Self {
        Self { poly, var }
    }
}

impl<'a, R: Ring> std::fmt::Display for DisplayPoly<'a, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self
            .poly
            .coeffs
            .iter()
            .enumerate()
            .rev()
            .filter(|(_, c)| !c.is_zero());

        fn write_term<E: RingElement>(
            f: &mut std::fmt::Formatter<'_>,
            e: usize,
            c: &E,
            var: &str,
        ) -> std::fmt::Result {
            if e == 0 {
                return write!(f, "{c}");
            }

            if !c.is_one() {
                write!(f, "{c}")?;
            }

            write!(f, "{var}")?;

            if e > 1 {
                write!(f, "^{e}")?;
            }

            Ok(())
        }

        match iter.next() {
            None => write!(f, "0")?,
            Some((e, c)) => write_term(f, e, c, self.var)?,
        };

        for (e, c) in iter {
            f.write_str(" + ")?;
            write_term(f, e, c, self.var)?;
        }

        Ok(())
    }
}

impl<R: Ring> std::fmt::Display for Poly<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display("x").fmt(f)
    }
}
