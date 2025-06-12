//! This module contains bitwise expressions ([`BExpr`]), which are just bitwise expressions
//! on n bits, and linear combinations of bitwise expressions ([`LBExpr`]).

use std::{collections::BTreeSet};
use crate::rings::{Ring, RingElement as _, BinaryRing};
use crate::valuation::Valuation;
use crate::{Symbol, ExprOp};

/// LBExpr is short for "Linear combination of Bitwise Expressions"
/// These are the expressions for which rewrite rules can be efficiently
/// generated.
#[derive(Clone, Debug)]
pub struct LBExpr<R: Ring>(pub Vec<(R::Element, BExpr)>);

impl<R: Ring> LBExpr<R> {
    /// Create the empty linear combination.
    /// It evaluates to 0.
    pub fn zero() -> Self {
        Self(Vec::new())
    }

    /// Creates an expression that equals a constant.
    pub fn constant(c: R::Element, r: &R) -> Self {
        Self(vec![(r.neg(c), BExpr::Ones)])
    }

    /// Creates an expression that equals a variable.
    pub fn var<T: Into<Symbol>>(name: T) -> Self {
        Self(vec![(R::one(), BExpr::Var(name.into()))])
    }

    /// Removes all terms with coefficient 0.
    pub fn remove_zero_terms(&mut self) {
        self.0.retain(|(i, _)| !i.is_zero());
    }

    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self) -> Vec<Symbol> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<Symbol>) {
        for (_, e) in &self.0 {
            e.vars_impl(v);
        }
    }

    /// Evaluate an expression with a valuation for the occurring variables.
    pub fn eval(&self, v: &mut Valuation<R>, r: &R) -> R::Element
    where
        R: BinaryRing
    {
        self.0.iter()
            .map(|(i, e)| r.mul(e.eval(v, r), i))
            .fold(R::zero(), |acc, x| r.add(acc, &x))
    }

    /// Returns a measure of the complexity of the expression.
    pub fn complexity(&self) -> u32
    where
        R: BinaryRing
    {
        // Complexity of a coefficient.
        let coeff_complexity = |i: &R::Element| {
            R::count_ones(i) / 2 + R::min_bits(i)
        };

        self.0.iter()
            .filter(|(c, _)| !c.is_zero())
            .map(|(c, e)| coeff_complexity(c) + e.complexity())
            .sum()
    }

    /// Parse a string to an expression.
    /// Note that this function is extremely limited
    /// and expects very specific syntax.
    /// It is used for convenience when testing things and
    /// not really meant to be used by something outside this crate.
    pub fn from_string(mut s: String, r: &R) -> Result<Self, String> {
        s.retain(|c| !matches!(c, ' ' | '\t' | '\n'));
        let mut it = s.chars().peekable();

        // This stores the current linear combination.
        let mut v = Vec::new();

        let mut neg = false;

        // Loop over the string/the summands.
        loop {
            // Are there still characters left?
            // If not then we're done.
            let mut c = match it.peek() {
                None => return Ok(Self(v)),
                Some(c) => *c,
            };

            if c == '-' {
                neg = !neg;
                it.next();
                c = *it.peek().ok_or("Unexpected end of input")?;
            }

            // If this is a digit then we expect num*BExpr.
            if c.is_ascii_digit() {
                // Parse the number.
                let mut num = r.parse_element(&mut it)
                    .ok_or("Expected number")?;

                // If the number is negative then negate it.
                if neg {
                    num = r.neg(num);
                }

                // Is it the expected '*'?
                let lookahead = it.peek();
                match lookahead {
                    Some('*') => {
                        it.next();

                        // Parse the BExpr.
                        let e = BExpr::parse(&mut it, 0)?;

                        // Push it.
                        v.push((num, e));
                    },

                    // If this is a different character then we push -num*(-1).
                    _ => v.push((r.neg(num), BExpr::Ones)),
                }
            } else {
                // We don't have a factor so just parse the BExpr.
                let e = BExpr::parse(&mut it, 0)?;

                let sign = match neg {
                    false => R::one(),
                    true => r.negative_one(),
                };

                // Push sign*e.
                v.push((sign, e));
            }

            // If the next character is not a plus or - then we are done.
            match it.peek() {
                Some('+') => { neg = false; it.next() },
                Some('-') => { neg = true; it.next() },
                _ => return Ok(Self(v)),
            };
        }
    }

    /// Converts an [`LBExpr`] to an [`ExprOp`].
    pub fn to_expr(&self, r: &R) -> ExprOp<R> {
        let mut it = self.0.iter();
        let s = match it.next() {
            None => return ExprOp::Const(R::zero()),
            Some(s) => s,
        };

        let from_summand = |(i, e): &(R::Element, BExpr)| {
            if let BExpr::Ones = e {
                ExprOp::Const(r.neg(i.clone()))
            } else {
                ExprOp::Mul(ExprOp::Const(i.clone()).into(), e.to_expr(r).into())
            }
        };

        let mut cur = from_summand(s);
        for s in it {
            cur = ExprOp::Add(cur.into(), from_summand(s).into());
        }

        cur
    }
}

impl<R: Ring> From<BExpr> for LBExpr<R> {
    fn from(e: BExpr) -> Self {
        Self(vec![(R::one(), e)])
    }
}

impl<R: Ring> std::fmt::Display for LBExpr<R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(crate::formatter::Formatter::C).fmt(f)
    }
}

/// Represents an expression whose value bit `i` depends only on the values of
/// bits `i` of the variables.
///
/// Note that the variant 'Ones' does not equal 1, but a 1 in every bit,
/// which is -1 in two's complement.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum BExpr {
    Ones,
    Var(Symbol),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Xor(Box<Self>, Box<Self>),
    Not(Box<Self>),
}

impl BExpr {
    pub fn var<T: Into<Symbol>>(name: T) -> Self {
        Self::Var(name.into())
    }

    pub fn and(l: Self, r: Self) -> Self {
        Self::And(l.into(), r.into())
    }

    pub fn or(l: Self, r: Self) -> Self {
        Self::Or(l.into(), r.into())
    }

    pub fn xor(l: Self, r: Self) -> Self {
        Self::Xor(l.into(), r.into())
    }

    #[allow(clippy::should_implement_trait)]
    pub fn not(e: Self) -> Self {
        Self::Not(e.into())
    }

    /// Is the topmost operator a unary operator?
    /// If there is no operator (i.e. it is [`BExpr::Ones`] or [`BExpr::Var`])
    /// then this return `true`.
    pub fn is_topmost_unary(&self) -> bool {
        matches!(self, Self::Ones | Self::Var(_) | Self::Not(_))
    }

    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self) -> Vec<Symbol> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<Symbol>) {
        use BExpr::*;
        match self {
            Ones => {},
            Var(c) => drop(v.insert(*c)),
            And(e1, e2) => { e1.vars_impl(v); e2.vars_impl(v) },
            Or(e1, e2) => { e1.vars_impl(v); e2.vars_impl(v) },
            Xor(e1, e2) => { e1.vars_impl(v); e2.vars_impl(v) },
            Not(e) => e.vars_impl(v),
        }
    }

    /// Evaluate an expression with a valuation for the occurring variables.
    pub fn eval<R: BinaryRing>(&self, v: &mut Valuation<R>, r: &R) -> R::Element {
        use BExpr::*;
        match self {
            Ones => r.negative_one(),
            Var(c) => v.value(*c, r).clone(),
            And(e1, e2) => R::and(e1.eval(v, r), &e2.eval(v, r)),
            Or(e1, e2) => R::or(e1.eval(v, r), &e2.eval(v, r)),
            Xor(e1, e2) => R::xor(e1.eval(v, r), &e2.eval(v, r)),
            Not(e) => r.not(e.eval(v, r)),
        }
    }

    /// Rename a variable.
    pub fn rename_var(&mut self, old: Symbol, new: Symbol) {
        use BExpr::*;
        match self {
            Ones => (),
            Var(v) => if *v == old { *v = new },
            And(l, r) => { l.rename_var(old, new); r.rename_var(old, new) },
            Or(l, r) => { l.rename_var(old, new); r.rename_var(old, new) },
            Xor(l, r) => { l.rename_var(old, new); r.rename_var(old, new) },
            Not(e) => e.rename_var(old, new),
        }
    }

    /// Returns some sort of complexity measure of the expression.
    pub fn complexity(&self) -> u32 {
        use BExpr::*;
        match self {
            Ones => 1,
            Var(_) => 1,
            And(l, r) => l.complexity() + r.complexity() + 1,
            Or(l, r) => l.complexity() + r.complexity() + 1,
            Xor(l, r) => l.complexity() + r.complexity() + 1,
            Not(e) => e.complexity() + 1,
        }
    }

    /// Parse a string to an expression.
    pub fn from_string(mut s: String) -> Result<Self, String> {
        s.retain(|c| !matches!(c, ' ' | '\t' | '\n'));
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    pub(self) fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize
    ) -> Result<Self, String> {
        use BExpr::*;

        let c = *it.peek().ok_or("Unexpected end of input")?;

        let mut e = if c == '(' {
            it.next();
            let e = Self::parse(it, 0)?;
            match it.next() {
                Some(')') => e,
                _ => return Err("Expected closing parenthesis".to_string()),
            }
        } else if c == '~' || c == '!' {
            it.next();
            let e = Self::parse(it, 15)?;
            Not(Box::new(e))
        } else if c.is_alphabetic() {
            it.next();
            let mut var = String::from(c);
            loop {
                let Some(c) = it.peek() else {
                    break
                };

                if !c.is_alphanumeric() {
                    break
                }

                var.push(*c);
                it.next();
            }

            Var(var.into())
        } else if c == '1' {
            it.next();
            Ones
        } else {
            return Err("Unexpected character".to_string());
        };

        loop {
            let c = match it.peek() {
                None => return Ok(e),
                Some(c) => *c,
            };

            let op_pre = match c {
                '|' => 1,
                '^' => 2,
                '&' => 3,
                _ => return Ok(e),
            };

            if op_pre <= pre {
                return Ok(e);
            }

            // If the current operators precedence is higher than
            // the one whose subexpressions we are currently parsing
            // then we need to finish this operator first.
            it.next();
            let rhs = Box::new(Self::parse(it, op_pre)?);
            let lhs = Box::new(e);
            e = match c {
                '&' => And(lhs, rhs),
                '|' => Or(lhs, rhs),
                '^' => Xor(lhs, rhs),
                _ => return Err("Unexpected operator".to_string()),
            };
        }
    }

    /// Converts a [`BExpr`] to an [`ExprOp`].
    pub fn to_expr<R: Ring>(&self, ring: &R) -> ExprOp<R> {
        use BExpr::*;
        match self {
            Ones => ExprOp::Const(ring.negative_one()),
            Var(c) => ExprOp::Var(*c),
            And(l, r) => ExprOp::And(
                l.to_expr(ring).into(),
                r.to_expr(ring).into()
            ),
            Or(l, r) => ExprOp::Or(
                l.to_expr(ring).into(),
                r.to_expr(ring).into()
            ),
            Xor(l, r) => ExprOp::Xor(
                l.to_expr(ring).into(),
                r.to_expr(ring).into()
            ),
            Not(e) => ExprOp::Not(e.to_expr(ring).into()),
        }
    }
}

impl std::fmt::Display for BExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display(crate::formatter::Formatter::C).fmt(f)
    }
}