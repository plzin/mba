//! This module contains uniform expressions ([UExpr]), which are just bitwise expressions
//! on n bits, and linear combinations of uniform expressions ([LUExpr]).

use std::{ops::Neg, collections::BTreeSet};
use rug::{Integer, Complete};
use crate::valuation::Valuation;
use crate::{Symbol, ExprOp, int_from_it};

/// LUExpr is short for "Linear combination of Uniform Expressions"
/// These are the expressions for which rewrite rules can be efficiently
/// generated.
#[derive(Clone, Debug)]
pub struct LUExpr(pub Vec<(Integer, UExpr)>);

impl LUExpr {
    /// Create the empty linear combination.
    /// It evaluates to 0.
    pub fn zero() -> Self {
        Self(Vec::new())
    }

    /// Creates an expression that equals a constant.
    pub fn constant(c: Integer) -> Self {
        Self(vec![(-c, UExpr::Ones)])
    }

    /// Creates an expression that equals a variable.
    pub fn var<T: Into<Symbol>>(name: T) -> Self {
        Self(vec![(1.into(), UExpr::Var(name.into()))])
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
    pub fn eval(&self, v: &mut Valuation, bits: u32) -> Integer {
        self.0.iter()
            .map(|(i, e)| i * e.eval(v))
            .fold(Integer::new(), |acc, x| (acc + x).keep_bits(bits))
    }

    /// Returns a measure of the complexity of the expression.
    pub fn complexity(&self, bits: u32) -> u32 {
        // Complexity of a coefficient.
        let coeff_complexity = |i: &Integer| {
            let abs = i.keep_signed_bits_ref(bits).complete().abs();
            abs.count_ones().unwrap() / 2 + abs.significant_bits()
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
    pub fn from_string<T: Into<String>>(s: T) -> Option<Self> {
        let mut s = s.into();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        // This stores the current linear combination.
        let mut v = Vec::new();

        let mut neg = false;

        // Loop over the string/the summands.
        loop {
            // Are there still characters left?
            // If not then we're done.
            let mut c = match it.peek() {
                None => return Some(Self(v)),
                Some(c) => *c,
            };

            if c == '-' {
                neg = true;
                it.next();
                c = *it.peek()?;
            }

            // If this is a digit then we expect num*UExpr.
            if c.is_ascii_digit() {
                // Parse the number.
                let mut num = int_from_it(&mut it)?;

                // If the number is negative then negate it.
                if neg {
                    num = num.neg();
                }

                // Is it the expected '*'?
                let lookahead = it.peek();
                match lookahead {
                    Some('*') => {
                        it.next();

                        // Parse the UExpr.
                        let e = UExpr::parse(&mut it, 0)?;

                        // Push it.
                        v.push((num, e));
                    },

                    // If this is a different character then we push -num*(-1).
                    _ => v.push((-num, UExpr::Ones)),
                }
            } else {
                // We don't have a factor so just parse the UExpr.
                let e = UExpr::parse(&mut it, 0)?;

                let sign = match neg {
                    false => 1,
                    true => -1,
                };

                // Push sign*e.
                v.push((sign.into(), e));
            }

            // If the next character is not a plus or - then we are done.
            match it.peek() {
                Some('+') => { neg = false; it.next() }, // Skip the +.
                Some('-') => { neg = true; it.next() },
                _ => return Some(Self(v)),
            };
        }
    }

    /// Converts an LUExpr to an Expr.
    pub fn to_expr(&self) -> ExprOp {
        let mut it = self.0.iter();
        let s = match it.next() {
            None => return ExprOp::Const(Integer::new()),
            Some(s) => s,
        };

        let from_summand = |(i, e): &(Integer, UExpr)| {
            if let UExpr::Ones = e {
                ExprOp::Const(-i.clone())
            } else {
                ExprOp::Mul(ExprOp::Const(i.clone()).into(), e.to_expr().into())
            }
        };

        let mut cur = from_summand(s);
        for s in it {
            cur = ExprOp::Add(cur.into(), from_summand(s).into());
        }

        cur
    }
}

impl From<UExpr> for LUExpr {
    fn from(e: UExpr) -> Self {
        Self(vec![(1.into(), e)])
    }
}

impl std::fmt::Display for LUExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.0.iter().filter(|(i, _)| *i != 0);
        let (i, e) = match iter.next() {
            Some(t) => t,
            None => return write!(f, "0"),
        };

        if *i == 1 {
            write!(f, "{}", e)?;
        } else {
            write!(f, "{}*({})", i, e)?;
        }

        for (i, e) in iter {
            write!(f, " + ")?;
            if *i == 1 {
                write!(f, "{}", e)?;
            } else {
                write!(f, "{}*({})", i, e)?;
            }
        }

        Ok(())
    }
}

/// Represents an expression that is uniform on all bits.
/// Note that the variant 'Ones' does not equal 1, but a 1 in every bit,
/// which is -1 in two's complement.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum UExpr {
    Ones,
    Var(Symbol),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Xor(Box<Self>, Box<Self>),
    Not(Box<Self>),
}

impl UExpr {
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

    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self) -> Vec<Symbol> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<Symbol>) {
        use UExpr::*;
        match self {
            Ones            => {},
            Var(c)          => drop(v.insert(*c)),
            And(e1, e2)     => { e1.vars_impl(v); e2.vars_impl(v) },
            Or(e1, e2)      => { e1.vars_impl(v); e2.vars_impl(v) },
            Xor(e1, e2)     => { e1.vars_impl(v); e2.vars_impl(v) },
            Not(e)          => e.vars_impl(v),
        }
    }

    /// Evaluate an expression with a valuation for the occurring variables.
    pub fn eval(&self, v: &mut Valuation) -> Integer {
        use UExpr::*;
        match self {
            Ones            => (-1).into(),
            Var(c)          => v.value(*c).clone(),
            And(e1, e2)     => e1.eval(v) & e2.eval(v),
            Or(e1, e2)      => e1.eval(v) | e2.eval(v),
            Xor(e1, e2)     => e1.eval(v) ^ e2.eval(v),
            Not(e)          => !e.eval(v),
        }
    }

    /// Rename a variable.
    pub fn rename_var(&mut self, old: Symbol, new: Symbol) {
        use UExpr::*;
        match self {
            Ones        => (),
            Var(v)      => if *v == old { *v = new },
            And(l, r)   => { l.rename_var(old, new); r.rename_var(old, new) },
            Or(l, r)    => { l.rename_var(old, new); r.rename_var(old, new) },
            Xor(l, r)   => { l.rename_var(old, new); r.rename_var(old, new) },
            Not(e)      => e.rename_var(old, new),
        }
    }

    /// Returns some sort of complexity measure of the expression.
    pub fn complexity(&self) -> u32 {
        use UExpr::*;
        match self {
            Ones => 1,
            Var(_) => 1,
            And(l, r) => l.complexity() + r.complexity() + 1,
            Or(l, r) => l.complexity() + r.complexity() + 1,
            Xor(l, r) => l.complexity() + r.complexity() + 1,
            Not(e) => e.complexity() + 1,
        }
    }

    fn write_safe(
        e1: &Self, e2: &Self, op: char, f: &mut std::fmt::Formatter<'_>
    ) -> std::fmt::Result {
        if let Self::Var(c) = e1 {
            write!(f, "{} {}", c, op)?;
        } else {
            write!(f, "({}) {}", e1, op)?;
        }

        if let Self::Var(c) = e2 {
            write!(f, " {}", c)
        } else {
            write!(f, " ({})", e2)
        }
    }

    /// Parse a string to an expression.
    pub fn from_string<T: Into<String>>(s: T) -> Option<Self> {
        let mut s = s.into();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    pub(self) fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize
    ) -> Option<Self> {
        use UExpr::*;

        let c = *it.peek()?;

        let mut e = if c == '(' {
            it.next();
            let e = Self::parse(it, 0)?;
            match it.next() {
                Some(')') => e,
                _ => return None,
            }
        } else if c == '~' {
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
            return None;
        };

        loop {
            let c = match it.peek() {
                None => return Some(e),
                Some(c) => *c,
            };

            let op_pre = match c {
                '|' => 1,
                '^' => 2,
                '&' => 3,
                _ => return Some(e),
            };

            if op_pre <= pre {
                return Some(e);
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
                _ => return None,
            };
        }
    }

    /// Converts a UExpr to an Expr.
    pub fn to_expr(&self) -> ExprOp {
        use UExpr::*;
        match self {
            Ones        => ExprOp::Const((-1).into()),
            Var(c)      => ExprOp::Var(*c),
            And(l, r)   => ExprOp::And(l.to_expr().into(), r.to_expr().into()),
            Or(l, r)    => ExprOp::Or(l.to_expr().into(), r.to_expr().into()),
            Xor(l, r)   => ExprOp::Xor(l.to_expr().into(), r.to_expr().into()),
            Not(e)      => ExprOp::Not(e.to_expr().into()),
        }
    }
}

impl std::fmt::Display for UExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use UExpr::*;
        match self {
            Ones => write!(f, "-1"),
            Var(c) => write!(f, "{}", c),
            And(e1, e2)   => Self::write_safe(e1, e2, '&', f),
            Or(e1, e2)    => Self::write_safe(e1, e2, '|', f),
            Xor(e1, e2)   => Self::write_safe(e1, e2, '^', f),
            Not(e) =>
                if let Self::Var(c) = e.as_ref() {
                    write!(f, "~{}", c)
                } else {
                    write!(f, "~({})", e)
                },
        }
    }
}