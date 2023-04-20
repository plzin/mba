#![allow(dead_code)]

use std::{ops::{Index, IndexMut, Neg}, collections::BTreeSet};

use rug::Integer;

use crate::{Expr, int_from_it};

/// LUExpr is short for "Linear combination of Uniform Expressions"
/// These are the expressions for which rewrite rules can be efficiently
/// generated.
#[derive(Clone, Debug)]
pub struct LUExpr(pub Vec<(Integer, UExpr)>);

impl LUExpr {
    /// Creates an expression that equals a constant.
    pub fn constant(c: Integer) -> Self {
        Self(vec![(-c, UExpr::Ones)])
    }

    /// Creates an expression that equals a variable.
    pub fn var(name: String) -> Self {
        Self(vec![(1.into(), UExpr::Var(name))])
    }

    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self) -> Vec<String> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<String>) {
        for (_, e) in &self.0 {
            e.vars_impl(v);
        }
    }

    /// Evaluate an expression with a valuation for the occurring variables.
    pub fn eval(&self, v: &Valuation) -> Integer {
        self.0.iter()
            .map(|(i, e)| i * e.eval(v))
            .sum()
    }

    /// Parse a string to an expression.
    /// Note that this function is extremely limited
    /// and expects very specific syntax.
    /// It is used for convenience when testing things and
    /// not really meant to be used by something outside this crate.
    pub(crate) fn from_string<T: ToString>(s: T) -> Option<Self> {
        let mut s = s.to_string();
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

            // If the next character is not a plus then we are done.
            match it.peek() {
                Some('+') => it.next(), // Skip the +.
                Some('-') => { neg = true; it.next() },
                _ => return Some(Self(v)),
            };
        }
    }

    /// Converts an LUExpr to an Expr.
    pub fn to_expr(&self) -> Expr {
        let mut it = self.0.iter();
        let s = match it.next() {
            None => return Expr::Const(Integer::new()),
            Some(s) => s,
        };

        let from_summand = |(i, e): &(Integer, UExpr)| {
            if let UExpr::Ones = e {
                Expr::Const(-i.clone())
            } else {
                Expr::Mul(Expr::Const(i.clone()).into(), e.to_expr().into())
            }
        };

        let mut cur = from_summand(s);
        for s in it {
            cur = Expr::Add(cur.into(), from_summand(s).into());
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
    Var(String),
    And(Box<Self>, Box<Self>),
    Or(Box<Self>, Box<Self>),
    Xor(Box<Self>, Box<Self>),
    Not(Box<Self>),
}

impl UExpr {
    pub fn var(c: String) -> Self {
        Self::Var(c)
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

    pub fn not(e: Self) -> Self {
        Self::Not(e.into())
    }

    /// Returns all variables in the expression.
    /// This will include duplicates.
    pub fn vars(&self) -> Vec<String> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<String>) {
        use UExpr::*;
        match self {
            Ones            => {},
            Var(c)          => drop(v.insert(c.clone())),
            And(e1, e2)     => { e1.vars_impl(v); e2.vars_impl(v) },
            Or(e1, e2)      => { e1.vars_impl(v); e2.vars_impl(v) },
            Xor(e1, e2)     => { e1.vars_impl(v); e2.vars_impl(v) },
            Not(e)          => e.vars_impl(v),
        }
    }

    /// Evaluate an expression with a valuation for the occurring variables.
    pub fn eval(&self, v: &Valuation) -> Integer {
        use UExpr::*;
        match self {
            Ones            => (-1).into(),
            Var(c)          => v[c].clone(),
            And(e1, e2)     => e1.eval(v) & e2.eval(v),
            Or(e1, e2)      => e1.eval(v) | e2.eval(v),
            Xor(e1, e2)     => e1.eval(v) ^ e2.eval(v),
            Not(e)          => !e.eval(v),
        }
    }

    /// Rename a variable.
    pub fn rename_var(&mut self, old: &str, new: &str) {
        use UExpr::*;
        match self {
            Ones        => (),
            Var(v)      => if *v == old { v.clear(); v.push_str(new) },
            And(l, r)   => { l.rename_var(old, new); r.rename_var(old, new) },
            Or(l, r)    => { l.rename_var(old, new); r.rename_var(old, new) },
            Xor(l, r)   => { l.rename_var(old, new); r.rename_var(old, new) },
            Not(e)      => e.rename_var(old, new),
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
    pub(crate) fn from_string<T: ToString>(s: T) -> Option<Self> {
        let mut s = s.to_string();
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

            Var(var)
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
    pub fn to_expr(&self) -> Expr {
        use UExpr::*;
        match self {
            Ones        => Expr::Const((-1).into()),
            Var(c)      => Expr::Var(c.clone()),
            And(l, r)   => Expr::And(l.to_expr().into(), r.to_expr().into()),
            Or(l, r)    => Expr::Or(l.to_expr().into(), r.to_expr().into()),
            Xor(l, r)   => Expr::Xor(l.to_expr().into(), r.to_expr().into()),
            Not(e)      => Expr::Not(e.to_expr().into()),
        }
    }
}

impl std::fmt::Display for UExpr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        use UExpr::*;
        match self {
            Ones => write!(f, "1"),
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

/// Stores values that should be substituted into variables.
#[derive(Debug)]
pub struct Valuation {
    /// The key value pairs are stored as a Vector
    /// because I doubt a hashmap/tree would be faster
    /// when there are so few variables.
    vals: Vec<(String, Integer)>,
}

impl Valuation {
    /// Initializes a valuation from a list of variables
    /// each of which will be Initialized to 0.
    pub fn zero(vars: Vec<String>) -> Self {
        let vals = vars.into_iter()
            .map(|c| (c, 0.into()))
            .collect();

        Self { vals }
    }
}

impl Index<&str> for Valuation {
    type Output = Integer;
    fn index(&self, index: &str) -> &Self::Output {
        &self.vals.iter()
            .find(|(name, _)| name == index)
            .unwrap().1
    }
}

impl IndexMut<&str> for Valuation {
    fn index_mut(&mut self, index: &str) -> &mut Self::Output {
        &mut self.vals.iter_mut()
            .find(|(name, _)| name == index)
            .unwrap().1
    }
}
