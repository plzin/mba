#![allow(dead_code)]

use std::{rc::Rc, collections::BTreeSet, ops::{Deref, DerefMut}, fmt::Display};

use num_traits::{Zero, One};
use rug::Integer;

use crate::valuation::Valuation;

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Expr(Rc<ExprOp>);

impl Expr {
    pub fn new(e: ExprOp) -> Self {
        Self(Rc::new(e))
    }

    pub fn as_ref(&self) -> &ExprOp {
        self.0.as_ref()
    }

    pub fn as_ptr(&self) -> *const ExprOp {
        Rc::as_ptr(&self.0)
    }

    pub fn strong_count(&self) -> usize {
        Rc::strong_count(&self.0)
    }

    pub fn ptr_eq(&self, other: &Self) -> bool {
        Rc::ptr_eq(&self.0, &other.0)
    }

    /// Simplify the expression.
    /// This does only trivial simplifications.
    /// E.g. multiplication by one or zero.
    pub fn simplify(&mut self) {
        let mut v = Vec::new();
        simplify_impl(self, &mut v);

        fn simplify_impl(e: &mut Expr, visited: &mut Vec<*const ExprOp>) {
            let ptr = e.as_ptr();
            if e.strong_count() > 1 && visited.contains(&ptr) {
                return
            }

            let m = unsafe { &mut *(ptr as *mut _) };

            match m {
                ExprOp::Const(_) | ExprOp::Var(_) => {},
                ExprOp::Add(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = r.clone();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Sub(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = ExprOp::Neg(r.clone()).into();
                        simplify_impl(e, visited);
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Mul(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() || r.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if l.is_one() {
                        *e = r.clone();
                    } else if r.is_one() {
                        *e = l.clone();
                    }
                },
                ExprOp::Div(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if r.is_one() {
                        *e = l.clone();
                    }
                },
                ExprOp::Neg(i) => {
                    simplify_impl(i, visited);
                    if i.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if let ExprOp::Neg(i) = i.as_ref() {
                        *e = i.clone();
                    }
                },
                ExprOp::And(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() || r.is_zero() {
                        *e = ExprOp::zero().into();
                    }
                },
                ExprOp::Or(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = r.clone();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Xor(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = r.clone();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Not(i) => {
                    simplify_impl(i, visited);
                },
                ExprOp::Shl(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Shr(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
                ExprOp::Sar(l, r) => {
                    simplify_impl(l, visited);
                    simplify_impl(r, visited);
                    if l.is_zero() {
                        *e = ExprOp::zero().into();
                    } else if r.is_zero() {
                        *e = l.clone();
                    }
                },
            }

            if e.strong_count() > 1 {
                visited.push(ptr);
            }
        }
    }

    /// Evaluate an expression.
    pub fn eval(&self, v: &mut Valuation, bits: u32) -> Integer {
        let mut cache = Vec::new();
        return eval_impl(self, v, bits, &mut cache);

        fn eval_impl(
            e: &Expr,
            v: &mut Valuation,
            bits: u32,
            cache: &mut Vec<(*const ExprOp, Integer)>
        ) -> Integer {
            if e.strong_count() > 1 {
                // This is a common subexpression.
                // We don't want to evaluate it twice.
                // So we look it up in the cache.
                let ptr = e.as_ptr();
                for (p, i) in cache.iter() {
                    if *p == ptr {
                        return i.clone();
                    }
                }
            }

            let v = match e.as_ref() {
                ExprOp::Const(n) => n.clone(),
                ExprOp::Var(name) => v.value(name).clone(),
                ExprOp::Add(l, r) => eval_impl(l, v, bits, cache) + eval_impl(r, v, bits, cache),
                ExprOp::Sub(l, r) => eval_impl(l, v, bits, cache) - eval_impl(r, v, bits, cache),
                ExprOp::Mul(l, r) => eval_impl(l, v, bits, cache) * eval_impl(r, v, bits, cache),
                ExprOp::Div(l, r) => {
                    let r = eval_impl(r, v, bits, cache);
                    // We don't want division by zero to panic,
                    // so we define it as zero.
                    if r.is_zero() {
                        Integer::from(0)
                    } else {
                        eval_impl(l, v, bits, cache) / r
                    }
                },
                ExprOp::Neg(i) => -eval_impl(i, v, bits, cache),
                ExprOp::And(l, r) => eval_impl(l, v, bits, cache) & eval_impl(r, v, bits, cache),
                ExprOp::Or(l, r) => eval_impl(l, v, bits, cache) | eval_impl(r, v, bits, cache),
                ExprOp::Xor(l, r) => eval_impl(l, v, bits, cache) ^ eval_impl(r, v, bits, cache),
                ExprOp::Not(i) => !eval_impl(i, v, bits, cache),
                ExprOp::Shl(l, r) => eval_impl(l, v, bits, cache)
                    << (eval_impl(r, v, bits, cache) % bits).to_u32_wrapping(),
                ExprOp::Shr(l, r) => eval_impl(l, v, bits, cache)
                    >> (eval_impl(r, v, bits, cache) % bits).to_u32_wrapping(),
                ExprOp::Sar(l, r) => eval_impl(l, v, bits, cache).keep_signed_bits(bits)
                    >> (eval_impl(r, v, bits, cache) % bits).to_u32_wrapping(),
            };

            let v = v.keep_bits(bits);

            if e.strong_count() > 1 {
                // This is a common subexpression.
                // We want to cache it.
                cache.push((e.as_ptr(), v.clone()));
            }

            v
        }
    }

    /// Substitutes an expression for a variable.
    /// If the Rc's in this expression are shared with
    /// other expressions then this will also substitute in those.
    pub fn substitute(&mut self, var: &str, s: &mut Expr) {
        let mut visited = Vec::new();
        substitute_impl(self, var, s, &mut visited);
        fn substitute_impl(
            e: &mut Expr,
            var: &str,
            s: &mut Expr,
            visited: &mut Vec<*const ExprOp>
        ) {
            let ptr = e.as_ptr();
            let recurse = if visited.contains(&ptr) {
                false
            } else {
                visited.push(ptr);
                true
            };

            use ExprOp::*;
            // SAFETY: This is okay because we make sure with extra logic
            // that this is never encountered twice.
            match unsafe { &mut *(ptr as *mut _) } {
                Const(_) => (),
                Var(v) => if *v == var { *e = s.clone() },
                Add(l, r) | Sub(l, r) | Mul(l, r) | Div(l, r)
                | And(l, r) | Or(l, r) | Xor(l, r) | Shl(l, r)
                | Shr(l, r) | Sar(l, r) => if recurse {
                    substitute_impl(l, var, s, visited);
                    substitute_impl(r, var, s, visited);
                },
                Neg(i) | Not(i) => if recurse {
                    substitute_impl(i, var, s, visited)
                },
            }
        }
    }

    /// Creates a copy of the expression where
    /// no subexpression is shared among the copies.
    pub fn deep_copy(&self) -> Self {
        let mut v = Vec::new();
        return deep_copy_impl(self, &mut v);

        fn deep_copy_impl(e: &Expr, v: &mut Vec<(*const ExprOp, Expr)>) -> Expr {
            if e.strong_count() > 1 {
                if let Some((_, r)) = v.iter().find(|(p, _)| *p == e.as_ptr()) {
                    return r.clone();
                }
            }

            let r = match e.as_ref() {
                ExprOp::Const(n) => ExprOp::Const(n.clone()),
                ExprOp::Var(name) => ExprOp::Var(name.clone()),
                ExprOp::Add(l, r) => ExprOp::Add(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Sub(l, r) => ExprOp::Sub(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Mul(l, r) => ExprOp::Mul(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Div(l, r) => ExprOp::Div(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Neg(i) => ExprOp::Neg(deep_copy_impl(i, v)),
                ExprOp::And(l, r) => ExprOp::And(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Or(l, r) => ExprOp::Or(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Xor(l, r) => ExprOp::Xor(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Not(i) => ExprOp::Not(deep_copy_impl(i, v)),
                ExprOp::Shl(l, r) => ExprOp::Shl(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Shr(l, r) => ExprOp::Shr(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Sar(l, r) => ExprOp::Sar(deep_copy_impl(l, v), deep_copy_impl(r, v)),
            };

            let r = Expr::new(r);

            if e.strong_count() > 1 {
                v.push((e.as_ptr(), r.clone()));
            }

            r
        }
    }

    /// Parse an expression from a string.
    pub fn from_string<T: Into<String>>(s: T) -> Option<Expr> {
        let mut s = s.into();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    // pre 0: parse as much as possible
    // ...
    // pre 15: parse as little as possible
    fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize
    ) -> Option<Expr> {
        use ExprOp::*;

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
            Expr::new(Not(e))
        } else if c == '-' {
            it.next();
            let e = Self::parse(it, 15)?;
            Expr::new(Neg(e))
        } else if c.is_alphabetic() {
            it.next();
            let mut var = String::from(c);
            loop {
                let Some(c) = it.peek() else {
                    break
                };

                if !c.is_alphanumeric() && *c != '_' {
                    break
                }

                var.push(*c);
                it.next();
            }

            Expr::new(Var(var))
        } else if c.is_ascii_digit() {
            let num = crate::int_from_it(it)?;
            Expr::new(Const(num))
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
                '<' | '>' => 4,
                '+' | '-' => 5,
                '*' | '/' => 6,
                ')' => return Some(e),
                _ => return None,
            };

            if op_pre <= pre {
                return Some(e);
            }

            // If the current operators precedence is higher than
            // the one whose subexpression we are currently parsing
            // then we need to finish this operator first.
            it.next();
            let mut rhs = |it| Self::parse(it, op_pre);
            let lhs = e;
            e = Expr::new(match c {
                '+' => Add(lhs, rhs(it)?),
                '-' => Sub(lhs, rhs(it)?),
                '*' => Mul(lhs, rhs(it)?),
                '/' => Div(lhs, rhs(it)?),
                '&' => And(lhs, rhs(it)?),
                '|' => Or(lhs, rhs(it)?),
                '^' => Xor(lhs, rhs(it)?),
                '<' => {
                    if it.next()? != '<' {
                        return None;
                    }

                    Shl(lhs, rhs(it)?)
                },
                '>' => {
                    if it.next()? != '>' {
                        return None;
                    }

                    if it.peek() == Some(&'>') {
                        it.next();
                        Sar(lhs, rhs(it)?)
                    } else {
                        Shr(lhs, rhs(it)?)
                    }
                },
                _ => unreachable!(),
            });
        }
    }

}

impl From<ExprOp> for Expr {
    fn from(e: ExprOp) -> Self {
        Self::new(e)
    }
}

impl Deref for Expr {
    type Target = ExprOp;

    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

impl Display for Expr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

/// A general expression.
/// We use Rc instead of Box in order to avoid copying common subexpressions.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum ExprOp {
    Const(Integer),
    Var(String),
    Add(Expr, Expr),
    Sub(Expr, Expr),
    Mul(Expr, Expr),
    /// This is an unsigned div for now.
    /// But we don't change this during obfuscation,
    /// so if you input an expression where this div is signed,
    /// it should still work.
    /// This is mostly used for evaluating the expression.
    Div(Expr, Expr),
    Neg(Expr),
    And(Expr, Expr),
    Or(Expr, Expr),
    Xor(Expr, Expr),
    Not(Expr),
    Shl(Expr, Expr),
    Shr(Expr, Expr),
    Sar(Expr, Expr),
}


impl ExprOp {
    /// Returns the zero constant.
    pub const fn zero() -> Self {
        Self::Const(Integer::new())
    }

    /// Is this the constant zero?
    pub fn is_zero(&self) -> bool {
        matches!(self, Self::Const(n) if n.is_zero())
    }

    /// Is this the constant one?
    pub fn is_one(&self) -> bool {
        matches!(self, Self::Const(n) if n.is_one())
    }

    /// Returns all variables in the expression.
    /// This can include duplicates.
    pub fn vars(&self) -> Vec<String> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<String>) {
        match self {
            ExprOp::Const(_) => {},
            ExprOp::Var(name) => drop(v.insert(name.clone())),
            ExprOp::Neg(e) | ExprOp::Not(e) => e.vars_impl(v),
            ExprOp::Add(l, r) | ExprOp::Sub(l, r) | ExprOp::Mul(l, r)
            | ExprOp::Div(l, r) | ExprOp::And(l, r) | ExprOp::Or(l, r)
            | ExprOp::Xor(l, r) | ExprOp::Shl(l, r) | ExprOp::Shr(l, r)
            | ExprOp::Sar(l, r) => {
                l.vars_impl(v);
                r.vars_impl(v);
            }
        }
    }

    fn print_simple_rc(
        e: &Expr,
        vars: &mut Vec<(*const ExprOp, char, String)>
    ) -> String {
        // If there is only one reference then just print it.
        if e.strong_count() == 1 {
            return e.print_simple_impl(vars);
        }

        // We don't want to assign a variable to a variable
        // so there is this shortcut here.
        if let ExprOp::Var(v) = e.as_ref() {
            return format!("{}", *v);
        }

        let ptr = e.as_ptr();

        // If the expression already has a variable then just print the variable.
        let var = vars.iter().find(|t| t.0 == ptr);
        if let Some(v) = var {
            v.1.to_string()
        } else {
            let v = match vars.last() {
                None => 'a',
                // TODO: If there are too many variables this is bad.
                Some(t) => (t.1 as u8 + 1) as char,
            };

            // Push everything.
            vars.push((ptr, v, String::new()));

            let idx = vars.len() - 1;

            // Get the initializer for the variable.
            vars[idx].2 = e.print_simple_impl(vars);

            // Return just the variable name.
            v.to_string()
        }
    }

    // Yes, this PERFORMANCE CRITICAL code could be more efficient...
    fn print_simple_impl(
        &self, vars: &mut Vec<(*const ExprOp, char, String)>
    ) -> String {
        let bin_op = |
            op: &str, l: &Expr, r: &Expr,
            vars: &mut Vec<(*const ExprOp, char, String)>
        | {
            let pred = self.precedence();

            let l = if pred > l.precedence() && l.strong_count() == 1 {
                format!("({})", ExprOp::print_simple_rc(l, vars))
            } else {
                format!("{}", ExprOp::print_simple_rc(l, vars))
            };

            let r = if pred > r.precedence() && r.strong_count() == 1 {
                format!("({})", ExprOp::print_simple_rc(r, vars))
            } else {
                format!("{}", ExprOp::print_simple_rc(r, vars))
            };

            format!("{} {} {}", l, op, r)
        };

        let un_op = |
            op: &str, i: &Expr,
            vars: &mut Vec<(*const ExprOp, char, String)>
        | {
            if self.precedence() > i.precedence() && i.strong_count() == 1 {
                format!("{}({})", op, ExprOp::print_simple_rc(i, vars))
            } else {
                format!("{}{}", op, ExprOp::print_simple_rc(i, vars))
            }
        };

        use ExprOp::*;
        match self {
            Const(i) => format!("{}", i),
            Var(n) => format!("{}", n),
            Add(l, r) => bin_op("+", l, r, vars),
            Sub(l, r) => bin_op("-", l, r, vars),
            Mul(l, r) => bin_op("*", l, r, vars),
            Div(l, r) => bin_op("/", l, r, vars),
            Neg(i) => un_op("-", i, vars),
            And(l, r) => bin_op("&", l, r, vars),
            Or(l, r) => bin_op("|", l, r, vars),
            Xor(l, r) => bin_op("^", l, r, vars),
            Not(i) => un_op("~", i, vars),
            Shl(l, r) => bin_op("<<", l, r, vars),
            Shr(l, r) => bin_op(">>", l, r, vars),
            Sar(l, r) => bin_op(">>>", l, r, vars),
        }
    }

    /// Returns the precedence of a binary operator.
    /// All operators are taken to be left associative.
    fn precedence(&self) -> usize {
        use ExprOp::*;
        match self {
            Or(_, _) => 1,
            Xor(_, _) => 2,
            And(_, _) => 3,
            Shl(_, _) | Shr(_, _) | Sar(_, _) => 4,
            Add(_, _) | Sub(_, _) => 5,
            Mul(_, _) | Div(_, _) => 6,
            Neg(_) | Not(_) => 15,
            Const(_) | Var(_) => 16,
        }
    }

    /// This does not currently share common subexpressions.
    #[cfg(feature = "z3")]
    pub fn to_z3_bv<'ctx>(&self, ctx: &'ctx z3::Context, width: u32) -> z3::ast::BV<'ctx> {
        use ExprOp::*;

        match self {
            Const(i) => crate::int_to_bv(ctx, width, i),
            Var(v) => z3::ast::BV::new_const(ctx, v.as_str(), width),
            Add(l, r) => l.to_z3_bv(ctx, width).bvadd(&r.to_z3_bv(ctx, width)),
            Sub(l, r) => l.to_z3_bv(ctx, width).bvsub(&r.to_z3_bv(ctx, width)),
            Mul(l, r) => l.to_z3_bv(ctx, width).bvmul(&r.to_z3_bv(ctx, width)),
            Div(l, r) => l.to_z3_bv(ctx, width).bvudiv(&r.to_z3_bv(ctx, width)),
            Neg(i) => i.to_z3_bv(ctx, width).bvneg(),
            And(l, r) => l.to_z3_bv(ctx, width).bvand(&r.to_z3_bv(ctx, width)),
            Or(l, r) => l.to_z3_bv(ctx, width).bvor(&r.to_z3_bv(ctx, width)),
            Xor(l, r) => l.to_z3_bv(ctx, width).bvxor(&r.to_z3_bv(ctx, width)),
            Not(i) => i.to_z3_bv(ctx, width).bvnot(),
            Shl(l, r) => l.to_z3_bv(ctx, width).bvshl(&r.to_z3_bv(ctx, width)),
            Shr(l, r) => l.to_z3_bv(ctx, width).bvlshr(&r.to_z3_bv(ctx, width)),
            Sar(l, r) => l.to_z3_bv(ctx, width).bvashr(&r.to_z3_bv(ctx, width)),
        }
    }
}

impl Display for ExprOp {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Stores a mapping of (sub)expressions to variables.
        let mut vars = Vec::new();

        let l = self.print_simple_impl(&mut vars);

        for (_, var, init) in vars.iter().rev() {
            writeln!(f, "{} = {}", var, init)?;
        }

        writeln!(f, "{}", l)
    }
}

#[test]
fn parse_expr_test() {
    let e = Expr::from_string("1 + 2 * 3").unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Add(
        ExprOp::Const(1.into()).into(),
        ExprOp::Mul(
            ExprOp::Const(2.into()).into(),
            ExprOp::Const(3.into()).into(),
        ).into(),
    ));

    let e = Expr::from_string("x << y + z").unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Shl(
        ExprOp::Var("x".into()).into(),
        ExprOp::Add(
            ExprOp::Var("y".into()).into(),
            ExprOp::Var("z".into()).into(),
        ).into(),
    ));

    let e = Expr::from_string("x >> 3 + z").unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Shr(
        ExprOp::Var("x".into()).into(),
        ExprOp::Add(
            ExprOp::Const(3.into()).into(),
            ExprOp::Var("z".into()).into(),
        ).into(),
    ));

    let e = Expr::from_string("x >>> (x + long_name) & 4").unwrap();
    assert_eq!(e.as_ref(), &ExprOp::And(
        ExprOp::Sar(
            ExprOp::Var("x".into()).into(),
            ExprOp::Add(
                ExprOp::Var("x".into()).into(),
                ExprOp::Var("long_name".into()).into(),
            ).into(),
        ).into(),
        ExprOp::Const(4.into()).into(),
    ));

    let e = Expr::from_string("x & y | z ^ 1").unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Or(
        ExprOp::And(
            ExprOp::Var("x".into()).into(),
            ExprOp::Var("y".into()).into(),
        ).into(),
        ExprOp::Xor(
            ExprOp::Var("z".into()).into(),
            ExprOp::Const(1.into()).into(),
        ).into(),
    ));
}

#[test]
fn expr_eval_test() {
    let e = Expr::from_string("1 + 2 * 3").unwrap();
    assert_eq!(e.eval(&mut Valuation::empty(), 8), 7);

    let e = Expr::from_string("(x << y) + z").unwrap();
    assert_eq!(e.eval(&mut Valuation::from_vec_panic(vec![
        ("x".into(), 1.into()),
        ("y".into(), 3.into()),
        ("z".into(), 7.into()),
    ]), 8), 15);
}
