//! General expressions.

use std::rc::Rc;
use std::collections::BTreeSet;
use std::ops::Deref;
use crate::rings::{Ring, RingElement as _, BinaryRing};
use crate::valuation::Valuation;
use crate::Symbol;

#[derive(Clone, PartialEq, Debug)]
pub struct Expr<R: Ring>(Rc<ExprOp<R>>);

impl<R: Ring> Expr<R> {
    pub fn new(e: ExprOp<R>) -> Self {
        Self(Rc::new(e))
    }

    #[allow(clippy::should_implement_trait)]
    pub fn as_ref(&self) -> &ExprOp<R> {
        self.0.as_ref()
    }

    pub fn as_ptr(&self) -> *const ExprOp<R> {
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

        fn simplify_impl<R: Ring>(
            e: &mut Expr<R>,
            visited: &mut Vec<*const ExprOp<R>>
        ) {
            let ptr = e.as_ptr();
            if e.strong_count() > 1 {
                if visited.contains(&ptr) {
                    return
                }
                visited.push(ptr);
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
            }
        }
    }

    /// Evaluate an expression.
    pub fn eval(&self, v: &mut Valuation<R>, r: &R) -> R::Element
    where
        R: BinaryRing
    {
        let mut cache = Vec::new();
        return eval_impl(self, v, r, &mut cache);

        fn eval_impl<R: BinaryRing>(
            e: &Expr<R>,
            v: &mut Valuation<R>,
            ring: &R,
            cache: &mut Vec<(*const ExprOp<R>, R::Element)>
        ) -> R::Element {
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
                ExprOp::Var(name) => v.value(*name, ring).clone(),
                ExprOp::Add(l, r) => ring.add(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Sub(l, r) => ring.sub(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Mul(l, r) => ring.mul(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Neg(i) => ring.neg(eval_impl(i, v, ring, cache)),
                ExprOp::And(l, r) => R::and(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Or(l, r) => R::or(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Xor(l, r) => R::xor(eval_impl(l, v, ring, cache), &eval_impl(r, v, ring, cache)),
                ExprOp::Not(i) => ring.not(eval_impl(i, v, ring, cache)),
            };

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
    pub fn substitute(&mut self, var: Symbol, s: &mut Expr<R>) {
        let mut visited = Vec::new();
        substitute_impl(self, var, s, &mut visited);
        fn substitute_impl<R: Ring>(
            e: &mut Expr<R>,
            var: Symbol,
            s: &mut Expr<R>,
            visited: &mut Vec<*const ExprOp<R>>
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
                Add(l, r) | Sub(l, r) | Mul(l, r)
                | And(l, r) | Or(l, r) | Xor(l, r) => if recurse {
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

        fn deep_copy_impl<R: Ring>(
            e: &Expr<R>,
            v: &mut Vec<(*const ExprOp<R>, Expr<R>)>
        ) -> Expr<R> {
            if e.strong_count() > 1
                && let Some((_, r)) = v.iter().find(|(p, _)| *p == e.as_ptr())
            {
                return r.clone();
            }

            let r = match e.as_ref() {
                ExprOp::Const(n) => ExprOp::Const(n.clone()),
                ExprOp::Var(name) => ExprOp::Var(*name),
                ExprOp::Add(l, r) => ExprOp::Add(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Sub(l, r) => ExprOp::Sub(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Mul(l, r) => ExprOp::Mul(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Neg(i) => ExprOp::Neg(deep_copy_impl(i, v)),
                ExprOp::And(l, r) => ExprOp::And(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Or(l, r) => ExprOp::Or(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Xor(l, r) => ExprOp::Xor(deep_copy_impl(l, v), deep_copy_impl(r, v)),
                ExprOp::Not(i) => ExprOp::Not(deep_copy_impl(i, v)),
            };

            let r = Expr::new(r);

            if e.strong_count() > 1 {
                v.push((e.as_ptr(), r.clone()));
            }

            r
        }
    }

    /// Parse an expression from a string.
    pub fn from_string(mut s: String, r: &R) -> Result<Expr<R>, String> {
        s.retain(|c| !matches!(c, ' ' | '\t' | '\n'));
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0, r)
    }

    // pre 0: parse as much as possible
    // ...
    // pre 15: parse as little as possible
    fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize,
        r: &R,
    ) -> Result<Expr<R>, String> {
        use ExprOp::*;

        let c = *it.peek().ok_or("Unexpected end of input")?;

        let mut e = if c == '(' {
            it.next();
            let e = Self::parse(it, 0, r)?;
            match it.next() {
                Some(')') => e,
                _ => return Err("Expected closing parenthesis".to_string()),
            }
        } else if c == '~' || c == '!' {
            it.next();
            let e = Self::parse(it, 15, r)?;
            Expr::new(Not(e))
        } else if c == '-' {
            it.next();
            let e = Self::parse(it, 15, r)?;
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

            Expr::new(Var(var.into()))
        } else if c.is_ascii_digit() {
            Expr::new(Const(r.parse_element(it).ok_or("Expected number")?))
        } else {
            return Err("Expected variable".to_string());
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
                '+' | '-' => 5,
                '*' => 6,
                ')' => return Ok(e),
                _ => return Err("Unexpected character".to_string()),
            };

            if op_pre <= pre {
                return Ok(e);
            }

            // If the current operators precedence is higher than
            // the one whose subexpression we are currently parsing
            // then we need to finish this operator first.
            it.next();
            let rhs = |it| Self::parse(it, op_pre, r);
            let lhs = e;
            e = Expr::new(match c {
                '+' => Add(lhs, rhs(it)?),
                '-' => Sub(lhs, rhs(it)?),
                '*' => Mul(lhs, rhs(it)?),
                '&' => And(lhs, rhs(it)?),
                '|' => Or(lhs, rhs(it)?),
                '^' => Xor(lhs, rhs(it)?),
                _ => unreachable!(),
            });
        }
    }

}

impl<R: Ring> From<ExprOp<R>> for Expr<R> {
    fn from(e: ExprOp<R>) -> Self {
        Self::new(e)
    }
}

impl<R: Ring> Deref for Expr<R> {
    type Target = ExprOp<R>;

    fn deref(&self) -> &Self::Target {
        self.as_ref()
    }
}

/// A general expression.
/// We use Rc instead of Box in order to avoid copying common subexpressions.
/// The semantics of this as used in [Expr::eval] are questionable.
#[derive(Clone, PartialEq, Debug)]
pub enum ExprOp<R: Ring> {
    Const(R::Element),
    Var(Symbol),
    Add(Expr<R>, Expr<R>),
    Sub(Expr<R>, Expr<R>),
    Mul(Expr<R>, Expr<R>),
    Neg(Expr<R>),
    And(Expr<R>, Expr<R>),
    Or(Expr<R>, Expr<R>),
    Xor(Expr<R>, Expr<R>),
    Not(Expr<R>),
}


impl<R: Ring> ExprOp<R> {
    /// Returns the zero constant.
    pub fn zero() -> Self {
        Self::Const(R::zero())
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
    pub fn vars(&self) -> Vec<Symbol> {
        let mut v = BTreeSet::new();
        self.vars_impl(&mut v);
        v.into_iter().collect()
    }

    pub(crate) fn vars_impl(&self, v: &mut BTreeSet<Symbol>) {
        match self {
            ExprOp::Const(_) => {},
            ExprOp::Var(name) => drop(v.insert(*name)),
            ExprOp::Neg(e) | ExprOp::Not(e) => e.vars_impl(v),
            ExprOp::Add(l, r) | ExprOp::Sub(l, r) | ExprOp::Mul(l, r)
            | ExprOp::And(l, r) | ExprOp::Or(l, r) | ExprOp::Xor(l, r) => {
                l.vars_impl(v);
                r.vars_impl(v);
            }
        }
    }
    /// Returns the precedence of a binary operator.
    /// All operators are taken to be left associative.
    pub(crate) fn precedence(&self) -> usize {
        use ExprOp::*;
        match self {
            Or(_, _) => 1,
            Xor(_, _) => 2,
            And(_, _) => 3,
            Add(_, _) | Sub(_, _) => 5,
            Mul(_, _) => 6,
            Neg(_) | Not(_) => 15,
            Const(_) | Var(_) => 16,
        }
    }

    /// This does not currently share common subexpressions.
    #[cfg(feature = "z3")]
    pub fn to_z3_bv<'ctx>(&self, ctx: &'ctx z3::Context, ring: &R) -> z3::ast::BV<'ctx>
    where
        R: BinaryRing
    {
        use ExprOp::*;

        match self {
            Const(i) => crate::int_to_bv(ctx, i, ring),
            Var(v) => z3::ast::BV::new_const(ctx, v.as_str(), ring.bits()),
            Add(l, r) => l.to_z3_bv(ctx, ring).bvadd(&r.to_z3_bv(ctx, ring)),
            Sub(l, r) => l.to_z3_bv(ctx, ring).bvsub(&r.to_z3_bv(ctx, ring)),
            Mul(l, r) => l.to_z3_bv(ctx, ring).bvmul(&r.to_z3_bv(ctx, ring)),
            Neg(i) => i.to_z3_bv(ctx, ring).bvneg(),
            And(l, r) => l.to_z3_bv(ctx, ring).bvand(&r.to_z3_bv(ctx, ring)),
            Or(l, r) => l.to_z3_bv(ctx, ring).bvor(&r.to_z3_bv(ctx, ring)),
            Xor(l, r) => l.to_z3_bv(ctx, ring).bvxor(&r.to_z3_bv(ctx, ring)),
            Not(i) => i.to_z3_bv(ctx, ring).bvnot(),
        }
    }
}

#[test]
fn parse_expr_test() {
    use crate::rings::U32;
    let e = Expr::from_string("1 + 2 * 3".to_owned(), &U32).unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Add(
        ExprOp::Const(1).into(),
        ExprOp::Mul(
            ExprOp::Const(2).into(),
            ExprOp::Const(3).into(),
        ).into(),
    ));

    let e = Expr::from_string("x & y + z".to_owned(), &U32).unwrap();
    assert_eq!(e.as_ref(), &ExprOp::And(
        ExprOp::Var("x".into()).into(),
        ExprOp::Add(
            ExprOp::Var("y".into()).into(),
            ExprOp::Var("z".into()).into(),
        ).into(),
    ));

    let e = Expr::from_string("x | 3 + z".to_owned(), &U32).unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Or(
        ExprOp::Var("x".into()).into(),
        ExprOp::Add(
            ExprOp::Const(3).into(),
            ExprOp::Var("z".into()).into(),
        ).into(),
    ));

    let e = Expr::from_string("x & (x + long_name) ^ 4".to_owned(), &U32)
        .unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Xor(
        ExprOp::And(
            ExprOp::Var("x".into()).into(),
            ExprOp::Add(
                ExprOp::Var("x".into()).into(),
                ExprOp::Var("long_name".into()).into(),
            ).into(),
        ).into(),
        ExprOp::Const(4).into(),
    ));

    let e = Expr::from_string("x & y | z ^ 1".to_owned(), &U32).unwrap();
    assert_eq!(e.as_ref(), &ExprOp::Or(
        ExprOp::And(
            ExprOp::Var("x".into()).into(),
            ExprOp::Var("y".into()).into(),
        ).into(),
        ExprOp::Xor(
            ExprOp::Var("z".into()).into(),
            ExprOp::Const(1).into(),
        ).into(),
    ));
}

#[test]
fn expr_eval_test() {
    use crate::rings::U32;
    let e = Expr::from_string("1 + 2 * 3".to_owned(), &U32).unwrap();
    assert_eq!(e.eval(&mut Valuation::empty(), &U32), 7);

    let e = Expr::from_string("(x * y) + z".to_owned(), &U32).unwrap();
    assert_eq!(e.eval(&mut Valuation::from_vec_panic(vec![
        ("x".into(), 1),
        ("y".into(), 3),
        ("z".into(), 7),
    ]), &U32), 10);
}
