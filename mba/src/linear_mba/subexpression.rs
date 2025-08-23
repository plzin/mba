//! Finds linear MBA subexpressions in an [`Expr`].

use crate::{
    Symbol,
    bitwise_expr::{BExpr, LBExpr},
    expr::{Expr, ExprOp},
    rings::BinaryRing,
};

/// A list of substitutions.
/// See e.g. [`expr_to_bexpr`] on what this is used for.
#[derive(Default)]
pub struct Subs<R: BinaryRing>(pub Vec<(Symbol, Expr<R>)>);

impl<R: BinaryRing> Subs<R> {
    /// Creates a new empty substitution list.
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Adds a new substitution.
    pub fn add(&mut self, expr: Expr<R>) -> Symbol {
        // Do we have this expression stored already?
        for (var, e) in &self.0 {
            if e == &expr {
                return *var;
            }
        }

        // Create a new substitution variable.
        let var = Symbol::from(format!("_sub_{}", self.0.len()));
        self.0.push((var, expr));
        var
    }
}

/// Converts part of an expression to a [`BExpr`], such that if you substituted
/// the Exprs in `subs` for the variables, you would get the original Expr.
///
/// It will generally try to make the [`LBExpr`] as big as possible.
///
/// If `force` is false, it will return [`None`] if the top-most operation is
/// not a [`BExpr`] operation. Otherwise, it will return a [`BExpr::Var`] whose
/// substitution is the original expression.
pub fn expr_to_bexpr<R: BinaryRing>(
    e: &Expr<R>,
    subs: &mut Subs<R>,
    force: bool,
) -> Option<BExpr> {
    // New substitution variable.
    let mut new_sub = || force.then(|| BExpr::Var(subs.add(e.clone())));

    if let ExprOp::Var(v) = e.as_ref() {
        return Some(BExpr::Var(*v));
    }

    // We don't try, when the expression is shared.
    if e.strong_count() > 1 {
        return new_sub();
    }

    match e.as_ref() {
        ExprOp::And(l, r) => Some(BExpr::and(
            expr_to_bexpr(l, subs, true).unwrap(),
            expr_to_bexpr(r, subs, true).unwrap(),
        )),
        ExprOp::Or(l, r) => Some(BExpr::or(
            expr_to_bexpr(l, subs, true).unwrap(),
            expr_to_bexpr(r, subs, true).unwrap(),
        )),
        ExprOp::Xor(l, r) => Some(BExpr::xor(
            expr_to_bexpr(l, subs, true).unwrap(),
            expr_to_bexpr(r, subs, true).unwrap(),
        )),
        ExprOp::Not(i) => {
            Some(BExpr::not(expr_to_bexpr(i, subs, true).unwrap()))
        },
        // Otherwise generate a new variable and add the substitution.
        _ => new_sub(),
    }
}

/// Tries to convert an expression into a factor and a [`BExpr`].
fn parse_term<R: BinaryRing>(
    e: &Expr<R>,
    subs: &mut Subs<R>,
    force: bool,
    ring: &R,
) -> Option<(R::Element, BExpr)> {
    if let ExprOp::Mul(l, r) = e.as_ref() {
        if let ExprOp::Const(i) = l.as_ref() {
            return expr_to_bexpr(r, subs, force).map(|u| (i.clone(), u));
        } else if let ExprOp::Const(i) = r.as_ref() {
            return expr_to_bexpr(l, subs, force).map(|u| (i.clone(), u));
        }
    } else if let ExprOp::Const(c) = e.as_ref() {
        return Some((ring.neg(c.clone()), BExpr::Ones));
    }

    expr_to_bexpr(e, subs, force).map(|u| (R::one(), u))
}

/// Converts part of an expression to an [`LBExpr`].
///
/// Returns the [`LBExpr`] and a list of substitutions, such that substituting
/// the expressions in the list into the variables, would give the original
/// expression.
///
/// It will generally try to make the [`LBExpr`] as big as possible.
///
/// If `force` is false, it will return [`None`] if the top-most operation is
/// not something a [`LBExpr`] can represent. Otherwise, it will return a new
/// variable, whose substitution is the original expression.
pub fn expr_to_lbexpr<R: BinaryRing>(
    e: &Expr<R>,
    force: bool,
    r: &R,
) -> Option<(LBExpr<R>, Subs<R>)> {
    let mut lu = LBExpr::zero();
    let mut subs = Subs::new();
    if expr_to_lbexpr_impl(e, &mut lu, &mut subs, false, force, r) {
        Some((lu, subs))
    } else {
        None
    }
}

fn expr_to_lbexpr_impl<R: BinaryRing>(
    e: &Expr<R>,
    lu: &mut LBExpr<R>,
    subs: &mut Subs<R>,
    negate: bool,
    force: bool,
    ring: &R,
) -> bool {
    //if e.strong_count() > 1 {
    //    if force
    //}
    match e.as_ref() {
        ExprOp::Add(l, r) => {
            expr_to_lbexpr_impl(l, lu, subs, negate, true, ring);
            expr_to_lbexpr_impl(r, lu, subs, negate, true, ring);
            true
        },

        ExprOp::Sub(l, r) => {
            expr_to_lbexpr_impl(l, lu, subs, negate, true, ring);
            expr_to_lbexpr_impl(r, lu, subs, !negate, true, ring);
            true
        },

        ExprOp::Neg(i) => {
            // Theoretically we could allow another whole
            // LBExpr in here but hopefully not too important.
            let c = if negate {
                R::one()
            } else {
                ring.negative_one()
            };
            lu.0.push((c, expr_to_bexpr(i, subs, true).unwrap()));
            true
        },

        // Otherwise parse the term from this expression.
        _ => {
            let Some((mut f, u)) = parse_term(e, subs, force, ring) else {
                return false;
            };

            if negate {
                f = ring.neg(f);
            }
            lu.0.push((f, u));
            true
        },
    }
}
