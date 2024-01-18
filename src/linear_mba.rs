//! Linear Mixed Boolean-Arithmetic.

use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{Zero, One, Euclid};

use crate::simplify_boolean::{simplify_from_truth_table, SimplificationConfig};
use crate::valuation::Valuation;
use crate::{diophantine, keep_signed_bits, keep_signed_bits_mut};
use crate::expr::{ExprOp, Expr};
use crate::lattice::{self, AffineLattice};
use crate::uniform_expr::*;
use crate::matrix::*;
use crate::vector::*;
use crate::Symbol;

/// Rewrite a linear combination of uniform expression
/// using a set of uniform expressions modulo `2^bits`.
pub fn rewrite(
    expr: &LUExpr,
    ops: &[LUExpr],
    bits: u32,
    randomize: bool
) -> Option<LUExpr> {
    // Find all variables we have access to.
    // This includes variables in the expression as well as potentially the ops.
    let mut v = std::collections::BTreeSet::new();
    expr.vars_impl(&mut v);
    ops.iter().for_each(|e| e.vars_impl(&mut v));
    let v: Vec<_> = v.into_iter().collect();

    // Solve the system.
    let l = solve_linear_system(expr, ops, &v, bits);

    // Does it have solutions?
    if l.is_empty() {
        return None;
    }

    // Sample a point from the lattice.
    let solution = if randomize {
        l.sample_point(bits)
    } else {
        l.offset
    };

    // Collect the solution into an LUExpr.
    Some(collect_solution(solution.view(), ops, bits))
}

/// Obfuscate any expression using linear MBA.
pub fn obfuscate(e: &mut Expr, bits: u32, cfg: &ObfuscationConfig) {
    // Find all variables we have access to.
    let mut vars = e.vars();
    for i in 0..(cfg.rewrite_vars - vars.len() as isize) {
        vars.push(format!("aux{}", i).into());
    }

    let mut v = Vec::new();
    obfuscate_impl(e, &mut v, &vars, bits, cfg);

    fn rewrite_random(
        e: &LUExpr, vars: &[Symbol], bits: u32, cfg: &ObfuscationConfig
    ) -> LUExpr {
        let mut vars = vars.to_vec();
        for v in e.vars() {
            if !vars.contains(&v) {
                vars.push(v);
            }
        }

        fn try_rewrite(
            e: &LUExpr, vars: &[Symbol], bits: u32, cfg: &ObfuscationConfig
        ) -> Option<LUExpr> {
            let mut ops = Vec::new();
            ops.push(LUExpr::from(UExpr::Ones));
            for _ in 0..cfg.rewrite_expr_count {
                ops.push(LUExpr::from(
                    random_bool_expr(vars, cfg.rewrite_expr_depth)
                ));
            }

            rewrite(e, &ops, bits, true)
        }

        match cfg.rewrite_tries {
            RewriteTries::Infinite => loop {
                if let Some(e) = try_rewrite(e, &vars, bits, cfg) {
                    return e;
                }
            },
            RewriteTries::FiniteFail(t) => {
                for _ in 0..t {
                    if let Some(e) = try_rewrite(e, &vars, bits, cfg) {
                        return e;
                    }
                }
                panic!("Failed to rewrite uniform expression.");
            },
            RewriteTries::FiniteOriginal(t) => {
                for _ in 0..t {
                    if let Some(e) = try_rewrite(e, &vars, bits, cfg) {
                        return e;
                    }
                }
                e.clone()
            },
        }
    }

    fn obfuscate_impl(
        er: &mut Expr,
        visited: &mut Vec<*const ExprOp>,
        vars: &[Symbol],
        bits: u32,
        cfg: &ObfuscationConfig
    ) {
        // Check if we have already visited this expression.
        let ptr = er.as_ptr();
        if er.strong_count() > 1 {
            if visited.contains(&ptr) {
                return;
            }
            visited.push(ptr);
        }

        let e = unsafe { &mut *(ptr as *mut ExprOp) };

        // Try to find the largest subexpression that is
        // linear MBA and obfuscate it on its own.
        if let Some((lu, mut subs)) = expr_to_luexpr(er, false) {
            *e = rewrite_random(&lu, vars, bits, cfg).to_expr();
            for (var, sub) in &mut subs.0 {
                obfuscate_impl(sub, visited, vars, bits, cfg);
                er.substitute(*var, sub);
            }
            return;
        }

        // If the expression isn't linear MBA, recurse on the subexpressions.
        match e {
            ExprOp::Mul(l, r) | ExprOp::Div(l, r) | ExprOp::Shl(l, r)
            | ExprOp::Shr(l, r) | ExprOp::Sar(l, r) => {
                obfuscate_impl(l, visited, vars, bits, cfg);
                obfuscate_impl(r, visited, vars, bits, cfg);
            }
            _ => panic!("Expression should be linear MBA, \
                but expr_to_luexpr failed ({:?}).", e),
        }
    }
}

pub fn deobfuscate(e: &mut Expr, bits: u32) {
    let cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::LeastComplexTerms,
        boolean: true,
    };

    let mut v = Vec::new();
    deobfuscate_impl(e, &mut v, bits, &cfg);
    e.simplify();

    fn deobfuscate_impl(
        er: &mut Expr,
        visited: &mut Vec<*const ExprOp>,
        bits: u32,
        cfg: &DeobfuscationConfig
    ) {
        // Check if we have already visited this expression.
        let ptr = er.as_ptr();
        if er.strong_count() > 1 {
            if visited.contains(&ptr) {
                return;
            }
            visited.push(ptr);
        }

        let e = unsafe { &mut *(ptr as *mut ExprOp) };

        // Try to find the largest subexpression that is
        // linear MBA and obfuscate it on its own.
        if let Some((lu, mut subs)) = expr_to_luexpr(er, false) {
            let lu = deobfuscate_luexpr(lu, bits, cfg);
            *e = lu.to_expr();

            // Deobfuscate all the subexpressions
            // and substitute them into the expression.
            for (var, sub) in &mut subs.0 {
                deobfuscate_impl(sub, visited, bits, cfg);
                er.substitute(*var, sub);
            }
            return;
        }

        // If the expression isn't linear MBA, recurse on the subexpressions.
        match e {
            ExprOp::Mul(l, r) | ExprOp::Div(l, r) | ExprOp::Shl(l, r)
            | ExprOp::Shr(l, r) | ExprOp::Sar(l, r) => {
                deobfuscate_impl(l, visited, bits, cfg);
                deobfuscate_impl(r, visited, bits, cfg);
            }
            _ => panic!("Expression should be linear MBA, \
                but expr_to_luexpr failed ({:?}).", e),
        }
    }
}

pub struct DeobfuscationConfig {
    /// The algorithm to use for finding a small solution.
    pub alg: SolutionAlgorithm,

    /// Detect whether the expression is purely boolean
    /// and find a simplified expression from the truth table.
    pub boolean: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionAlgorithm {
    /// Just use the solution from the linear system solver.
    /// This should not really be used as
    /// [SolutionAlgorithm::LeastComplexTerms] is almost as
    /// fast and much better.
    Fast,

    /// Uses the least complex rewrite operations possible.
    /// This goes from most complex to least complex rewrite operation
    /// and tries to remove the rewrite operation from the solution.
    LeastComplexTerms,

    /// Not recommended as it is slow and can result in
    /// linear combinations with lots of terms all of which
    /// have small coefficients.
    /// Finds the shortest solution on the affine lattice
    /// of solutions to the linear system where the entries
    /// are weighted by the complexity of their corresponding
    /// rewrite operations.
    ShortVector,
}

/// Deobfuscate a linear MBA expression.
pub fn deobfuscate_luexpr(e: LUExpr, bits: u32, cfg: &DeobfuscationConfig) -> LUExpr {
    use SolutionAlgorithm::*;

    // Get all the variables.
    let vars = e.vars();

    // This could be optimized.
    // We compute the values of the expression twice.
    // Once here and once in `solve_linear_system`.
    if cfg.boolean {
        // Allocate a valuation.
        let mut val = Valuation::zero();

        let entries = 1 << vars.len();
        let mut v = Vec::with_capacity(entries);
        for i in 0..entries {
            // Initialize the valuation.
            for (j, c) in vars.iter().enumerate() {
                val.set_value(*c, -BigInt::from((i >> j) & 1));
            }

            // Evaluate the expression.
            let r = keep_signed_bits(&e.eval(&mut val, bits), bits);
            let r = if r.is_zero() {
                false
            } else if r.sign() == num_bigint::Sign::Minus
                && r.magnitude().is_one() {
                true
            } else {
                break
            };

            // Write the desired result into the vector.
            v.push(r);
        }

        if v.len() == entries {
            let cfg = SimplificationConfig::default();
            let u = simplify_from_truth_table(&vars, v.iter().cloned(), &cfg);
            return u.into();
        }
    }

    //
    // Insert some default operations.
    //
    let mut ops = Vec::new();

    // Constants
    ops.push(UExpr::Ones);

    // Variables
    for v in &vars {
        ops.push(UExpr::Var(*v));
    }

    // Binary operations (and, or, xor) for all variable pairs.
    for v in &vars {
        for w in &vars {
            if std::ptr::eq(v, w) {
                break;
            }

            ops.push(UExpr::and(UExpr::Var(*w), UExpr::Var(*v)));
            ops.push(UExpr::or(UExpr::Var(*w), UExpr::Var(*v)));
            ops.push(UExpr::xor(UExpr::Var(*w), UExpr::Var(*v)));
        }
    }

    // Insert the original operations so we always have a solution.
    for u in &e.0 {
        if !ops.contains(&u.1) {
            ops.push(u.1.clone());
        }
    }

    if cfg.alg == LeastComplexTerms {
        // Sort by descending complexity.
        ops.sort_by_cached_key(|e| -(e.complexity() as i32));
    }

    // Convert them to LUExprs.
    // We could template this if performance is actually an issue.
    let ops: Vec<_> = ops.iter().map(|e| LUExpr::from(e.clone())).collect();

    // Solve the system.
    let mut l = solve_linear_system(&e, &ops, &vars, bits);

    // If I did everything correctly,
    // this system should always have a solution.
    assert!(l.offset.dim() == ops.len());

    if cfg.alg == LeastComplexTerms {
        for i in 0..l.offset.dim() {
            assert!(!l.lattice.basis[i][i].is_zero());
            let q = l.offset[i].clone().div_euclid(&l.lattice.basis[i][i]);
            if !q.is_zero() {
                l.offset -= &(&l.lattice.basis[i] * &q);
            }
        }
    }

    // If this should be fast, just use the particular solution
    // from the solver. In practice, this seems to work well.
    //if fast {
    if matches!(cfg.alg, Fast | LeastComplexTerms) {
        return collect_solution(l.offset.view(), &ops, bits);
    }

    //
    // Otherwise we actually find the smallest solution in the lattice.
    // Since the LLL/CVP algorithms expect the L^2 norm,
    // but we want different coordinates (corresponding to different
    // rewrite operations) to have different weights,
    // we scale the lattice basis and offset accordingly.
    //

    // Compute the complexity measure for each coordinate.
    let complexity: Vec<_> = ops.iter()
        .map(|e| e.complexity(bits))
        .collect();

    // Scale the lattice basis...
    for v in l.lattice.basis.rows_mut() {
        for (e, c) in v.iter_mut().zip(&complexity) {
            *e *= *c;
        }
    }

    // ... and offset.
    for (e, c) in l.offset.iter_mut().zip(&complexity) {
        *e *= *c;
    }

    // Run LLL because it improves the speed of the CVP algorithm.
    //l.lattice.lll(&0.9921875, lattice::F64);
    l.lattice.lll(&0.75, lattice::F64);

    // Compute the shortest vector in the affine lattice via CVP
    // (used internally).
    let mut solution = l.svp(None, lattice::F64).unwrap();

    // Divide the solution by the complexity.
    for (e, c) in solution.iter_mut().zip(&complexity) {
        debug_assert!(e.is_multiple_of(&BigInt::from(*c)));
        *e /= *c;
    }

    // Collect it into an LUExpr.
    collect_solution(solution.view(), &ops, bits)
}

/// Creates the linear system corresponding to the target expression
/// and rewrite operations.
fn solve_linear_system(
    expr: &LUExpr,
    ops: &[LUExpr],
    vars: &[Symbol],
    bits: u32
) -> AffineLattice {
    // If you want more than 64 vars, get a lot of RAM
    // and change the iterators to use Integer instead of usize.
    assert!(vars.len() < 64, "More than 63 variables are currently \
        not supported (You wouldn't be able to run this anyways).");

    // Allocate a valuation.
    let mut val = Valuation::zero();

    let rows = 1 << vars.len();
    let cols = ops.len();

    let mut a = Matrix::zero(rows, cols);
    let mut b = OwnedVector::zero(rows);

    // Build up the matrix.
    for i in 0..rows {
        let row = a.row_mut(i);

        // Initialize the valuation.
        for (j, c) in vars.iter().enumerate() {
            val.set_value(*c, -BigInt::from((i >> j) & 1));
        }

        // Write the values of the operations into this row of the matrix.
        for (j, e) in ops.iter().enumerate() {
            row[j] = keep_signed_bits(&e.eval(&mut val, bits), bits);
        }

        // Write the desired result into the vector.
        b[i] = keep_signed_bits(&expr.eval(&mut val, bits), bits);
    }

    diophantine::solve_modular(a.view(), b.view(), &(BigInt::one() << bits))
}

/// Converts a solution to the linear system and the operations
/// to a single `LUExpr` with no duplicate `UExpr`s.
fn collect_solution(
    solution: &IVectorView,
    ops: &[LUExpr],
    bits: u32
) -> LUExpr {
    // Put it in a LUExpr.
    let mut l = LUExpr::zero();
    for (c, o) in solution.iter().zip(ops.iter()) {
        for (d, e) in &o.0 {
            // Is the UExpr already in the linear combination?
            match l.0.iter_mut().find(|(_, f)| f == e) {
                Some((f, _)) => *f += c * d,
                None => l.0.push((c * d, e.clone())),
            }
        }
    }

    for (i, _) in &mut l.0 {
        keep_signed_bits_mut(i, bits);
    }

    l.remove_zero_terms();
    l
}

/// A list of substitutions.
/// See e.g. [expr_to_uexpr] on what this is used for.
#[derive(Default)]
struct Subs(pub Vec<(Symbol, Expr)>);

impl Subs {
    /// Creates a new empty substitution list.
    pub fn new() -> Self {
        Self(Vec::new())
    }

    /// Adds a new substitution.
    pub fn add(&mut self, expr: Expr) -> Symbol {
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

/// Converts part of an expression to a UExpr,
/// such that if you substituted the Exprs in `subs` for the variables,
/// you would get the original Expr.
/// It will generally try to make the LUExpr as big as possible.
/// If `force` is false, it will return [None] if the top-most operation
/// is not a [UExpr] operation.
/// Otherwise, it will return a [UExpr::Var] whose substitution
/// is the original expression.
fn expr_to_uexpr(
    e: &Expr, subs: &mut Subs, force: bool
) -> Option<UExpr> {
    // New substitution variable.
    let mut new_sub = || force.then(|| UExpr::Var(subs.add(e.clone())));

    if let ExprOp::Var(v) = e.as_ref() {
        return Some(UExpr::Var(*v));
    }

    // We don't try, when the expression is shared.
    if e.strong_count() > 1 {
        return new_sub();
    }

    match e.as_ref() {
        ExprOp::And(l, r) => Some(UExpr::and(
            expr_to_uexpr(l, subs, true).unwrap(),
            expr_to_uexpr(r, subs, true).unwrap()
        )),
        ExprOp::Or(l, r) => Some(UExpr::or(
            expr_to_uexpr(l, subs, true).unwrap(),
            expr_to_uexpr(r, subs, true).unwrap()
        )),
        ExprOp::Xor(l, r) => Some(UExpr::xor(
            expr_to_uexpr(l, subs, true).unwrap(),
            expr_to_uexpr(r, subs, true).unwrap()
        )),
        ExprOp::Not(i) => Some(UExpr::not(
            expr_to_uexpr(i, subs, true).unwrap()
        )),
        // Otherwise generate a new variable and add the substitution.
        _ => new_sub(),
    }
}

/// Tries to convert an expression into a factor and a UExpr.
fn parse_term(
    e: &Expr, subs: &mut Subs, force: bool
) -> Option<(BigInt, UExpr)> {
    if let ExprOp::Mul(l, r) = e.as_ref() {
        if let ExprOp::Const(i) = l.as_ref() {
            return expr_to_uexpr(r, subs, force).map(|u| (i.clone(), u));
        } else if let ExprOp::Const(i) = r.as_ref() {
            return expr_to_uexpr(l, subs, force).map(|u| (i.clone(), u));
        }
    } else if let ExprOp::Const(c) = e.as_ref() {
        return Some((-c.clone(), UExpr::Ones));
    }

    expr_to_uexpr(e, subs, force).map(|u| (One::one(), u))
}

/// Converts part of an expression to an [LUExpr].
/// Returns the [LUExpr] and a list of substitutions,
/// such that substituting the expressions in the list
/// into the variables, would give the original expression.
/// It will generally try to make the LUExpr as big as possible.
/// If `force` is false, it will return [None] if the top-most operation
/// is not something a [LUExpr] can represent (e.g. [Expr::Div]).
/// Otherwise, it will return a new variable, whose substitution
/// is the original expression.
fn expr_to_luexpr(
    e: &Expr, force: bool
) -> Option<(LUExpr, Subs)> {
    let mut lu = LUExpr::zero();
    let mut subs = Subs::new();
    if expr_to_luexpr_impl(e, &mut lu, &mut subs, false, force) {
        Some((lu, subs))
    } else {
        None
    }
}

fn expr_to_luexpr_impl(
    e: &Expr,
    lu: &mut LUExpr,
    subs: &mut Subs,
    negate: bool,
    force: bool,
) -> bool {
    //if e.strong_count() > 1 {
    //    if force
    //}
    match e.as_ref() {
        ExprOp::Add(l, r) => {
            expr_to_luexpr_impl(l, lu, subs, negate, true);
            expr_to_luexpr_impl(r, lu, subs, negate, true);
            true
        },

        ExprOp::Sub(l, r) => {
            expr_to_luexpr_impl(l, lu, subs, negate, true);
            expr_to_luexpr_impl(r, lu, subs, !negate, true);
            true
        },

        ExprOp::Neg(i) => {
            // Theoretically we could allow another whole
            // LUExpr in here but hopefully not too important.
            let c = if negate { 1 } else { -1 };
            lu.0.push((c.into(), expr_to_uexpr(i, subs, true).unwrap()));
            true
        },

        // Otherwise parse the term from this expression.
        _ => {
            let Some((mut f, u)) = parse_term(e, subs, force) else {
                return false
            };

            if negate {
                f = -f;
            }
            lu.0.push((f, u));
            true
        }
    }
}

/// Generates a random boolean expression.
/// It would be very desirable to make this smarter.
/// Currently it generates a lot of non-sense expressions,
/// which simplify to zero or one easily.
fn random_bool_expr(vars: &[Symbol], max_depth: usize) -> UExpr {
    assert!(!vars.is_empty(), "There needs to be \
        at least one variable for the random expression.");

    let rand_var = || UExpr::Var(
        vars[rand::random::<usize>() % vars.len()]
    );

    if max_depth == 0 {
        return rand_var();
    }

    // Generate one of the four variants uniformly at random.
    let d = max_depth - 1;
    match rand::random::<usize>() % 5 {
        0 => rand_var(),
        1 => UExpr::not(random_bool_expr(vars, d)),
        2 => UExpr::and(random_bool_expr(vars, d), random_bool_expr(vars, d)),
        3 => UExpr::or(random_bool_expr(vars, d), random_bool_expr(vars, d)),
        4 => UExpr::xor(random_bool_expr(vars, d), random_bool_expr(vars, d)),
        _ => unreachable!(),
    }
}


/// More obscure settings for internals of the obfuscation.
/// Use `ObfuscationConfig::default()` to get a reasonable default.
pub struct ObfuscationConfig {
    /// The number of variables that appear in the obfuscated expression.
    /// The default is 4.
    pub rewrite_vars: isize,

    /// The maximum depth of the rewrite operations.
    /// The default is 3.
    pub rewrite_expr_depth: usize,

    /// The number of rewrite operations.
    /// The default is 24.
    pub rewrite_expr_count: usize,

    /// We randomly generate rewrite operations,
    /// which means that it might not be possible to
    /// rewrite and LUExpr with them.
    /// In particular, if the other settings are too weak,
    /// it might be impossible for any rewrite operations
    /// to rewrite the expression.
    /// The default is `FiniteFail(128)`,
    /// which will panic after 128 attempts.
    pub rewrite_tries: RewriteTries,
}

impl Default for ObfuscationConfig {
    /// Returns a reasonable default configuration.
    /// The the documentation of the members for the default values.
    fn default() -> Self {
        Self {
            rewrite_vars: 4,
            rewrite_expr_depth: 3,
            rewrite_expr_count: 24,
            rewrite_tries: RewriteTries::FiniteFail(128),
        }
    }
}

/// Used internally in [ObfuscationConfig].
pub enum RewriteTries {
    /// Try forever. (Not a good idea).
    /// It is possible for an LUExpr to be too complex
    /// for the current settings, so that this will never terminate.
    Infinite,

    /// Try for this many times and panic if the rewrite fails.
    FiniteFail(usize),

    /// Try for this many times and return the
    /// original (unobfuscated) expression if the rewrite fails.
    FiniteOriginal(usize),
}

#[test]
fn linear_obfuscate_test() {
    let cfg = ObfuscationConfig::default();
    let bits = 8;
    let e = Expr::from_string("x + y * x").unwrap();
    let mut o = e.deep_copy();
    obfuscate(&mut o, bits, &cfg);
    for _ in 0..100 {
        let mut v = Valuation::random(bits);
        assert_eq!(e.eval(&mut v, bits), o.eval(&mut v, bits));
    }
}

#[test]
fn deobfuscate_linear_test() {
    let cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::LeastComplexTerms,
        boolean: true,
    };
    let x = Symbol::from("x");
    let y = Symbol::from("y");
    let e = UExpr::xor(UExpr::not(UExpr::var(x)), UExpr::var(y));
    let d = deobfuscate_luexpr(e.clone().into(), 4, &cfg);
    println!("{d}");

    let mut v = Valuation::zero();
    for xv in 0..15 {
        *v.value(x) = xv.into();
        for yv in 0..15 {
            *v.value(y) = yv.into();
            assert_eq!(keep_signed_bits(&e.eval(&mut v), 4),
                keep_signed_bits(&d.eval(&mut v, 4), 4));
        }
    }
}

#[test]
fn deobfuscate_test() {
    let bits = 8;
    //let e = Expr::from_string("x + y * x").unwrap();
    let e = Expr::from_string("179 * x - 6 + y * x * 4").unwrap();

    // Create an example expression.
    println!("Original:\n{e}");

    // Obfuscate it.
    let mut o = e.deep_copy();
    let cfg = ObfuscationConfig {
        //rewrite_vars: 8,
        rewrite_expr_depth: 2,
        rewrite_expr_count: 48,
        ..ObfuscationConfig::default()
    };
    obfuscate(&mut o, bits, &cfg);
    println!("Obfuscated:\n{o}");

    for _ in 0..100 {
        let mut v = Valuation::random(bits);
        assert_eq!(e.eval(&mut v, bits), o.eval(&mut v, bits));
    }

    // Deobfuscate it.
    deobfuscate(&mut o, bits);
    println!("Deobfuscated:\n{o}");

    for _ in 0..100 {
        let mut v = Valuation::random(bits);
        assert_eq!(e.eval(&mut v, bits), o.eval(&mut v, bits));
    }
}