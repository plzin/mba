//! Obfuscation of an arbitrary [`Expr`].
//!
//! The current algorithm is not at all meant to be good or interesting. All it
//! does is find the largest possible subexpressions that are linear MBA and
//! rewrites them using randomly generated rewrite operations.

use rand::Rng;

use super::{rewrite::rewrite, subexpression::expr_to_lbexpr};
use crate::{
    Symbol,
    bitwise_expr::{BExpr, LBExpr},
    expr::{Expr, ExprOp},
    rings::BinaryRing,
};

/// More obscure settings for internals of the obfuscation.
/// Use `ObfuscationConfig::default()` to get a reasonable default.
pub struct ObfuscationConfig {
    /// The number of additional variables that are added in the obfuscated
    /// expression.
    ///
    /// The default is 2.
    pub auxiliary_vars: usize,

    /// The maximum depth of the rewrite operations.
    ///
    /// The default is 3.
    pub rewrite_expr_depth: usize,

    /// The number of rewrite operations.
    ///
    /// The default is 24.
    pub rewrite_expr_count: usize,

    /// We randomly generate rewrite operations, which means that it might not
    /// be possible to rewrite and [`LBExpr`] with them. In particular, if the
    /// other settings are too weak, it might be impossible for any rewrite
    /// operations to rewrite the expression.
    ///
    /// The default is `FiniteFail(128)`, which will panic after 128 attempts.
    pub rewrite_tries: RewriteTries,
}

impl Default for ObfuscationConfig {
    /// Returns a reasonable default configuration.
    /// The the documentation of the members for the default values.
    fn default() -> Self {
        Self {
            auxiliary_vars: 2,
            rewrite_expr_depth: 3,
            rewrite_expr_count: 24,
            rewrite_tries: RewriteTries::FiniteFail(128),
        }
    }
}

/// Used internally in [ObfuscationConfig].
pub enum RewriteTries {
    /// Try forever. (Not a good idea).
    /// It is possible for an [`LBExpr`] to be too complex
    /// for the current settings, so that this will never terminate.
    Infinite,

    /// Try for this many times and panic if the rewrite fails.
    FiniteFail(usize),

    /// Try for this many times and return the
    /// original (unobfuscated) expression if the rewrite fails.
    FiniteOriginal(usize),
}

/// Obfuscate any expression using linear MBA.
pub fn obfuscate<R: BinaryRing, Rand: Rng>(
    e: &mut Expr<R>,
    cfg: &ObfuscationConfig,
    rng: &mut Rand,
    ring: &R,
) {
    // Find all variables we have access to.
    let mut vars = e.vars();
    for i in 0..cfg.auxiliary_vars {
        vars.push(format!("aux{i}").into());
    }

    let mut v = Vec::new();
    obfuscate_impl(e, &mut v, &vars, cfg, rng, ring);

    fn rewrite_random<R: BinaryRing, Rand: Rng>(
        e: &LBExpr<R>,
        vars: &[Symbol],
        cfg: &ObfuscationConfig,
        rng: &mut Rand,
        ring: &R,
    ) -> LBExpr<R> {
        let mut vars = vars.to_vec();
        for v in e.vars() {
            if !vars.contains(&v) {
                vars.push(v);
            }
        }

        fn try_rewrite<R: BinaryRing, Rand: Rng>(
            e: &LBExpr<R>,
            vars: &[Symbol],
            cfg: &ObfuscationConfig,
            rng: &mut Rand,
            ring: &R,
        ) -> Option<LBExpr<R>> {
            let mut ops = Vec::new();
            ops.push(LBExpr::from(BExpr::Ones));
            for _ in 0..cfg.rewrite_expr_count {
                ops.push(LBExpr::from(random_bool_expr(
                    vars,
                    cfg.rewrite_expr_depth,
                    rng,
                )));
            }

            rewrite(e, &ops, Some(rng), ring)
        }

        match cfg.rewrite_tries {
            RewriteTries::Infinite => loop {
                if let Some(e) = try_rewrite(e, &vars, cfg, rng, ring) {
                    return e;
                }
            },
            RewriteTries::FiniteFail(t) => {
                for _ in 0..t {
                    if let Some(e) = try_rewrite(e, &vars, cfg, rng, ring) {
                        return e;
                    }
                }
                panic!("Failed to rewrite bitwise expression.");
            },
            RewriteTries::FiniteOriginal(t) => {
                for _ in 0..t {
                    if let Some(e) = try_rewrite(e, &vars, cfg, rng, ring) {
                        return e;
                    }
                }
                e.clone()
            },
        }
    }

    fn obfuscate_impl<R: BinaryRing, Rand: Rng>(
        er: &mut Expr<R>,
        visited: &mut Vec<*const ExprOp<R>>,
        vars: &[Symbol],
        cfg: &ObfuscationConfig,
        rng: &mut Rand,
        ring: &R,
    ) {
        // Check if we have already visited this expression.
        let ptr = er.as_ptr();
        if er.strong_count() > 1 {
            if visited.contains(&ptr) {
                return;
            }
            visited.push(ptr);
        }

        let e = unsafe { &mut *(ptr as *mut ExprOp<R>) };

        // Try to find the largest subexpression that is
        // linear MBA and obfuscate it on its own.
        if let Some((lu, mut subs)) = expr_to_lbexpr(er, false, ring) {
            *e = rewrite_random(&lu, vars, cfg, rng, ring).to_expr(ring);
            for (var, sub) in &mut subs.0 {
                obfuscate_impl(sub, visited, vars, cfg, rng, ring);
                er.substitute(*var, sub);
            }
            return;
        }

        // If the expression isn't linear MBA, recurse on the subexpressions.
        match e {
            ExprOp::Mul(l, r) => {
                obfuscate_impl(l, visited, vars, cfg, rng, ring);
                obfuscate_impl(r, visited, vars, cfg, rng, ring);
            },
            _ => panic!(
                "Expression should be linear MBA, but expr_to_lbexpr failed \
                 ({e:?})."
            ),
        }
    }
}

/// Generates a random boolean expression.
/// It would be very desirable to make this smarter.
/// Currently it generates a lot of non-sense expressions,
/// which simplify to zero or one easily.
fn random_bool_expr<Rand: Rng>(
    vars: &[Symbol],
    max_depth: usize,
    rng: &mut Rand,
) -> BExpr {
    assert!(
        !vars.is_empty(),
        "There needs to be at least one variable for the random expression."
    );
    let num_vars = vars.len() as u32;
    assert_eq!(
        num_vars as usize,
        vars.len(),
        "Not more than `u32::MAX` variables allowed."
    );

    if max_depth == 0 {
        return BExpr::Var(vars[(rng.random::<u32>() % num_vars) as usize]);
    }

    // Generate one of the four variants uniformly at random.
    let d = max_depth - 1;
    match rng.random::<u8>() % 5 {
        0 => BExpr::Var(vars[(rng.random::<u32>() % num_vars) as usize]),
        1 => BExpr::not(random_bool_expr(vars, d, rng)),
        2 => BExpr::and(
            random_bool_expr(vars, d, rng),
            random_bool_expr(vars, d, rng),
        ),
        3 => BExpr::or(
            random_bool_expr(vars, d, rng),
            random_bool_expr(vars, d, rng),
        ),
        4 => BExpr::xor(
            random_bool_expr(vars, d, rng),
            random_bool_expr(vars, d, rng),
        ),
        _ => unreachable!(),
    }
}

#[test]
fn linear_obfuscate_test() {
    use rand::{SeedableRng as _, rngs::StdRng};

    use crate::{rings::U8, valuation::Valuation};

    let mut rng = StdRng::seed_from_u64(0);
    let cfg = ObfuscationConfig::default();
    let r = &U8;
    let e = Expr::<U8>::from_string("x + y * x".to_owned(), r).unwrap();
    let mut o = e.deep_copy();
    obfuscate(&mut o, &cfg, &mut rng, r);
    for _ in 0..100 {
        let mut v = Valuation::random_seeded(0);
        assert_eq!(e.eval(&mut v, r), o.eval(&mut v, r));
    }
}
