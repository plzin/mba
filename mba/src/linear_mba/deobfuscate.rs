//! Deobfuscation of an [`Expr`].
//!
//! Note that this is very limited. All it does is find the largest possible
//! linear MBA subexpressions and deobfuscates them using a set of common
//! rewrite operations.

use num_bigint::BigInt;
use num_integer::Integer;

use crate::rings::{BinaryRing, Ring as _, RingElement as _, F64, Z};
use crate::expr::{Expr, ExprOp};
use crate::bitwise_expr::{BExpr, LBExpr};
use crate::simplify_boolean::{simplify_from_truth_table, SimplificationConfig};
use crate::valuation::Valuation;
use crate::matrix::Matrix;
use crate::lattice::{AffineLattice, Lattice};
use super::{rewrite::{solve_linear_system, collect_solution}, subexpression::expr_to_lbexpr};

/// Configuration for deobfuscation.
pub struct DeobfuscationConfig {
    /// The algorithm to use for finding a small solution.
    pub alg: SolutionAlgorithm,

    /// Detect whether the expression is purely boolean
    /// and find a simplified expression from the truth table.
    pub boolean: bool,
}

/// The algorithm to use for finding a small solution.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolutionAlgorithm {
    /// Just use the solution from the linear system solver.
    /// This should not really be used as
    /// [`SolutionAlgorithm::LeastComplexTerms`] is almost as
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


pub fn deobfuscate<R: BinaryRing>(
    e: &mut Expr<R>,
    cfg: &DeobfuscationConfig,
    ring: &R,
) {
    let mut v = Vec::new();
    deobfuscate_impl(e, &mut v, cfg, ring);
    e.simplify();

    fn deobfuscate_impl<R: BinaryRing>(
        er: &mut Expr<R>,
        visited: &mut Vec<*const ExprOp<R>>,
        cfg: &DeobfuscationConfig,
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
            let lu = deobfuscate_lbexpr(lu, cfg, ring);
            *e = lu.to_expr(ring);

            // Deobfuscate all the subexpressions
            // and substitute them into the expression.
            for (var, sub) in &mut subs.0 {
                deobfuscate_impl(sub, visited, cfg, ring);
                er.substitute(*var, sub);
            }
            return;
        }

        // If the expression isn't linear MBA, recurse on the subexpressions.
        match e {
            ExprOp::Mul(l, r) => {
                deobfuscate_impl(l, visited, cfg, ring);
                deobfuscate_impl(r, visited, cfg, ring);
            }
            _ => panic!("Expression should be linear MBA, \
                but expr_to_lbexpr failed ({e:?})."),
        }
    }
}

/// Deobfuscate a linear MBA expression.
pub fn deobfuscate_lbexpr<R: BinaryRing>(
    e: LBExpr<R>,
    cfg: &DeobfuscationConfig,
    ring: &R,
) -> LBExpr<R> {
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
                val.set_value(*c, ring.neg(ring.element_from_usize((i >> j) & 1)));
            }

            // Evaluate the expression.
            let r = e.eval(&mut val, ring);
            let r = if r.is_zero() {
                false
            } else if ring.neg(r).is_one() {
                true
            } else {
                break;
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
    ops.push(BExpr::Ones);

    // Variables
    for v in &vars {
        ops.push(BExpr::Var(*v));
    }

    // Binary operations (and, or, xor) for all variable pairs.
    for v in &vars {
        for w in &vars {
            if std::ptr::eq(v, w) {
                break;
            }

            ops.push(BExpr::and(BExpr::Var(*w), BExpr::Var(*v)));
            ops.push(BExpr::or(BExpr::Var(*w), BExpr::Var(*v)));
            ops.push(BExpr::xor(BExpr::Var(*w), BExpr::Var(*v)));
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
        ops.sort_by_cached_key(|e| std::cmp::Reverse(e.complexity()));
    }

    // Convert them to LBExprs.
    // We could template this if performance is actually an issue.
    let ops: Vec<_> = ops.into_iter().map(LBExpr::from).collect();

    // Solve the system.
    let mut l = solve_linear_system(&e, &ops, &vars, ring);

    // If I did everything correctly,
    // this system should always have a solution.
    assert_eq!(l.offset.dim(), ops.len());

    if cfg.alg == LeastComplexTerms {
        // It might be possible to make this more efficient by not doing the
        // Hermite normal form using `BigInt`s, but you can't just use the
        // ring elements.
        let lc = l.canonicalize(ring);
        l.offset = lc.offset.transform(|e| ring.element_from_bigint(e));
        l.offset.reduce(&l.lattice.basis, ring);
    }

    // If this should be fast, just use the particular solution
    // from the solver. In practice, this seems to work well.
    if matches!(cfg.alg, Fast | LeastComplexTerms) {
        return collect_solution(l.offset.view(), &ops, ring);
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
        .map(|e| e.complexity() as usize)
        .collect();

    // We need to create an integer lattice.
    let mut generating_set = Matrix::zero(
        l.lattice.rank() + l.lattice.ambient_dim(),
        l.lattice.ambient_dim(),
    );

    for i in 0..l.lattice.ambient_dim() {
        generating_set[(l.lattice.rank() + i, i)] = Z::one() << ring.bits();
    }

    let mut offset = l.offset.transform::<Z, _>(
        |i| BigInt::from(R::to_representative(i))
    );

    // Scale the lattice basis...
    for (old_row, new_row) in l.lattice.basis.rows()
        .zip(generating_set.rows_mut())
    {
        for ((old, new), &c) in old_row.iter()
            .zip(new_row)
            .zip(&complexity)
        {
            *new = BigInt::from(R::to_representative(old) * c);
        }
    }

    // ... and offset.
    for (e, &c) in offset.iter_mut().zip(&complexity) {
        *e *= c;
    }

    let mut int_lattice = AffineLattice {
        lattice: Lattice::from_generating_set(generating_set, &Z),
        offset,
    };

    // Run LLL because it improves the speed of the CVP algorithm.
    //int_lattice.lattice.lll(&0.9921875, &F64, &Z);
    int_lattice.lattice.lll(&0.75, &F64, &Z);

    // Compute the shortest vector in the affine lattice via CVP
    // (used internally).
    let mut solution = int_lattice.svp(None, &F64, &Z).unwrap();

    // Divide the solution by the complexity.
    for (e, c) in solution.iter_mut().zip(&complexity) {
        debug_assert!(e.is_multiple_of(&BigInt::from(*c)));
        *e /= *c;
    }

    let s = solution.transform(|n| ring.element_from_bigint(n));

    // Collect it into an LBExpr.
    collect_solution(s.view(), &ops, ring)
}

#[test]
fn deobfuscate_linear_test() {
    use crate::Symbol;
    let cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::LeastComplexTerms,
        boolean: true,
    };
    let r = crate::rings::BinaryBigInt::new(4);
    let x = Symbol::from("x");
    let y = Symbol::from("y");
    let e = BExpr::xor(BExpr::not(BExpr::var(x)), BExpr::var(y));
    let d = deobfuscate_lbexpr(e.clone().into(), &cfg, &r);
    println!("{d}");

    let mut v = Valuation::zero();
    for xv in 0usize..15 {
        *v.value(x, &r) = xv.into();
        for yv in 0usize..15 {
            *v.value(y, &r) = yv.into();
            assert_eq!(e.eval(&mut v, &r), d.eval(&mut v, &r));
        }
    }
}

#[test]
fn deobfuscate_test() {
    use rand::{SeedableRng as _, rngs::StdRng};
    use crate::rings::U8;
    //use crate::formatter::Formatter;
    use crate::linear_mba::obfuscate::{obfuscate, ObfuscationConfig};

    let mut rng = StdRng::seed_from_u64(0);
    let r = &U8;

    // Create an example expression.
    //let e = Expr::from_string("x + y * x").unwrap();
    let e = Expr::from_string("179 * x - 6 + y * x * 4".to_owned(), &U8)
        .unwrap();

    //println!("Original:\n{}", e.display(Formatter::Rust, 0, r));

    // Obfuscate it.
    let mut o = e.deep_copy();
    let obf_cfg = ObfuscationConfig {
        //rewrite_vars: 8,
        rewrite_expr_depth: 2,
        rewrite_expr_count: 48,
        ..ObfuscationConfig::default()
    };
    obfuscate(&mut o, &obf_cfg, &mut rng, r);
    //println!("Obfuscated:\n{}", o.display(Formatter::Rust, 0, r));

    for _ in 0..100 {
        let mut v = Valuation::random_seeded(0);
        assert_eq!(e.eval(&mut v, r), o.eval(&mut v, r));
    }

    // Deobfuscate it.
    let deobf_cfg = DeobfuscationConfig {
        alg: SolutionAlgorithm::LeastComplexTerms,
        boolean: true,
    };
    deobfuscate(&mut o, &deobf_cfg, r);
    //println!("Deobfuscated:\n{}", o.display(Formatter::Rust, 0, r));

    for _ in 0..100 {
        let mut v = Valuation::random_seeded(0);
        assert_eq!(e.eval(&mut v, r), o.eval(&mut v, r));
    }
}