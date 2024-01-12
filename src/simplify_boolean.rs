//! Simplify boolean expressions using egg.
//! There are some obvious improvements you could make:
//! - Given a truth table, find a minimal expression
//!  using e.g. the Quine-McCluskey algorithm, before handing
//!  it to egg.
//! - Use a custom RewriteSchedular that is aware of
//!  the cost function, i.e. tries cost reducing rewrites
//!  more aggressively.
//! - Maybe a more refined cost function.
//! - More refined rewrite rules.

use std::cell::RefCell;
use egg::*;
use crate::{uniform_expr::UExpr, valuation::Valuation};

/// The way we simplify is by running equality saturation
/// on the expression. The egraph for the expression would be
/// prohibitively large, so that full equality saturation
/// would take too long. egg has a lot of parameters that
/// change the behavior of when it stops.
/// See [https://docs.rs/egg/0.9.5/egg/struct.Runner.html](here)
/// for the documentation of most of these fields.
/// When equality saturation stops without the actual
/// saturation being complete, we start again with the
/// current best expression. If we do not find a better
/// expression in the next `retry_same_cost` iterations,
/// we will stop. We will do this for `retry_limit` times
/// and return the best expression we found.
pub struct SimplificationConfig {
    /// The maximum number of times to start again
    /// with the current best expression.
    /// This is not an egg parameter.
    /// The default is 10, though you might want to
    /// adjust this, depending on the number of variables
    /// and other settings.
    retry_limit: usize,
    /// The number of times to try to find a better
    /// expression even though the cost hasn't improved.
    /// If the expression is the exact same, then
    /// we return it.
    /// Say we found an expression with cost 12.
    /// If we find an expression with cost 12 again,
    /// we will try again for this number of times,
    /// as long as the expression changes.
    /// Tje default is 3.
    retry_same_cost: usize,
    /// The maximum number of nodes in the egraph.
    /// The default is 10'000.
    node_limit: usize,
    /// The maximum number of iterations for the runner.
    /// The default is 30.
    iter_limit: usize,
    /// The maximum time to run for.
    /// The default is 1s.
    time_limit: std::time::Duration,
}

impl SimplificationConfig {
    /// Returns the default configuration.
    pub fn default() -> Self {
        Self {
            retry_limit: 10,
            retry_same_cost: 3,
            node_limit: 10_000,
            iter_limit: 30,
            time_limit: std::time::Duration::from_millis(1000),
        }
    }
}

/// Language for boolean expressions.
define_language! {
    enum BooleanLanguage {
        "0" = Zero,
        "1" = One,
        "not" = Not(Id),
        "and" = And([Id; 2]),
        "or" = Or([Id; 2]),
        "xor" = Xor([Id; 2]),
        Symbol(Symbol),
    }
}

type BExpr = RecExpr<BooleanLanguage>;

// Helper function that checks for equality of two BExprs.
// This is not optimized for common subexpressions.
fn bexpr_eq(l: &BExpr, r: &BExpr) -> bool {
    fn bexpr_eq_impl(l: &BExpr, lid: Id, r: &BExpr, rid: Id) -> bool {
        match (&l[lid], &r[rid]) {
            (BooleanLanguage::Zero, BooleanLanguage::Zero) => true,
            (BooleanLanguage::One, BooleanLanguage::One) => true,
            (BooleanLanguage::Not(lid), BooleanLanguage::Not(rid)) =>
                bexpr_eq_impl(l, *lid, r, *rid),
            (BooleanLanguage::And([ll, lr]), BooleanLanguage::And([rl, rr])) =>
                bexpr_eq_impl(l, *ll, r, *rl) && bexpr_eq_impl(l, *lr, r, *rr),
            (BooleanLanguage::Or([ll, lr]), BooleanLanguage::Or([rl, rr])) =>
                bexpr_eq_impl(l, *ll, r, *rl) && bexpr_eq_impl(l, *lr, r, *rr),
            (BooleanLanguage::Xor([ll, lr]), BooleanLanguage::Xor([rl, rr])) =>
                bexpr_eq_impl(l, *ll, r, *rl) && bexpr_eq_impl(l, *lr, r, *rr),
            (BooleanLanguage::Symbol(l), BooleanLanguage::Symbol(r)) => l == r,
            _ => false,
        }
    }

    bexpr_eq_impl(
        l, Id::from(l.as_ref().len() - 1),
        r, Id::from(r.as_ref().len() - 1)
    )
}

/// The set of equivalences.
lazy_static::lazy_static! {
    static ref RULES: Vec<Rewrite<BooleanLanguage, ()>> = make_rules();
}

/// Returns a simplified boolean expression for a given truth table.
/// The truth table is an iterator over the value of the rows.
/// This function iterates once over the table to determine
/// whether the DNF or the CNF is smaller and then calls
/// `simplify_from_truth_table_dnf` or `simplify_from_truth_table_cnf`.
/// Thus the iterator needs to be cloneable.
pub fn simplify_from_truth_table<Table>(
    vars: &[Symbol], table: Table, cfg: &SimplificationConfig
) -> UExpr
where
    Table: Iterator<Item = bool> + Clone,
{
    let (size, ones) = table.clone().fold((0, 0),
        |(size, ones), v| (size + 1usize, ones + v as usize));

    if ones >= size / 2 {
        simplify_from_truth_table_cnf(
            vars,
            table.enumerate().filter(|(_, v)| !*v).map(|(i, _)| i),
            cfg
        )
    } else {
        simplify_from_truth_table_dnf(
            vars,
            table.enumerate().filter(|(_, v)| *v).map(|(i, _)| i),
            cfg
        )
    }
}

/// Returns a simplified boolean expression for a given truth table.
/// The truth table is an iterator over the indices of the rows
/// that have a value of 1.
/// This constructs the DNF and then simplifies it using algebra.
pub fn simplify_from_truth_table_dnf<Table>(
    vars: &[Symbol], table: Table, cfg: &SimplificationConfig
) -> UExpr
where
    Table: Iterator<Item = usize>,
{
    let mut nodes = Vec::new();
    for v in vars {
        nodes.push(BooleanLanguage::Symbol(*v));
        nodes.push(BooleanLanguage::Not(Id::from(nodes.len() - 1)));
    }

    let mut nodes = RefCell::new(nodes);

    fn from_1(
        vars: &[Symbol], i: usize, nodes: &RefCell<Vec<BooleanLanguage>>
    ) -> usize {
        let mut iter = (0..vars.len()).map(|j|
            j * 2 + (1 - ((i >> j) & 1))
        );

        let first = iter.next().unwrap();
        iter.fold(first, |acc, x| {
            let mut nodes = nodes.borrow_mut();
            nodes.push(BooleanLanguage::And([Id::from(acc), Id::from(x)]));
            nodes.len() - 1
        })
    }

    let mut iter = table.map(|i| from_1(vars, i, &nodes));
    let Some(first) = iter.next() else {
        return UExpr::not(UExpr::Ones);
    };
    let root = iter.fold(first, |acc, x| {
        let mut nodes = nodes.borrow_mut();
        nodes.push(BooleanLanguage::Or([Id::from(acc), Id::from(x)]));
        nodes.len() - 1
    });

    // If there is only one variable `a` and the expression is `a`,
    // then the variable entry will the root. But the vector already
    // contains `not a`, so we need to remove it, so that `a` is the root.
    nodes.borrow_mut().truncate(root + 1);

    let dnf = RecExpr::from(nodes.into_inner());
    bexpr_to_uexpr(&simplify_bexpr(&dnf, cfg))
}

/// Returns a simplified boolean expression for a given truth table.
/// The truth table is an iterator over the indices of the rows
/// that have a value of 0.
/// This constructs the DNF and then simplifies it using algebra.
pub fn simplify_from_truth_table_cnf<Table>(
    vars: &[Symbol], table: Table, cfg: &SimplificationConfig
) -> UExpr
where
    Table: Iterator<Item = usize>,
{
    let mut nodes = Vec::new();
    for v in vars {
        nodes.push(BooleanLanguage::Symbol(*v));
        nodes.push(BooleanLanguage::Not(Id::from(nodes.len() - 1)));
    }

    let mut nodes = RefCell::new(nodes);

    fn from_0(
        vars: &[Symbol], i: usize, nodes: &RefCell<Vec<BooleanLanguage>>
    ) -> usize {
        let mut iter = (0..vars.len()).map(|j|
            j * 2 + ((i >> j) & 1)
        );

        let first = iter.next().unwrap();
        iter.fold(first, |acc, x| {
            let mut nodes = nodes.borrow_mut();
            nodes.push(BooleanLanguage::Or([Id::from(acc), Id::from(x)]));
            nodes.len() - 1
        })
    }

    let mut iter = table.map(|i| from_0(vars, i, &nodes));
    let Some(first) = iter.next() else {
        return UExpr::Ones;
    };
    let root = iter.fold(first, |acc, x| {
        let mut nodes = nodes.borrow_mut();
        nodes.push(BooleanLanguage::And([Id::from(acc), Id::from(x)]));
        nodes.len() - 1
    });

    // If there is only one variable `a` and the expression is `a`,
    // then the variable entry will the root. But the vector already
    // contains `not a`, so we need to remove it, so that `a` is the root.
    nodes.borrow_mut().truncate(root + 1);

    let dnf = RecExpr::from(nodes.into_inner());
    bexpr_to_uexpr(&simplify_bexpr(&dnf, cfg))
}

/// Simplify a boolean expression.
pub fn simplify(e: &UExpr, cfg: &SimplificationConfig) -> UExpr {
    let e = uexpr_to_bexpr(e);
    let best = simplify_bexpr(&e, cfg);
    bexpr_to_uexpr(&best)
}

/// Simplifies a BExpr.
fn simplify_bexpr(e: &BExpr, cfg: &SimplificationConfig) -> BExpr {
    // The current expression with the best cost and its cost.
    log::trace!("Simplifying {e}");
    let mut best_cost = usize::MAX;
    let mut best = e.clone();
    let mut best_cost_iter = 0usize;
    for i in 0.. {
        let mut runner = Runner::default()
            .with_iter_limit(cfg.iter_limit)
            .with_node_limit(cfg.node_limit)
            .with_time_limit(cfg.time_limit)
            .with_expr(&best)
            .run(&*RULES);

        let stop_reason = runner.stop_reason.unwrap();
        log::trace!("{:?}", stop_reason);

        let extractor = Extractor::new(&runner.egraph, AstSize);

        let root = runner.roots[0];
        let (cost, new_best) = extractor.find_best(root);

        log::trace!("{new_best} with cost {cost}");
        assert!(cost <= best_cost, "Cost should not increase");

        if cost == best_cost {
            if bexpr_eq(&best, &new_best) {
                return new_best;
            } else {
                best_cost_iter += 1;
            }
        } else {
            best_cost_iter = 0;
        }

        if best_cost_iter >= cfg.retry_same_cost
            || matches!(stop_reason, StopReason::Saturated)
            || i == cfg.retry_limit - 1 {
            return new_best;
        }
        best = new_best;
        best_cost = cost;
    }

    unreachable!();
}

/// Converts a uniform expression to a egg boolean expression.
fn uexpr_to_bexpr(e: &UExpr) -> BExpr {
    fn uexpr_to_bexpr_impl(e: &UExpr, nodes: &mut Vec<BooleanLanguage>) -> Id {
        match e {
            UExpr::Var(v) => {
                nodes.push(BooleanLanguage::Symbol(*v));
            },
            UExpr::Ones => nodes.push(BooleanLanguage::One),
            UExpr::Not(e) => {
                let e = uexpr_to_bexpr_impl(e, nodes);
                nodes.push(BooleanLanguage::Not(e));
            },
            UExpr::And(l, r) => {
                let l = uexpr_to_bexpr_impl(l, nodes);
                let r = uexpr_to_bexpr_impl(r, nodes);
                nodes.push(BooleanLanguage::And([l, r]));
            },
            UExpr::Or(l, r) => {
                let l = uexpr_to_bexpr_impl(l, nodes);
                let r = uexpr_to_bexpr_impl(r, nodes);
                nodes.push(BooleanLanguage::Or([l, r]));
            },
            UExpr::Xor(l, r) => {
                let l = uexpr_to_bexpr_impl(l, nodes);
                let r = uexpr_to_bexpr_impl(r, nodes);
                nodes.push(BooleanLanguage::Xor([l, r]));
            },
        }

        Id::from(nodes.len() - 1)
    }

    let mut nodes = Vec::new();
    uexpr_to_bexpr_impl(e, &mut nodes);
    RecExpr::from(nodes)
}

/// Converts an egg boolean expression into a uniform expression.
fn bexpr_to_uexpr(e: &BExpr) -> UExpr {
    fn bexpr_to_uexpr_impl(e: &BExpr, idx: Id) -> UExpr {
        match &e[idx] {
            BooleanLanguage::Zero => UExpr::not(UExpr::Ones),
            BooleanLanguage::One => UExpr::Ones,
            BooleanLanguage::Not(i) => UExpr::not(bexpr_to_uexpr_impl(e, *i)),
            BooleanLanguage::And([i, j]) =>
                UExpr::and(bexpr_to_uexpr_impl(e, *i), bexpr_to_uexpr_impl(e, *j)),
            BooleanLanguage::Or([i, j]) =>
                UExpr::or(bexpr_to_uexpr_impl(e, *i), bexpr_to_uexpr_impl(e, *j)),
            BooleanLanguage::Xor([i, j]) =>
                UExpr::xor(bexpr_to_uexpr_impl(e, *i), bexpr_to_uexpr_impl(e, *j)),
            BooleanLanguage::Symbol(s) => UExpr::var(s.to_string()),
        }
    }

    bexpr_to_uexpr_impl(e, Id::from(e.as_ref().len() - 1))
}

/// The equivalences.
/// Some are definitely redundant,
/// but I think it is better to have more than less.
fn make_rules() -> Vec<Rewrite<BooleanLanguage, ()>> {
    let mut rules = vec![
        rewrite!("and-0"; "(and ?x 0)" => "0"),
        rewrite!("or-1"; "(or ?x 1)" => "1"),
        rewrite!("xor-self"; "(xor ?x ?x)" => "0"),

        // Complementation
        rewrite!("and-not"; "(and ?x (not ?x))" => "0"),
        rewrite!("or-not"; "(or ?x (not ?x))" => "1"),
        rewrite!("xor-not"; "(xor ?x (not ?x))" => "1"),

        // Absorption.
        rewrite!("and-absorb"; "(and ?x (or ?x ?y))" => "?x"),
        rewrite!("or-absorb"; "(or ?x (and ?x ?y))" => "?x"),
    ];

    rules.extend(vec![
        // Associativity
        rewrite!("and-assoc"; "(and ?x (and ?y ?z))" <=> "(and (and ?x ?y) ?z)"),
        rewrite!("or-assoc"; "(or ?x (or ?y ?z))" <=> "(or (or ?x ?y) ?z)"),
        rewrite!("xor-assoc"; "(xor ?x (xor ?y ?z))" <=> "(xor (xor ?x ?y) ?z)"),

        // Commutativity
        rewrite!("and-comm"; "(and ?x ?y)" <=> "(and ?y ?x)"),
        rewrite!("or-comm"; "(or ?x ?y)" <=> "(or ?y ?x)"),
        rewrite!("xor-comm"; "(xor ?x ?y)" <=> "(xor ?y ?x)"),

        rewrite!("and-1"; "(and ?x 1)" <=> "?x"),
        rewrite!("and-idemp"; "(and ?x ?x)" <=> "?x"),

        rewrite!("or-0"; "(or ?x 0)" <=> "?x"),
        rewrite!("or-idemp"; "(or ?x ?x)" <=> "?x"),

        rewrite!("xor-0"; "(xor ?x 0)" <=> "?x"),
        rewrite!("xor-1"; "(xor ?x 1)" <=> "(not ?x)"),

        // Distributivity
        rewrite!("and-or"; "(and ?x (or ?y ?z))" <=> "(or (and ?x ?y) (and ?x ?z))"),
        rewrite!("or-and"; "(or ?x (and ?y ?z))" <=> "(and (or ?x ?y) (or ?x ?z))"),
        rewrite!("and-xor"; "(and ?x (xor ?y ?z))" <=> "(xor (and ?x ?y) (and ?x ?z))"),

        // Double negation
        rewrite!("not-not"; "(not (not ?x))" <=> "?x"),

        // De Morgan's laws
        rewrite!("not-and"; "(not (and ?x ?y))" <=> "(or (not ?x) (not ?y))"),
        rewrite!("not-or"; "(not (or ?x ?y))" <=> "(and (not ?x) (not ?y))"),
        rewrite!("not-xor"; "(not (xor ?x ?y))" <=> "(xor (not ?x) ?y)"),

        rewrite!("and-xor-1"; "(and ?x (xor ?x ?y))" <=> "(and ?x (not ?y))"),
        rewrite!("or-xor-1"; "(or ?x (xor ?x ?y))" <=> "(or ?x ?y)"),
        rewrite!("xor-and-1"; "(xor ?x (and ?x ?y))" <=> "(and ?x (not ?y))"),
        rewrite!("xor-or-1"; "(xor ?x (or ?x ?y))" <=> "(and (not ?x) ?y)"),
    ].concat());

    rules
}

fn simplify_explanation_test() {
    let expr: BExpr = "(and ?x (xor ?x ?y))".parse().unwrap();
    println!("{expr}");

    let mut runner = Runner::default()
        .with_explanations_enabled()
        .with_expr(&expr)
        .run(&*RULES);

    // Print the total egraph size.
    println!("Total egraph size: {}", runner.egraph.total_size());

    let root = runner.roots[0];

    let extractor = Extractor::new(&runner.egraph, AstSize);
    let (best_cost, best) = extractor.find_best(root);
    println!("Simplified {} to {} with cost {}", expr, best, best_cost);

    println!("{}", runner.explain_equivalence(&expr, &best));
}

/// Verifies that all the equivalences are indeed equivalences.
#[test]
fn verify_boolean_equivalences() {
    fn ast_to_uexpr(e: &PatternAst<BooleanLanguage>) -> UExpr {
        fn ast_to_uexpr_impl(e: &PatternAst<BooleanLanguage>, idx: Id) -> UExpr {
            let r = match &e[idx] {
                ENodeOrVar::Var(v) => return UExpr::Var(v.to_string().into()),
                ENodeOrVar::ENode(e) => e,
            };

            match r {
                BooleanLanguage::Zero => UExpr::not(UExpr::Ones),
                BooleanLanguage::One => UExpr::Ones,
                BooleanLanguage::Not(i) => UExpr::not(ast_to_uexpr_impl(e, *i)),
                BooleanLanguage::And([i, j]) =>
                    UExpr::and(ast_to_uexpr_impl(e, *i), ast_to_uexpr_impl(e, *j)),
                BooleanLanguage::Or([i, j]) =>
                    UExpr::or(ast_to_uexpr_impl(e, *i), ast_to_uexpr_impl(e, *j)),
                BooleanLanguage::Xor([i, j]) =>
                    UExpr::xor(ast_to_uexpr_impl(e, *i), ast_to_uexpr_impl(e, *j)),
                BooleanLanguage::Symbol(s) => panic!("There should be no free variables."),
            }
        }

        ast_to_uexpr_impl(e, Id::from(e.as_ref().len() - 1))
    }

    fn verify(lhs: &PatternAst<BooleanLanguage>, rhs: &PatternAst<BooleanLanguage>) {
        let lhs = ast_to_uexpr(lhs);
        let rhs = ast_to_uexpr(rhs);
        println!("Verifying equivalence: {lhs} == {rhs}");
        let mut vars = std::collections::BTreeSet::new();
        lhs.vars_impl(&mut vars);
        rhs.vars_impl(&mut vars);
        let mut vars: Vec<_> = vars.into_iter().collect();
        let mut val = Valuation::zero();

        assert!(vars.len() <= 64);

        for i in 0..(1 << vars.len()) {
            for (j, v) in vars.iter().enumerate() {
                val.set_value(*v, -rug::Integer::from((i >> j) & 1));
            }

            assert_eq!(lhs.eval(&mut val), rhs.eval(&mut val),
                "Invalid boolean equivalence: {lhs} != {rhs} for {val:?}");
        }
    }

    for rule in &*RULES {
        let Some(lhs) = rule.searcher.get_pattern_ast() else {
            continue;
        };

        let Some(rhs) = rule.applier.get_pattern_ast() else {
            continue;
        };

        verify(lhs, rhs);

    }
}

fn n_vars(n: usize) -> Vec<Symbol> {
    if n <= 26 {
        return (b'a'..=b'z').take(n)
            .map(|c| Symbol::from(std::str::from_utf8(&[c]).unwrap()))
            .collect();
    }

    (1..=n).map(|i| Symbol::from(format!("x{}", i))).collect()
}

// Used during development and I didn't want to delete it.
fn short_from_rand_truth_table() {
    let n = rand::random::<usize>() % 5 + 1;
    let vars = n_vars(n);

    let truth_table: Vec<_> = std::iter::repeat_with(rand::random::<bool>)
        .take(1 << n)
        .collect();
    let cfg = SimplificationConfig {
        node_limit: 1_000_000,
        iter_limit: 100,
        time_limit: std::time::Duration::from_secs(60),
        ..SimplificationConfig::default()
    };
    let e = simplify_from_truth_table(
        &vars,
        truth_table.iter().cloned(),
        &cfg
    );

    println!("{e}");
}

// Used during development and I didn't want to delete it.
fn list_minimal_expressions() {
    let vars = n_vars(3);
    let mut truth_table = vec![false; 1 << vars.len()];
    let cfg = SimplificationConfig::default();

    loop {
        let e = simplify_from_truth_table(&vars, truth_table.iter().cloned(), &cfg);
        for v in &truth_table {
            print!("{}", *v as usize);
        }
        println!(": {e}");

        // Is it all ones?
        if truth_table.iter().all(|v| *v) {
            break;
        }

        // Binary addition.
        for v in &mut truth_table {
            *v = !*v;
            if *v {
                break;
            }
        }
    }
}