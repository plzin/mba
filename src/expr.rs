#![allow(dead_code)]

use std::{rc::Rc, collections::BTreeSet};

use rug::Integer;

/// Used internally to convert an integer from and iterator over characters.
pub(crate) fn int_from_it(
    it: &mut std::iter::Peekable<std::str::Chars>
) -> Option<Integer> {
    let c = it.next()?;
    let mut num: Integer = c.to_digit(10)?.into();

    let mut lookahead = it.peek();
    while let Some(c) = lookahead {
        if !c.is_ascii_digit() {
            break;
        }

        num *= 10;
        num += c.to_digit(10).unwrap();
        it.next();
        lookahead = it.peek();
    }

    Some(num)
}

/// A general expression.
/// We use Rc instead of Box in order to avoid copying common subexpressions.
#[derive(Clone, Debug)]
pub enum Expr {
    Const(Integer),
    Var(String),
    Add(Rc<Expr>, Rc<Expr>),
    Sub(Rc<Expr>, Rc<Expr>),
    Mul(Rc<Expr>, Rc<Expr>),
    Div(Rc<Expr>, Rc<Expr>),
    Neg(Rc<Expr>),
    And(Rc<Expr>, Rc<Expr>),
    Or(Rc<Expr>, Rc<Expr>),
    Xor(Rc<Expr>, Rc<Expr>),
    Not(Rc<Expr>),
}

impl Expr {
    /// Returns the zero constant.
    pub const fn zero() -> Self {
        Self::Const(Integer::new())
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
            Expr::Const(_) => {},
            Expr::Var(name) => drop(v.insert(name.clone())),
            Expr::Neg(e) | Expr::Not(e) => e.vars_impl(v),
            Expr::Add(l, r) | Expr::Sub(l, r) | Expr::Mul(l, r)
            | Expr::Div(l, r) | Expr::And(l, r) | Expr::Or(l, r)
            | Expr::Xor(l, r) => {
                l.vars_impl(v);
                r.vars_impl(v);
            }
        }
    }

    /// Substitutes an expression for a variable.
    /// If the Rc's in this expression are shared with
    /// other expressions then this will also substitute in those.
    pub fn substitute(mut self, var: &str, e: &mut Rc<Expr>) -> Self {
        self.substitute_mut(var, e);
        self
    }

    /// Substitutes an expression for a variable.
    /// If the Rc's in this expression are shared with
    /// other expressions then this will also substitute in those.
    pub fn substitute_mut(&mut self, var: &str, e: &mut Rc<Expr>) {
        let mut visited = Vec::new();

        use Expr::*;
        match self {
            Const(_) => (),
            Var(v) => if v == var { *self = (**e).clone() },
            Add(l, r) | Sub(l, r) | Mul(l, r) | Div(l, r)
                | And(l, r) | Or(l, r) | Xor(l, r) => {
                    Self::substitute_impl(l, var, e, &mut visited);
                    Self::substitute_impl(r, var, e, &mut visited);
            },
            Neg(i) | Not(i)
                => Self::substitute_impl(i, var, e, &mut visited),
        }
    }

    fn substitute_impl(
        this: &mut Rc<Expr>, var: &str,
        e: &mut Rc<Expr>, visited: &mut Vec<*const Expr>
    ) {
        let ptr = Rc::as_ptr(this);
        let recurse = if visited.contains(&ptr) {
            false
        } else {
            visited.push(ptr);
            true
        };


        use Expr::*;
        // SAFETY: This is okay because we make sure with extra logic
        // that this is never encountered twice.
        match unsafe { &mut *(ptr as *mut _) } {
            Const(_) => (),
            Var(v) => if *v == var { *this = e.clone() },
            Add(l, r) | Sub(l, r) | Mul(l, r) | Div(l, r)
                | And(l, r) | Or(l, r) | Xor(l, r) => if recurse {
                    Self::substitute_impl(l, var, e, visited);
                    Self::substitute_impl(r, var, e, visited);
            },
            Neg(i) | Not(i) => if recurse {
                Self::substitute_impl(i, var, e, visited)
            },
        }
    }

    /// Prints the expression while avoiding to reprint
    /// common subexpressions by assigning them to variables.
    /// This only works if the Rc's used in the expression
    /// are not shared with other Exprs.
    pub fn print_simple(&self) {
        // Stores a mapping of (sub)expressions to variables.
        let mut vars = Vec::new();

        let l = self.print_simple_impl(&mut vars);

        for (_, var, init) in vars.iter().rev() {
            println!("{} = {}", var, init);
        }

        println!("{}", l);
    }

    fn print_simple_rc(
        e: &Rc<Expr>,
        vars: &mut Vec<(*const Expr, char, String)>
    ) -> String {
        // If there is only one reference then just print it.
        if Rc::strong_count(e) == 1 {
            return e.print_simple_impl(vars);
        }

        // We don't want to assign a variable to a variable
        // so there is this shortcut here.
        if let Expr::Var(v) = &**e {
            return format!("{}", *v);
        }

        let ptr = Rc::as_ptr(e);

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
        &self, vars: &mut Vec<(*const Expr, char, String)>
    ) -> String {
        let bin_op = |
            op: char, l: &Rc<Expr>, r: &Rc<Expr>,
            vars: &mut Vec<(*const Expr, char, String)>
        | {
            let pred = self.precedence();

            let l = if pred > l.precedence() && Rc::strong_count(l) == 1 {
                format!("({})", Expr::print_simple_rc(l, vars))
            } else {
                format!("{}", Expr::print_simple_rc(l, vars))
            };

            let r = if pred > r.precedence() && Rc::strong_count(r) == 1 {
                format!("({})", Expr::print_simple_rc(r, vars))
            } else {
                format!("{}", Expr::print_simple_rc(r, vars))
            };

            format!("{} {} {}", l, op, r)
        };

        let un_op = |
            op: char, i: &Rc<Expr>,
            vars: &mut Vec<(*const Expr, char, String)>
        | {
            if self.precedence() > i.precedence() && Rc::strong_count(i) == 1 {
                format!("{}({})", op, Expr::print_simple_rc(i, vars))
            } else {
                format!("{}{}", op, Expr::print_simple_rc(i, vars))
            }
        };

        use Expr::*;
        match self {
            Const(i) => format!("{}", i),
            Var(n) => format!("{}", n),
            Add(l, r) => bin_op('+', l, r, vars),
            Sub(l, r) => bin_op('-', l, r, vars),
            Mul(l, r) => bin_op('*', l, r, vars),
            Div(l, r) => bin_op('/', l, r, vars),
            Neg(i) => un_op('-', i, vars),
            And(l, r) => bin_op('&', l, r, vars),
            Or(l, r) => bin_op('|', l, r, vars),
            Xor(l, r) => bin_op('^', l, r, vars),
            Not(i) => un_op('~', i, vars),
        }
    }

    pub fn from_string<T: ToString>(s: T) -> Option<Expr> {
        let mut s = s.to_string();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    /// Returns the precedence of a binary operator.
    /// All operators are taken to be left associative.
    fn precedence(&self) -> usize {
        use Expr::*;
        match self {
            Or(_, _) => 1,
            Xor(_, _) => 2,
            And(_, _) => 3,
            Add(_, _) | Sub(_, _) => 4,
            Mul(_, _) | Div(_, _) => 5,
            Neg(_) | Not(_) => 15,
            Const(_) | Var(_) => 16,
        }
    }

    // pre 0: parse as much as possible
    // ...
    // pre 15: parse as little as possible
    fn parse(
        it: &mut std::iter::Peekable<std::str::Chars>,
        pre: usize
    ) -> Option<Expr> {
        use Expr::*;

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
            Not(Rc::new(e))
        } else if c == '-' {
            it.next();
            let e = Self::parse(it, 15)?;
            Neg(Rc::new(e))
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
        } else if c.is_ascii_digit() {
            let num = int_from_it(it)?;
            Const(num)
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
                '+' | '-' => 4,
                '*' | '/' => 5,
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
            let rhs = Rc::new(Self::parse(it, op_pre)?);
            let lhs = Rc::new(e);
            e = match c {
                '+' => Add(lhs, rhs),
                '-' => Sub(lhs, rhs),
                '*' => Mul(lhs, rhs),
                '/' => Div(lhs, rhs),
                '&' => And(lhs, rhs),
                '|' => Or(lhs, rhs),
                '^' => Xor(lhs, rhs),
                _ => unreachable!(),
            };
        }
    }
}
