#![allow(dead_code)]

use std::rc::Rc;

use rug::Integer;

/// Used internally to convert an integer from and iterator over characters.
pub(crate) fn int_from_it(
    it: &mut std::iter::Peekable<std::str::Chars>
) -> Option<Integer> {
    let c = it.next()?;
    let mut num: Integer = c.to_digit(10).unwrap().into();

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
    Var(char),
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
    pub fn from_string<T: ToString>(s: T) -> Option<Expr> {
        let mut s = s.to_string();
        s.retain(|c| !c.is_whitespace());
        let mut it = s.chars().peekable();

        Self::parse(&mut it, 0)
    }

    /// Returns the precedence of a binary operator.
    /// All operators are taken to be left associative.
    fn precedence(c: char) -> usize {
        match c {
            '|' => 1,
            '^' => 2,
            '&' => 3,
            '+' | '-' => 4,
            '*' | '/' => 5,
            _ => 0
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
            Var(c)
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

            let op_pre = Self::precedence(c);
            if op_pre <= pre {
                return Some(e);
            }

            // If the current operators precedence is higher than
            // the one whos subexpression we are currently parsing
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
                _ => return None,
            };
        }
    }
}
