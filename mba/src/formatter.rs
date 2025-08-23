//! Formatting for expressions.

use std::fmt::{Display, Write as _};

use crate::{
    Symbol,
    bitwise_expr::{BExpr, LBExpr},
    expr::{Expr, ExprOp},
    rings::{Ring, RingElement},
};

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Formatter {
    C,
    Rust,
    Tex,
    LLVM,
}

impl Formatter {
    pub fn add_op(&self) -> &str {
        "+"
    }

    pub fn sub_op(&self) -> &str {
        "-"
    }

    pub fn mul_op(&self) -> &str {
        match self {
            Formatter::C | Formatter::Rust => "*",
            Formatter::Tex => "\\cdot",
            Formatter::LLVM => panic!("LLVM does not support multiplication"),
        }
    }

    pub fn neg_op(&self) -> &str {
        "-"
    }

    pub fn and_op(&self) -> &str {
        match self {
            Formatter::C | Formatter::Rust => "&",
            Formatter::Tex => "\\land",
            Formatter::LLVM => panic!("LLVM does not support bitwise and"),
        }
    }

    pub fn or_op(&self) -> &str {
        match self {
            Formatter::C | Formatter::Rust => "|",
            Formatter::Tex => "\\lor",
            Formatter::LLVM => panic!("LLVM does not support bitwise or"),
        }
    }

    pub fn xor_op(&self) -> &str {
        match self {
            Formatter::C | Formatter::Rust => "^",
            Formatter::Tex => "\\oplus",
            Formatter::LLVM => panic!("LLVM does not support bitwise xor"),
        }
    }

    pub fn not_op(&self) -> &str {
        match self {
            Formatter::C => "~",
            Formatter::Rust => "!",
            Formatter::Tex => "\\neg",
            Formatter::LLVM => panic!("LLVM does not support bitwise not"),
        }
    }
}

impl BExpr {
    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// `formatter` to format the expression.
    pub fn display(&self, formatter: Formatter) -> DisplayableBExpr<'_> {
        DisplayableBExpr { expr: self, formatter }
    }

    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// formatter to format the expression. It will wrap the output in
    /// parentheses if the expression is not a topmost unary. This is basically
    /// a helper function.
    pub fn display_wrapped(
        &self,
        formatter: Formatter,
    ) -> DisplayableBExprWrapped<'_> {
        DisplayableBExprWrapped { expr: self, formatter }
    }
}

pub struct DisplayableBExpr<'a> {
    expr: &'a BExpr,
    formatter: Formatter,
}

impl<'a> DisplayableBExpr<'a> {
    fn handle_binop(
        &self,
        l: &BExpr,
        r: &BExpr,
        op: &str,
        f: &mut std::fmt::Formatter<'_>,
    ) -> std::fmt::Result {
        write!(
            f,
            "{} {op} {}",
            l.display_wrapped(self.formatter),
            r.display_wrapped(self.formatter),
        )
    }
}

impl<'a> Display for DisplayableBExpr<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.expr {
            BExpr::Ones => f.write_str("-1"),
            BExpr::Var(n) => n.fmt(f),
            BExpr::Not(e) => match self.formatter {
                Formatter::C => {
                    write!(f, "~{}", e.display_wrapped(self.formatter))
                },
                Formatter::Rust => {
                    write!(f, "!{}", e.display_wrapped(self.formatter))
                },
                Formatter::Tex => {
                    write!(f, "\\neg{{{}}}", e.display(self.formatter))
                },
                Formatter::LLVM => {
                    write!(f, "~{}", e.display_wrapped(Formatter::C))
                },
            },
            BExpr::And(l, r) => {
                self.handle_binop(l, r, self.formatter.and_op(), f)
            },
            BExpr::Or(l, r) => {
                self.handle_binop(l, r, self.formatter.or_op(), f)
            },
            BExpr::Xor(l, r) => {
                self.handle_binop(l, r, self.formatter.xor_op(), f)
            },
        }
    }
}

pub struct DisplayableBExprWrapped<'a> {
    expr: &'a BExpr,
    formatter: Formatter,
}

impl<'a> Display for DisplayableBExprWrapped<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.expr.is_topmost_unary() {
            self.expr.display(self.formatter).fmt(f)
        } else {
            write!(f, "({})", self.expr.display(self.formatter))
        }
    }
}

impl<R: Ring> LBExpr<R> {
    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// `formatter` to format the expression.
    pub fn display<'a>(
        &'a self,
        formatter: Formatter,
    ) -> LBExprFormatter<'a, R> {
        LBExprFormatter { expr: self, formatter }
    }

    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// `formatter` to format the expression as a function with the given
    /// `function_name`. This doesn't make sense for [`Formatter::Tex`] so it is
    /// the same as as [`Self::display`] in that case, except that it is slower
    /// because it computes the list of variables on construction.
    pub fn display_function<'a>(
        &'a self,
        formatter: Formatter,
        function_name: &'a str,
        r: &'a R,
    ) -> FunctionFormatter<'a, R, LBExprFormatter<'a, R>> {
        FunctionFormatter::new(
            self.vars(),
            self.display(formatter),
            formatter,
            function_name,
            r,
        )
    }
}

pub struct LBExprFormatter<'a, R: Ring> {
    expr: &'a LBExpr<R>,
    formatter: Formatter,
}

impl<'a, R: Ring> LBExprFormatter<'a, R> {
    fn write_term(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        c: &R::Element,
        e: &BExpr,
    ) -> std::fmt::Result {
        if !c.is_one() {
            match self.formatter {
                Formatter::Rust => write!(f, "Wrapping({c}) * ")?,
                Formatter::C => write!(f, "{c} * ")?,
                Formatter::Tex => write!(f, "{c}\\cdot")?,
                Formatter::LLVM => write!(f, "{c} * ")?,
            }
        }

        write!(f, "{}", e.display_wrapped(self.formatter))
    }
}

impl<'a, R: Ring> Display for LBExprFormatter<'a, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.expr.0.iter().filter(|(c, _)| !c.is_zero());

        let Some((c, e)) = iter.next() else {
            return f.write_str("0");
        };

        // There is one important special case where there is only a single term
        // whose coefficient is 1. In that case we never want to use outer
        // parentheses. This could be handled more generally by remembering the
        // current operator precedence, but that would mean more code and
        // usually we want parentheses for clarity even if they are not needed.
        // So I will just special-case it with some ugly code here.
        if c.is_one() {
            // If there is another term we just use normal formatting for both.
            if let Some((next_c, next_e)) = iter.next() {
                // We don't need to use `write_term` here, we already know the
                // coefficient is 1.
                write!(f, "{}", e.display_wrapped(self.formatter))?;
                f.write_str(" + ")?;
                self.write_term(f, next_c, next_e)?;
            }
            // If not, we don't use `display_wrapped`.
            else {
                write!(f, "{}", e.display(self.formatter))?;
            }
        } else {
            self.write_term(f, c, e)?;
        }

        for (c, e) in iter {
            f.write_str(" + ")?;
            self.write_term(f, c, e)?;
        }

        Ok(())
    }
}

impl<R: Ring> Expr<R> {
    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// `formatter` to format the expression. Each line will be indented by
    /// `tabs` spaces. Common subexpressions will be assigned to variables
    /// with names `v1`, `v2`, etc. on their own lines. The final result will
    /// be prefixed with `prefix`, which will usually be something like
    /// `"return "` or `"uint8_t result = "`. There will be no trailing newline.
    pub fn display<'a>(
        &'a self,
        formatter: Formatter,
        tabs: usize,
        prefix: &'a str,
        r: &'a R,
    ) -> ExprFormatter<'a, R> {
        let mut formatter = ExprFormatter {
            expr: self,
            formatter,
            buf: String::new(),
            subs: Vec::new(),
            tabs,
            prefix,
            r,
            llvm_add: 0,
            llvm_sub: 0,
            llvm_mul: 0,
            llvm_and: 0,
            llvm_or: 0,
            llvm_xor: 0,
        };

        formatter.format();

        formatter
    }

    /// Creates a wrapper struct that implements [`Display`] and uses the given
    /// `formatter` to format the expression as a function with the given
    /// `function_name`. This doesn't make sense for [`Formatter::Tex`] so it is
    /// the same as as [`Self::display`] in that case, except that it is slower
    /// because it computes the list of variables on construction.
    pub fn display_function<'a>(
        &'a self,
        formatter: Formatter,
        function_name: &'a str,
        r: &'a R,
    ) -> FunctionFormatter<'a, R, ExprFormatter<'a, R>> {
        let prefix = match formatter {
            Formatter::C => "return ",
            _ => "",
        };
        FunctionFormatter::new(
            self.vars(),
            self.display(formatter, 1, prefix, r),
            formatter,
            function_name,
            r,
        )
    }
}

struct CommonSubExpr<R: Ring> {
    /// The pointer to the subexpression that has a strong count > 1.
    ptr: *const ExprOp<R>,

    /// The variable name that is used to represent the subexpression.
    var: Symbol,

    /// The initializer for the variable.
    init: String,
}

pub struct ExprFormatter<'a, R: Ring> {
    expr: &'a Expr<R>,
    formatter: Formatter,
    buf: String,
    subs: Vec<CommonSubExpr<R>>,
    tabs: usize,
    prefix: &'a str,
    r: &'a R,
    // LLVM SSA naming counters by op mnemonic
    llvm_add: usize,
    llvm_sub: usize,
    llvm_mul: usize,
    llvm_and: usize,
    llvm_or: usize,
    llvm_xor: usize,
}

impl<'a, R: Ring> ExprFormatter<'a, R> {
    fn format(&mut self) {
        if self.formatter == Formatter::LLVM {
            let result = self.llvm_emit_value(self.expr);
            self.buf.clear();
            self.buf.push_str(&result);
            return;
        }
        // If there is only one reference then just print it.
        if self.expr.strong_count() == 1 {
            return self.format_op();
        }

        // We don't want to assign a variable to a variable
        // so there is this shortcut here.
        if let ExprOp::Var(v) = self.expr.as_ref() {
            write!(&mut self.buf, "{v}").unwrap();
            return;
        }

        let ptr = self.expr.as_ptr();

        let sub = self.subs.iter().find(|t| t.ptr == ptr);
        if let Some(sub) = sub {
            // If the expression already has a variable then just write the
            // variable name.
            write!(&mut self.buf, "{}", sub.var).unwrap();
        } else {
            // Otherwise, we need to create a new variable.
            // We can avoid the allocation by using our `self.buf` which we need
            // to write the variable name into anyways, but subslicing it does
            // involve some UTF-8 logic, which should be cheaper than allocating
            // though.
            let name_start = self.buf.len();
            write!(&mut self.buf, "v{}", self.subs.len() + 1).unwrap();
            let v = Symbol::new(&self.buf[name_start..]);

            let idx = self.subs.len();

            // Push it first so that we don't try to create a variable for the
            // same subexpression again.
            self.subs.push(CommonSubExpr { ptr, var: v, init: String::new() });

            // Format the subexpression into a new buffer.
            let mut cur_buf = std::mem::take(&mut self.buf);

            // Format the subexpression into the new buffer.
            self.format_op();

            // Swap the buffers back. `cur_buf` now contains the formatted
            // subexpression.
            std::mem::swap(&mut self.buf, &mut cur_buf);

            // Set the initializer for the variable.
            self.subs[idx].init = cur_buf;
        }
    }

    fn format_op(&mut self) {
        macro_rules! format_bin_op {
            ($op:ident, $l:ident, $r:ident) => {{
                if self.formatter == Formatter::LLVM {
                    unreachable!()
                }
                let pred = self.expr.precedence();

                let old = std::mem::replace(&mut self.expr, $l);
                if pred > $l.precedence() && $l.strong_count() == 1 {
                    self.buf.push_str("(");
                    // We know `l` has a strong count == 1, so we can use
                    // `format_op`.
                    self.format_op();
                    self.buf.push_str(")");
                } else {
                    self.format();
                }

                write!(&mut self.buf, " {} ", self.formatter.$op()).unwrap();

                self.expr = $r;
                if pred > $r.precedence() && $r.strong_count() == 1 {
                    self.buf.push_str("(");
                    self.format_op();
                    self.buf.push_str(")");
                } else {
                    self.format();
                }

                self.expr = old;
            }};
        }

        macro_rules! format_un_op {
            ($op:ident, $e:ident) => {{
                let pred = self.expr.precedence();

                let old = std::mem::replace(&mut self.expr, $e);

                self.buf.push_str(self.formatter.$op());
                if pred > $e.precedence() && $e.strong_count() == 1 {
                    self.buf.push_str("(");
                    self.format_op();
                    self.buf.push_str(")");
                } else {
                    self.format();
                }

                self.expr = old;
            }};
        }

        match self.expr.as_ref() {
            ExprOp::Const(i) => {
                if self.formatter == Formatter::Rust {
                    write!(&mut self.buf, "Wrapping({i})").unwrap()
                } else {
                    write!(&mut self.buf, "{i}").unwrap()
                }
            },
            ExprOp::Var(v) => write!(&mut self.buf, "{v}").unwrap(),
            ExprOp::Add(l, r) => format_bin_op!(add_op, l, r),
            ExprOp::Sub(l, r) => format_bin_op!(sub_op, l, r),
            ExprOp::Mul(l, r) => format_bin_op!(mul_op, l, r),
            ExprOp::Neg(i) => format_un_op!(neg_op, i),
            ExprOp::And(l, r) => format_bin_op!(and_op, l, r),
            ExprOp::Or(l, r) => format_bin_op!(or_op, l, r),
            ExprOp::Xor(l, r) => format_bin_op!(xor_op, l, r),
            ExprOp::Not(i) => format_un_op!(not_op, i),
        }
    }

    fn llvm_type_name(&self) -> impl Display + use<'a, R> {
        self.r.data_type_name(Formatter::LLVM)
    }

    fn llvm_fresh_var(&mut self, op: &str) -> Symbol {
        let name = match op {
            "or" => {
                self.llvm_or += 1;
                format!("or{}", self.llvm_or)
            },
            "and" => {
                self.llvm_and += 1;
                format!("and{}", self.llvm_and)
            },
            "xor" => {
                self.llvm_xor += 1;
                format!("xor{}", self.llvm_xor)
            },
            "add" => {
                self.llvm_add += 1;
                format!("add{}", self.llvm_add)
            },
            "sub" => {
                self.llvm_sub += 1;
                format!("sub{}", self.llvm_sub)
            },
            "mul" => {
                self.llvm_mul += 1;
                format!("mul{}", self.llvm_mul)
            },
            _ => format!("v{}", self.subs.len() + 1),
        };
        Symbol::new(&name)
    }

    fn find_existing_var(&self, ptr: *const ExprOp<R>) -> Option<Symbol> {
        self.subs.iter().find(|s| s.ptr == ptr).map(|s| s.var)
    }

    fn llvm_emit_value(&mut self, e: &Expr<R>) -> String {
        use crate::expr::ExprOp::*;

        match e.as_ref() {
            Const(i) => format!("{i}"),
            Var(v) => format!("%{v}"),
            Add(l, r) => self.llvm_emit_bin("add", l, r, e.as_ptr()),
            Sub(l, r) => self.llvm_emit_bin("sub", l, r, e.as_ptr()),
            Mul(l, r) => self.llvm_emit_bin("mul", l, r, e.as_ptr()),
            And(l, r) => self.llvm_emit_bin("and", l, r, e.as_ptr()),
            Or(l, r) => self.llvm_emit_bin("or", l, r, e.as_ptr()),
            Xor(l, r) => self.llvm_emit_bin("xor", l, r, e.as_ptr()),
            Neg(i) => {
                let ty = self.llvm_type_name();
                let rhs = self.llvm_emit_value(i);
                let var =
                    self.find_existing_var(e.as_ptr()).unwrap_or_else(|| {
                        let v = self.llvm_fresh_var("sub");
                        self.subs.push(CommonSubExpr {
                            ptr: e.as_ptr(),
                            var: v,
                            init: format!("sub {ty} 0, {rhs}"),
                        });
                        v
                    });
                format!("%{}", var)
            },
            Not(i) => {
                let ty = self.llvm_type_name();
                let val = self.llvm_emit_value(i);
                let var =
                    self.find_existing_var(e.as_ptr()).unwrap_or_else(|| {
                        let v = self.llvm_fresh_var("xor");
                        self.subs.push(CommonSubExpr {
                            ptr: e.as_ptr(),
                            var: v,
                            init: format!("xor {ty} {val}, -1"),
                        });
                        v
                    });
                format!("%{}", var)
            },
        }
    }

    fn llvm_emit_bin(
        &mut self,
        op: &str,
        l: &Expr<R>,
        r: &Expr<R>,
        ptr: *const ExprOp<R>,
    ) -> String {
        if let Some(v) = self.find_existing_var(ptr) {
            return format!("%{}", v);
        }
        let ty = self.llvm_type_name();
        let lv = self.llvm_emit_value(l);
        let rv = self.llvm_emit_value(r);
        let v = self.llvm_fresh_var(op);
        self.subs.push(CommonSubExpr {
            ptr,
            var: v,
            init: format!("{op} {ty} {lv}, {rv}"),
        });
        format!("%{}", v)
    }
}

impl<R: Ring> Display for ExprFormatter<'_, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.formatter {
            Formatter::C => {
                let ty = self.r.data_type_name(self.formatter).to_string();
                for sub in self.subs.iter().rev() {
                    writeln!(
                        f,
                        "{:\t>tabs$}{ty} {} = {};\n",
                        "",
                        sub.var,
                        sub.init,
                        tabs = self.tabs
                    )?;
                }
            },
            Formatter::Rust => {
                for sub in self.subs.iter().rev() {
                    writeln!(
                        f,
                        "{:\t>tabs$}let {} = {};\n",
                        "",
                        sub.var,
                        sub.init,
                        tabs = self.tabs
                    )?;
                }
            },
            Formatter::Tex => {
                for sub in self.subs.iter().rev() {
                    writeln!(f, "{} = {}\\\\", sub.var, sub.init)?;
                }
                return write!(f, "{}{}", self.prefix, self.buf);
            },
            Formatter::LLVM => {
                // Generate LLVM code with dependencies listed first
                for sub in self.subs.iter() {
                    writeln!(
                        f,
                        "{:\t>tabs$}%{} = {}",
                        "",
                        sub.var,
                        sub.init,
                        tabs = self.tabs
                    )?;
                }
                return Ok(());
            },
        }

        write!(
            f,
            "{:\t>tabs$}{}{}",
            "",
            self.prefix,
            self.buf,
            tabs = self.tabs
        )
    }
}

pub struct FunctionFormatter<'a, R: Ring, D> {
    vars: Vec<Symbol>,
    inner: D,
    formatter: Formatter,
    function_name: &'a str,
    r: &'a R,
}

impl<'a, R: Ring, D> FunctionFormatter<'a, R, D> {
    pub fn new(
        vars: Vec<Symbol>,
        inner: D,
        formatter: Formatter,
        function_name: &'a str,
        r: &'a R,
    ) -> Self {
        let mut vars = vars;
        vars.sort_by_cached_key(|v| {
            let s = v.as_str();
            (s.len(), s)
        });

        Self { vars, inner, formatter, function_name, r }
    }
}

impl<'a, R: Ring, D: Display> Display for FunctionFormatter<'a, R, D> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self.formatter {
            Formatter::C => {
                let ty = self.r.data_type_name(self.formatter).to_string();
                write!(f, "{ty} {}(", self.function_name)?;

                let mut vars = self.vars.iter();
                if let Some(v) = vars.next() {
                    write!(f, "{ty} {v}")?;

                    for v in vars {
                        write!(f, ", {ty} {v}")?;
                    }
                }

                write!(f, ") {{\n{};\n}}", self.inner)
            },
            Formatter::Rust => {
                let ty = self.r.data_type_name(self.formatter).to_string();
                write!(f, "fn {}(", self.function_name)?;

                let mut vars = self.vars.iter();
                if let Some(v) = vars.next() {
                    write!(f, "{v}: {ty}")?;

                    for v in vars {
                        write!(f, ", {v}: {ty}")?;
                    }
                }

                write!(f, ") -> {ty} {{\n{}\n}}", self.inner)
            },
            Formatter::Tex => self.inner.fmt(f),
            Formatter::LLVM => {
                let ty = self.r.data_type_name(self.formatter).to_string();
                write!(f, "define {ty} @{}(", self.function_name)?;

                let mut vars = self.vars.iter();
                if let Some(v) = vars.next() {
                    write!(f, "{ty} %{v}")?;
                    for v in vars {
                        write!(f, ", {ty} %{v}")?;
                    }
                }

                writeln!(f, ") {{")?;
                writeln!(f, "entry:")?;
                let instructions = format!("{}", self.inner);
                if !instructions.is_empty() {
                    if instructions.ends_with('\n') {
                        write!(f, "{}", instructions)?;
                    } else {
                        writeln!(f, "{}", instructions)?;
                    }
                }
                // Compute final result value
                let last_value = if let Some(sub) =
                    instructions.lines().rev().find_map(|line| {
                        let line = line.trim();
                        if line.starts_with('%') {
                            line.find('=')
                                .map(|eq_pos| line[..eq_pos].trim().to_string())
                        } else {
                            None
                        }
                    }) {
                    sub
                } else {
                    // No instructions; fallback to a constant zero of correct
                    // type
                    "0".to_string()
                };

                writeln!(f, "{:\t>tabs$}ret {ty} {last_value}", "", tabs = 1)?;
                write!(f, "}}")
            },
        }
    }
}
