//! Generates LaTeX code for various objects.

use crate::lattice::{AffineLattice, Lattice};
use crate::matrix::{Matrix, MatrixStorage};
use crate::poly::Poly;
use crate::rings::{Ring, RingElement};
use crate::vector::{Vector, VectorStorage};


pub fn bold<T>(t: T) -> TexBold<T> {
    TexBold(t)
}

pub struct TexBold<T>(T);

impl<T: std::fmt::Display> std::fmt::Display for TexBold<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\\mathbf{{{}}}", self.0)
    }
}

pub fn parens<T>(t: T) -> TexParens<T> {
    TexParens(t)
}

pub struct TexParens<T>(T);


impl<T: std::fmt::Display> std::fmt::Display for TexParens<T> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\\left({}\\right)", self.0)
    }
}

pub fn underbrace<T, U>(inner: T, label: U) -> TexUnderbrace<T, U> {
    TexUnderbrace { inner, label }
}

pub struct TexUnderbrace<T, U> {
    inner: T,
    label: U,
}

impl<T, U> std::fmt::Display for TexUnderbrace<T, U>
where
    T: std::fmt::Display,
    U: std::fmt::Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\\underbrace{{{}}}_{{{}}}", self.inner, self.label)
    }
}

impl<R: Ring> Poly<R> {
    pub fn to_tex(&self) -> TexPoly<R> {
        TexPoly { poly: self }
    }
}

pub struct TexPoly<'a, R: Ring> {
    poly: &'a Poly<R>,
}

impl<'a, R: Ring> std::fmt::Display for TexPoly<'a, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.poly.coeffs.iter()
            .enumerate()
            .rev()
            .filter(|(_, e)| !e.is_zero());

        fn write_term<E: RingElement>(
            f: &mut std::fmt::Formatter<'_>,
            e: usize,
            c: &E
        ) -> std::fmt::Result {
            if e == 0 {
                return write!(f, "{c}");
            }

            if !c.is_one() {
                write!(f, "{c}")?;
            }

            write!(f, "X")?;

            if e > 1 {
                write!(f, "^{{{e}}}")?;
            }

            Ok(())
        }

        match iter.next() {
            None => return write!(f, "0"),
            Some((i, e)) => write_term(f, i, e)?,
        }

        for (i, e) in iter {
            write!(f, "+")?;
            write_term(f, i, e)?;
        }

        Ok(())
    }
}



impl<R: Ring, S: VectorStorage<R> + ?Sized> Vector<R, S> {
    pub fn to_tex(&self) -> TexVector<R, S> {
        TexVector { vector: self }
    }
}

pub struct TexVector<'a, R: Ring, S: VectorStorage<R> + ?Sized> {
    vector: &'a Vector<R, S>,
}

impl<'a, R: Ring, S: VectorStorage<R> + ?Sized> std::fmt::Display for TexVector<'a, R, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\\left[\\begin{{array}}{{}}")?;
        for e in self.vector.iter() {
            write!(f, "{e}\\\\")?;
        }
        write!(f, "\\end{{array}}\\right]")
    }
}

impl<R: Ring, S: MatrixStorage<R> + ?Sized> Matrix<R, S> {
    pub fn to_tex(&self) -> TexMatrix<R, S> {
        TexMatrix { matrix: self }
    }
}

pub struct TexMatrix<'a, R: Ring, S: MatrixStorage<R> + ?Sized> {
    matrix: &'a Matrix<R, S>,
}

impl<'a, R: Ring, S: MatrixStorage<R> + ?Sized> std::fmt::Display for TexMatrix<'a, R, S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "\\left[\\begin{{array}}{{}}")?;
        for row in self.matrix.rows() {
            let mut iter = row.iter();
            if let Some(e) = iter.next() {
                write!(f, "{e}")?;
            }

            for e in iter {
                write!(f, " & {e}")?;
            }

            write!(f, "\\\\")?;
        }
        write!(f, "\\end{{array}}\\right]")
    }
}

impl<R: Ring> Lattice<R> {
    pub fn to_tex(&self) -> TexLattice<R> {
        TexLattice { lattice: self }
    }
}

pub struct TexLattice<'a, R: Ring> {
    lattice: &'a Lattice<R>,
}

impl<'a, R: Ring> std::fmt::Display for TexLattice<'a, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut iter = self.lattice.basis.rows().enumerate();
        if let Some((i, e)) = iter.next() {
            write!(f, "a_{}{}", i + 1, e.to_tex())?;
        }

        for (i, e) in iter {
            write!(f, "+a_{}{}", i + 1, e.to_tex())?;
        }

        Ok(())
    }
}


impl<R: Ring> AffineLattice<R> {
    pub fn to_tex(&self) -> TexAffineLattice<R> {
        TexAffineLattice { lattice: self }
    }
}

pub struct TexAffineLattice<'a, R: Ring> {
    lattice: &'a AffineLattice<R>,
}

impl<'a, R: Ring> std::fmt::Display for TexAffineLattice<'a, R> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.lattice.offset.to_tex())?;

        for (i, e) in self.lattice.lattice.basis.rows().enumerate() {
            write!(f, "+a_{}{}", i + 1, e.to_tex())?;
        }

        Ok(())
    }
}