use mba::{
    choose_binary_ring,
    lattice::AffineLattice,
    matrix::Matrix,
    rings::{BinaryRing, RingElement as _},
    solver::{modular_diagonalize, solve_scalar_congruence},
    tex,
    vector::Vector,
};
use wasm_bindgen::prelude::*;

/// Stores the intermediate results during the computation.
#[wasm_bindgen(getter_with_clone)]
pub struct SolveTrace {
    /// The diagonalization of the matrix.
    pub diag: String,

    /// The resulting diagonal system.
    #[wasm_bindgen(js_name = "scalarSystem")]
    pub scalar_system: String,

    /// Linear congruences and solutions.
    #[wasm_bindgen(js_name = "linearSolutions")]
    pub linear_solutions: String,

    /// Vector form of the solution.
    #[wasm_bindgen(js_name = "vectorSolution")]
    pub vector_solution: String,

    /// The final solution.
    #[wasm_bindgen(js_name = "finalSolution")]
    pub final_solution: String,
}

#[wasm_bindgen(js_name = "solveLinearSystem")]
pub fn solve_linear_system(
    #[wasm_bindgen(js_name = "matrixString")] matrix_string: String,
    bits: u32,
) -> Result<SolveTrace, String> {
    // Parse the matrix here so the code isn't duplicated for each ring.
    let mut a = Vec::new();
    let mut b = Vec::new();

    let mut rows = 0;
    let mut cols = 0;

    // Each line contains a row.
    for line in matrix_string.lines() {
        // Each element is a space-separated string.
        for element in line.split_ascii_whitespace() {
            a.push(element);
        }

        let Some(rhs) = a.pop() else {
            return Err("The system must have a right-hand side".into());
        };

        b.push(rhs);

        // The first row determines the number of columns.
        if rows == 0 {
            cols = a.len();
            if cols == 0 {
                return Err("The system must have at least one columns".into());
            }
        } else {
            // Check that each row has the same number of entries.
            let row_entries = a.len() - rows * cols;
            if row_entries != cols {
                return Err(format!(
                    "Row {} has a different number of entries ({row_entries}) \
                     than the first row ({cols})",
                    rows + 1
                ));
            }
        }

        rows += 1;
    }

    assert_eq!(rows * cols, a.len());
    assert_eq!(rows, b.len());

    choose_binary_ring!(
        solve_linear_system_impl(a, b, rows, cols, &r),
        r = bits
    )
}

fn solve_linear_system_impl<R: BinaryRing>(
    m: Vec<&str>,
    n: Vec<&str>,
    rows: usize,
    cols: usize,
    r: &R,
) -> Result<SolveTrace, String> {
    // Parse the entries.
    let mut a = Matrix::zero(rows, cols);
    let mut b = Vector::zero(rows);

    for i in 0..rows {
        for j in 0..cols {
            a[(i, j)] = match r.element_from_string(m[i * cols + j]) {
                Some(e) => e,
                None => {
                    return Err(format!(
                        "Failed to parse element at ({i}, {j})"
                    ));
                },
            };
        }
    }

    for i in 0..rows {
        b[i] = match r.element_from_string(n[i]) {
            Some(e) => e,
            None => {
                return Err(format!(
                    "Failed to parse element at ({i}, {cols})"
                ));
            },
        };
    }

    // Diagonalize the matrix.
    let mut d = a.clone();
    let (s, t) = modular_diagonalize(d.view_mut(), r);
    let diag = format!(
        "{}={}{}{}",
        tex::underbrace(d.to_tex(), tex::bold("D")),
        tex::underbrace(s.to_tex(), tex::bold("S")),
        tex::underbrace(a.to_tex(), tex::bold("A")),
        tex::underbrace(t.to_tex(), tex::bold("T")),
    );

    // Transform the right-hand side.
    let b_new = s.mul_vec_post(&b, r);

    let scalar_system = format!(
        "{}\\mathbf{{x'}}={}{}={}",
        tex::underbrace(d.to_tex(), tex::bold("D")),
        tex::underbrace(s.to_tex(), tex::bold("S")),
        tex::underbrace(b.to_tex(), tex::bold("b")),
        b_new.to_tex()
    );

    let b = b_new;

    let min_dim = d.min_dim();

    // Check if any row trivially has no solution.
    if let Some(i) = b.iter().skip(min_dim).position(|e| !e.is_zero()) {
        let i = min_dim + i;
        let linear_solutions = format!(
            "\\text{{Row {}: }} 0={}\\implies \\text{{No solution}}",
            i + 1,
            b[i]
        );

        return Ok(SolveTrace {
            diag,
            scalar_system,
            linear_solutions,
            vector_solution: String::new(),
            final_solution: String::new(),
        });
    }

    // Some solution to the system.
    let mut offset = Vector::zero(d.num_cols());

    // The basis of the kernel.
    let mut basis = Matrix::zero(0, d.num_cols());

    let mut linear_solutions = "\\begin{align}".to_owned();

    // Variable index for basis.
    let mut j = 1;

    // Solve the scalar linear congruences.
    for i in 0..min_dim {
        let a = &d[(i, i)];
        let b = &b[i];

        linear_solutions +=
            &format!("{}x'_{{{}}}&={} &\\implies ", a, i + 1, b);

        let Some((x, kern)) = solve_scalar_congruence(a, b, r) else {
            linear_solutions += "\\text{No solution!}&\\end{align}";
            return Ok(SolveTrace {
                diag,
                scalar_system,
                linear_solutions,
                vector_solution: String::new(),
                final_solution: String::new(),
            });
        };

        if kern.is_zero() {
            linear_solutions += &format!("x'_{{{}}}&={}\\\\", i + 1, x);
        } else {
            linear_solutions +=
                &format!("x'_{{{}}}&={}+{}a_{{{}}}\\\\", i + 1, x, kern, j);
            j += 1;
        }

        // The particular solution is an entry is
        // the particular solution of the whole system.
        offset[i] = x;

        // If the kernel is zero, then the vector is zero for sure.
        if !kern.is_zero() {
            let row = basis.num_rows();
            basis.append_zero_rows(1);
            basis[(row, i)] = kern;
        }
    }

    // If there are more variables then equations
    // then there are no restrictions on the variables
    // from index d.rows
    for i in d.num_rows()..d.num_cols() {
        let row = basis.num_rows();
        basis.append_zero_rows(1);
        basis[(row, i)] = R::one();

        linear_solutions += &format!("&&x'_{{{}}}&=a_{{{}}}\\\\", i + 1, j);
        j += 1;
    }

    linear_solutions += "\\end{align}";

    let mut solution = AffineLattice::from_offset_basis(offset, basis);

    let x_old = tex::parens(solution.to_tex()).to_string();
    let vector_solution = format!("\\mathbf{{x'}}={x_old}");

    solution.offset = t.mul_vec_post(&solution.offset, r);
    solution.lattice.basis = solution.lattice.basis.mul(t.transpose(), r);

    let final_solution = format!(
        "\\mathbf{{x}}={}{}={}",
        tex::underbrace(t.to_tex(), tex::bold("T")),
        tex::underbrace(x_old, tex::bold("x'")),
        solution.to_tex()
    );

    Ok(SolveTrace {
        diag,
        scalar_system,
        linear_solutions,
        vector_solution,
        final_solution,
    })
}
