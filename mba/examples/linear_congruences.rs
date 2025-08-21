use mba::matrix::Matrix;
use mba::solver;
use mba::vector::Vector;

// Solve a system of linear congruences.
fn main() {
    let solution = solver::solve_modular_via_integer_hnf(
        Matrix::from_array([[3, 5], [4, 2]]).view(),
        Vector::from_slice(&[2.into(), 0.into()]),
        &256u32.into(),
    );

    if solution.is_empty() {
        println!("No solution");
    } else {
        println!("Off: {:?}\nBasis:", solution.offset);
        for b in solution.lattice.basis.rows() {
            println!("{b:?}");
        }
    }
}
