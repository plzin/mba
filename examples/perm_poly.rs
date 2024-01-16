use mba::perm_poly;

// This generates an 8-bit permutation polynomial of degree 3 and its inverse.
fn main() {
    // Initialize the zero ideal, which allows us to determine
    // when two polynomial expressions induce the same function.
    // It is needed internally by almost all the functions in `perm_poly`.
    let zi = perm_poly::ZeroIdeal::init(8);

    // Generate a random permutation polynomial of degree 3
    // and its inverse.
    let (p, q) = perm_poly::perm_pair(&zi, 3);

    println!("p(x) = {p}");
    println!("q(x) = {q}");

    // Check that p and q are inverses.
    let pq = perm_poly::compose(&p, &q, &zi).simplified(&zi);
    let qp = perm_poly::compose(&q, &p, &zi).simplified(&zi);

    println!("p(q(x)) = {pq}");
    println!("q(p(x)) = {qp}");
}