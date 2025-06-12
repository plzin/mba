#![feature(duration_millis_float)]

use std::time::{Duration, Instant};
use mba::solver::{solve_modular_via_integer_hnf, solve_via_integer_diagonalize, solve_via_modular_diagonalize};
use mba::rings::BinaryRing;
use mba::{matrix::Matrix, vector::Vector};
use num_bigint::BigUint;
use num_traits::One;
use rand::distr::Distribution as _;
use rand::{rngs::StdRng, SeedableRng, distr::Uniform};

fn main() {
    let args: Vec<_> = std::env::args().collect();

    let bits = match args.get(1) {
        Some(arg) => match arg.parse::<u32>() {
            Ok(bits) => bits,
            Err(_) => {
                eprintln!("Invalid argument: {arg}");
                return;
            },
        },
        None => {
            eprintln!("Usage: {} <bits> <iter>", args[0]);
            return;
        },
    };

    let iter = match args.get(2) {
        Some(arg) => match arg.parse::<usize>() {
            Ok(iter) => iter,
            Err(_) => {
                eprintln!("Invalid argument: {arg}");
                return;
            },
        },
        None => 10000,
    };

    mba::choose_binary_ring!(bench(&r, iter), r = bits);
}

fn bench<R: BinaryRing>(r: &R, iter: usize) {
    let rng = &mut StdRng::seed_from_u64(0);
    let dist = Uniform::new(0, 20).unwrap();

    let mut int_hnf_duration = Duration::ZERO;
    let mut int_diag_duration = Duration::ZERO;
    let mut mod_diag_duration = Duration::ZERO;

    let modulus = BigUint::one() << r.bits();

    for _ in 0..iter {
        let dim = dist.sample(rng);
        let cols = dist.sample(rng);
        let a = Matrix::random(dim, cols, r, rng);
        let b = Vector::random(dim, r, rng);


        {
            let a_int = a.transform(|e| R::to_representative(e).into());
            let b_int = b.transform(|e| R::to_representative(e).into());
            let start = Instant::now();
            std::hint::black_box(solve_modular_via_integer_hnf(
                std::hint::black_box(a_int.view()),
                std::hint::black_box(b_int.view()),
                std::hint::black_box(&modulus),
            ));
            int_hnf_duration += start.elapsed();
        }

        {
            let a = a.clone();
            let b = b.clone();
            let start = Instant::now();
            std::hint::black_box(solve_via_integer_diagonalize(
                std::hint::black_box(a),
                std::hint::black_box(b),
                std::hint::black_box(r),
            ));
            int_diag_duration += start.elapsed();
        }

        {
            let start = Instant::now();
            std::hint::black_box(solve_via_modular_diagonalize(
                std::hint::black_box(a),
                std::hint::black_box(b),
                std::hint::black_box(r),
            ));
            mod_diag_duration += start.elapsed();
        }
    }

    println!("int_hnf_duration: {:.2}ms", int_hnf_duration.as_millis_f64());
    println!("int_diag_duration: {:.2}ms", int_diag_duration.as_millis_f64());
    println!("mod_diag_duration: {:.2}ms", mod_diag_duration.as_millis_f64());
}