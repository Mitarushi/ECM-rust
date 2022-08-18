use std::str::FromStr;
use std::time::Instant;

use ibig::{ubig, UBig};
use rand::{Rng, SeedableRng, thread_rng};
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::addition_chain::compute_optimal_hint;
use crate::ecm::ecm;
use crate::elliptic_curve::{EllipticCurve, EllipticPoint};
use crate::miller_rabin::miller_rabin;
use crate::poly::MultipointEvaluation;
use crate::small_factor::{pollard_rho, trial_division};
use crate::utils::{bit_length, clean_divisor, eratosthenes, mod_inv, pow_check, pow_less_than};

mod elliptic_curve;
mod poly;
mod utils;
mod addition_chain;
mod miller_rabin;
mod small_factor;
mod ecm;

fn factorize_sub(n: &UBig, b1: u64, b2: u64, d: u64, k: u64, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>,
                 thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    if n == &ubig!(1) {
        return vec![];
    }

    let (p, pow) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100, thread_num, rng) {
        vec![p; pow as usize]
    } else {
        let q = ecm(&p, b1, b2, d, k, step1_mul, step1_hint, thread_num, rng);

        let mut result = Vec::new();
        for p in q {
            result.append(&mut factorize_sub(&p, b1, b2, d, k, &step1_mul, &step1_hint, thread_num, rng));
        }

        let mut results = result.clone();
        for _ in 0..pow - 1 {
            results.append(&mut result.clone());
        }
        results
    }
}

fn factorize(n: &UBig, b1: u64, b2: u64, d: u64, k: u64, thread_num: usize, seed: Option<u64>) -> Vec<UBig> {
    let seed = seed.unwrap_or_else(|| thread_rng().gen::<u64>());
    println!("seed: {}", seed);
    let mut rng = StdRng::seed_from_u64(seed);

    let (mut result, n) = trial_division(n, 10000);
    let result_pollard = pollard_rho(&n, 100000, 12, thread_num, &mut rng);

    let prime_less_than_b1 = eratosthenes(b1);
    let step1_mul = prime_less_than_b1.into_iter().map(|p| pow_less_than(p, b1)).collect::<Vec<_>>();
    let step1_hint = step1_mul.clone().into_par_iter().map(|x| compute_optimal_hint(x)).collect::<Vec<_>>();

    for i in result_pollard.into_iter() {
        result.append(&mut factorize_sub(&i, b1, b2, d, k, &step1_mul, &step1_hint, thread_num, &mut rng));
    }
    result.sort();
    result
}


fn main() {
    // println!("Hello, world!");
    // println!("{:?}", eratosthenes(100));
    let start_time = Instant::now();
    println!("result: {:?}", factorize(&UBig::from_str("283598282799012588354313727318318100165490374946550831678436461954855068456871761675152071482710347887068874127489").unwrap(),
                                       200000, 10000000, 2048, 12, 12, Some(1234)));
    println!("time: {:?}", start_time.elapsed());

    // println!("{:?}", modinv(&BigInt::from(3456757u64), &BigInt::from(5567544567843u64)));
}
