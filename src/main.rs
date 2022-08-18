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


fn print_factors(n: &UBig, factors: &Vec<UBig>, is_prime: &Vec<bool>) {
    println!("{} =", n);
    for (idx, factor) in factors.iter().enumerate() {
        print!("\t{}", if idx == 0 { "   " } else { " × " });
        let f = factor.to_string();
        print!("{} ({} digits", f, f.len());
        print!(", {}", if is_prime[idx] { "prime" } else { "composite" });
        println!(")");
    }
}


fn factorize_sub(n: &UBig, b1: u64, b2: u64, d: u64, k: u64, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>,
                 thread_num: usize, rng: &mut StdRng, attempt_count: &mut usize) -> Vec<UBig> {
    if n == &ubig!(1) {
        return vec![];
    }

    let (p, pow) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100, thread_num, rng) {
        vec![p; pow as usize]
    } else {
        println!("factorizing : {}", n);

        let mut attempt_count_local = 0;
        let (sigma, factor) = loop {
            *attempt_count += thread_num;
            attempt_count_local += thread_num;

            let start = Instant::now();
            let t = ecm(&p, b1, b2, d, k, step1_mul, step1_hint, thread_num, rng);
            let elapsed = start.elapsed();

            println!("factorize attempt count : {} (total : {})  time: {:?}", attempt_count_local, *attempt_count, elapsed);

            if let Some(t) = t {
                break t;
            }
        };

        let is_prime = factor.iter().map(|f| miller_rabin(f, 100, thread_num, rng)).collect::<Vec<_>>();

        println!("factors found : ");
        print_factors(n, &factor, &is_prime);
        println!("sigma : {:?}", sigma);
        println!();

        let mut result = Vec::new();
        for (p, is_prime) in factor.iter().zip(is_prime.iter()) {
            if *is_prime {
                result.push(p.clone());
            } else {
                result.append(&mut factorize_sub(p, b1, b2, d, k, step1_mul, step1_hint, thread_num, rng, attempt_count));
            }
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
    println!("seed : {}", seed);
    println!();

    let mut rng = StdRng::seed_from_u64(seed);

    let (mut result, n_remain) = trial_division(n, 10000);
    let result_pollard = pollard_rho(&n_remain, 100000, 12, thread_num, &mut rng);

    let result_print = vec![result.clone(), result_pollard.clone()].concat();
    let result_print = clean_divisor(&result_print, n);
    if result_print.len() > 1 {
        let is_prime = result_print.iter().map(|f| miller_rabin(f, 100, thread_num, &mut rng)).collect::<Vec<_>>();
        println!("small factors found : ");
        print_factors(&n, &result_print, &is_prime);
        println!();
    }

    let prime_less_than_b1 = eratosthenes(b1);
    let step1_mul = prime_less_than_b1.into_iter().map(|p| pow_less_than(p, b1)).collect::<Vec<_>>();
    let step1_hint = step1_mul.clone().into_par_iter().map(|x| compute_optimal_hint(x)).collect::<Vec<_>>();

    let mut attempt_count = 0;

    for i in result_pollard.into_iter() {
        result.append(&mut factorize_sub(&i, b1, b2, d, k, &step1_mul, &step1_hint, thread_num, &mut rng, &mut attempt_count));
    }
    result.sort();
    result
}


fn main() {
    let start_time = Instant::now();

    let n = UBig::from_str("283598282799012588354313727318318100165490374946550831678436461954855068456871761675152071482710347887068874127489").unwrap();

    let result = factorize(&n, 200000, 10000000, 2048, 12, 12, Some(1234));

    println!("result : ");
    print_factors(&n, &result, &vec![true; result.len()]);
    println!("process time: {:?}", start_time.elapsed());
}
