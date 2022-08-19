use std::time::Instant;

use ibig::{ubig, UBig};
use rand::{Rng, SeedableRng, thread_rng};
use rand::rngs::StdRng;
use rayon::prelude::*;

use clap::{App, Arg};

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
    println!("{}", n);
    for (idx, factor) in factors.iter().enumerate() {
        print!("\t{}", if idx == 0 { " = " } else { " x " });
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
    let is_prime = miller_rabin(&n, 100, thread_num, rng);

    if pow > 1 {
        println!("{} th root found", pow);
        print_factors(n, &vec![p.clone(); pow], &vec![is_prime; pow]);
        println!();
    }

    if is_prime {
        vec![p; pow as usize]
    } else {
        println!("factorizing : {}", p);

        let mut attempt_count_local = 0;
        let (sigma, factor) = loop {
            *attempt_count += thread_num;
            attempt_count_local += thread_num;

            let start = Instant::now();
            let t = ecm(&p, b1, b2, d, k, step1_mul, step1_hint, thread_num, rng);
            let elapsed = start.elapsed();

            println!("factorize attempt count : {} (total : {})  time : {:?}", attempt_count_local, *attempt_count, elapsed);

            if let Some(t) = t {
                break t;
            }
        };

        let is_prime = factor.iter().map(|f| miller_rabin(f, 100, thread_num, rng)).collect::<Vec<_>>();

        println!("factors found : ");
        print_factors(&p, &factor, &is_prime);
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
    let args = App::new("ecm-factorize")
        .version("1.0.0")
        .about("factorize a number using elliptic curve method")
        .arg(Arg::new("n")
            .help("the number to factorize")
            .required(true)
            .takes_value(true)
            .short('n')
            .long("n"))
        .arg(Arg::new("b1")
            .help("the bound of primes in the step 1")
            .required(true)
            .takes_value(true)
            .short('b')
            .long("b1"))
        .arg(Arg::new("b2")
            .help("the bound of primes in the step 2")
            .required(true)
            .takes_value(true)
            .short('B')
            .long("b2"))
        .arg(Arg::new("d")
            .help("the chunk size of the step 2")
            .required(true)
            .takes_value(true)
            .short('d')
            .long("d"))
        .arg(Arg::new("k")
            .help("the degree of power used in the step 2")
            .required(true)
            .takes_value(true)
            .short('k')
            .long("k"))
        .arg(Arg::new("thread_num")
            .help("the number of threads to use")
            .required(true)
            .takes_value(true)
            .short('t')
            .long("thread_num"))
        .arg(Arg::new("seed")
            .help("the seed to use")
            .takes_value(true)
            .short('s')
            .long("seed"));

    let matches = args.get_matches();

    let n = matches.value_of("n").unwrap().parse::<UBig>().unwrap();
    let b1 = matches.value_of("b1").unwrap().parse::<u64>().unwrap();
    let b2 = matches.value_of("b2").unwrap().parse::<u64>().unwrap();
    let d = matches.value_of("d").unwrap().parse::<u64>().unwrap();
    let k = matches.value_of("k").unwrap().parse::<u64>().unwrap();
    let thread_num = matches.value_of("thread_num").unwrap().parse::<usize>().unwrap();
    let seed = matches.value_of("seed").unwrap_or("").parse::<u64>().ok();

    let start_time = Instant::now();
    let result = factorize(&n, b1, b2, d, k, thread_num, seed);
    let elapsed = start_time.elapsed();

    println!("result : ");
    print_factors(&n, &result, &vec![true; result.len()]);
    println!("process time: {:?}", elapsed);
}
