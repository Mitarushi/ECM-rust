use std::mem::swap;
use std::str::FromStr;
use std::time::Instant;

use ibig::{ubig, UBig};
use ibig::modular::{Modulo, ModuloRing};
use rand::{Rng, SeedableRng, thread_rng};
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::addition_chain::compute_optimal_hint;
use crate::elliptic_curve::{EllipticCurve, EllipticPoint};
use crate::poly::MultipointEvaluation;
use crate::utils::{gcd, mod_inv};

mod elliptic_curve;
mod poly;
mod utils;
mod addition_chain;

fn eratosthenes(n: u64) -> Vec<u64> {
    let mut primes = vec![true; n as usize];
    primes[0] = false;
    primes[1] = false;
    for i in 2..n {
        if primes[i as usize] {
            for j in (i * i..n).step_by(i as usize) {
                primes[j as usize] = false;
            }
        }
    }
    primes.iter()
        .enumerate()
        .filter(|&(_, &p)| p)
        .map(|(i, _)| i as u64)
        .collect()
}

fn pow_less_than(p: u64, n: u64) -> u64 {
    let mut r = 1;
    while r * p <= n {
        r *= p;
    }
    r
}

fn cube<'a>(x: &'a Modulo<'a>) -> Modulo<'a> {
    x * x * x
}

fn clean_divisor(a: &Vec<UBig>, n: &UBig) -> Vec<UBig> {
    let mut result = vec![n.clone()];
    result.sort();
    result.dedup();

    for x in a.iter() {
        let mut idx = 0;
        while idx < result.len() {
            let g = x.gcd(&result[idx]);

            if &g != &ubig!(1) && &g != &result[idx] {
                result[idx] /= &g;
                result.push(g);
            } else {
                idx += 1;
            }
        }
    }

    result
}

fn ecm_sub(n: &UBig, b1: u64, b2: u64, d: u64, sigma: usize, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>) -> Option<UBig> {
    let ring = ModuloRing::new(&n);

    assert!(sigma >= 6);
    let sigma = ring.from(sigma);
    // let sigma = ring.from(UBig::from_str("8170945836124664").unwrap());
    println!("curve sigma: {:?}", sigma);
    let u = &sigma * &sigma - ring.from(5);
    let v = &sigma * ring.from(4);
    let c_b = cube(&u) * &v * ring.from(4);
    let tmp = &v - &u;
    let c_a = cube(&tmp) * (&u * ring.from(3) + &v) - (&c_b + &c_b);
    let g = n.gcd(&c_b.residue());
    if &g != &ubig!(1) {
        return if &g == n {
            None
        } else {
            Some(g)
        };
    }
    let c = &c_a * mod_inv(&c_b, &ring);

    let mut q = EllipticPoint::new(cube(&u), cube(&v));
    let curve = EllipticCurve::new(c, &ring);

    for (p, hint) in step1_mul.iter().zip(step1_hint.iter()) {
        q = curve.mul_with_hint(&q, *p, *hint, &ring);
    }

    let g = n.gcd(&q.z.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
    } else if &g == n {
        return None;
    }

    println!("hey!!!");

    let mut s = q.clone();
    let mut t = curve.double_h(&s);

    let mut a = vec![-s.affine_x(&ring)];
    for i in 2..d {
        if gcd(i, d) == 1 {
            let g = n.gcd(&t.z.residue());
            if &ubig!(1) == &g {
                a.push(-t.affine_x(&ring));
            } else if &g != n {
                return Some(g);
            }
        }
        s = curve.add_h(&t, &q, &s);
        swap(&mut s, &mut t);
    }

    let k = 1 << (64 - a.len().leading_zeros());
    a.resize(k, ring.from(1));

    let q2 = t;
    let mut i = b1 / d + 1;

    let stride = curve.double_h(&q2);
    let mut s = curve.mul(&q2, i - 2, &ring);
    let mut t = curve.mul(&q2, i, &ring);

    let mut h = ring.from(1);

    let multi_eval = MultipointEvaluation::new(&a, &ring);

    while i * d < b2 {
        let mut b = Vec::new();
        for _ in 0..k {
            let g = n.gcd(&s.z.residue());
            if &ubig!(1) == &g {
                b.push(t.affine_x(&ring));
            } else if &g != n {
                return Some(g);
            }
            s = curve.add_h(&t, &stride, &s);
            swap(&mut s, &mut t);
            i += 2;
        }
        b.resize(k, ring.from(1));

        h *= multi_eval.eval(&b);
    }

    let g = n.gcd(&h.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
    }
    None
}

fn ecm(n: &UBig, b1: u64, b2: u64, d: u64, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>,
       thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    loop {
        let sigma_vec = (0..thread_num).into_iter().map(|_| rng.gen_range(6..1usize << 63)).collect::<Vec<_>>();

        let result_tmp = sigma_vec.into_par_iter().map(|sigma| {
            ecm_sub(n, b1, b2, d, sigma, step1_mul, step1_hint)
        }).collect::<Vec<_>>();

        let is_end = result_tmp.iter().any(|x| x.is_some());
        if is_end {
            let result = clean_divisor(&result_tmp.into_iter().filter_map(|x| x).collect(), n);
            for x in result.iter() {
                println!("divisor: {:?}", x);
            }
            return result;
        }
    }
}


fn trial_division(n: &UBig, k: u64) -> (Vec<UBig>, UBig) {
    let mut primes = Vec::new();
    let mut n = n.clone();
    for p in 2..k + 1 {
        let p = UBig::from(p);
        while &n % &p == ubig!(0) {
            n /= &p;
            primes.push(p.clone());
        }
    }
    (primes, n)
}

fn pollard_rho(n: &UBig, k: u64, trial_num: u64, thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    let ring = ModuloRing::new(&n);

    let mut result = Vec::new();

    fn f<'a>(x: &Modulo<'a>, c: &Modulo<'a>) -> Modulo<'a> {
        x * x + c
    }

    let mut gen_params = |k: usize| {
        (0..k).into_iter().map(|_| ring.from(rng.gen_range(ubig!(1)..n - 1))).collect::<Vec<_>>()
    };

    let test_run = (trial_num / 10).max(1);
    let test = |a: Modulo, c: Modulo| {
        let mut result = Vec::new();

        let mut a = a;
        let mut b = a.clone();

        for _ in 0..k {
            a = f(&a, &c);
            b = f(&f(&b, &c), &c);

            let g = n.gcd(&(&a - &b).residue());
            if &ubig!(1) < &g && &g < n {
                result.push(g);
            }
        }

        result.sort();
        result.dedup();
        result
    };

    let main = |a: Modulo, c: Modulo| {
        let mut result = Vec::new();

        let mut a = a;
        let mut b = a.clone();

        let mut g = ring.from(1);
        for _ in 0..k {
            a = f(&a, &c);
            b = f(&f(&b, &c), &c);

            g *= &a - &b;
        }

        let g = n.gcd(&g.residue());
        if &ubig!(1) < &g && &g < n {
            result.push(g);
        }

        result.sort();
        result.dedup();
        result
    };

    for i in 0..1 {
        let num = if i == 0 {
            test_run
        } else {
            k
        };

        for _ in 0..(num as usize + thread_num - 1) / thread_num {
            let a = gen_params(thread_num);
            let c = gen_params(thread_num);

            let ac = a.into_iter().zip(c.into_iter()).collect::<Vec<_>>();

            let mut result_tmp = ac.into_par_iter().map(|(a, c)| {
                if i == 0 {
                    test(a, c)
                } else {
                    main(a, c)
                }
            }).collect::<Vec<_>>();

            for i in result_tmp.iter_mut() {
                result.append(i);
            }
        }
    }

    println!("{:?}", result);
    clean_divisor(&result, n)
}

fn miller_rabin(n: &UBig, k: u64, thread_num: usize, rng: &mut StdRng) -> bool {
    let ring = ModuloRing::new(&n);
    let d = (n - &ubig!(1)).trailing_zeros().unwrap();
    let s = (n - &ubig!(1)) >> d;

    let main = |a: Modulo| {
        let mut x = a.pow(&s);

        if &x.residue() == &ubig!(1) || &(&x.residue() + &ubig!(1)) == n {
            return true;
        }

        for _ in 0..d - 1 {
            x = &x * &x;
            if &(&x.residue() + &ubig!(1)) == n {
                return true;
            }
        }
        false
    };

    for _ in 0..(k as usize + thread_num - 1) / thread_num {
        let a_vec = (0..thread_num).into_iter().map(|_| ring.from(rng.gen_range(ubig!(2)..n.clone()))).collect::<Vec<_>>();
        let result = a_vec.into_iter().map(|a| main(a)).collect::<Vec<_>>();

        if result.into_iter().any(|x| !x) {
            return false;
        }
    }
    true
}

fn nth_root(x: &UBig, k: usize) -> UBig {
    let mut y = ubig!(1) << (x.bit_len() / k + 1);
    let mut q = x / &y.pow(k - 1);
    loop {
        y = (&y * UBig::from(k) + &q - &y) / UBig::from(k);
        let z = y.pow(k - 1);
        q = x / &z;
        if &q >= &y {
            return y;
        }
    }
}

fn pow_check(n: &UBig) -> Option<(UBig, usize)> {
    let max_k = n.bit_len();

    for k in (1..max_k).rev() {
        let root = nth_root(n, k);
        if &root.pow(k) == n {
            return Some((root, k));
        }
    }
    None
}

fn factorize_sub(n: &UBig, b1: u64, b2: u64, d: u64, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>,
                 thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    if n == &ubig!(1) {
        return vec![];
    }

    let (p, k) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100, thread_num, rng) {
        vec![p; k as usize]
    } else {
        let q = ecm(&p, b1, b2, d, step1_mul, step1_hint, thread_num, rng);

        let mut result = Vec::new();
        for p in q {
            result.append(&mut factorize_sub(&p, b1, b2, d, &step1_mul, &step1_hint, thread_num, rng));
        }

        let mut results = result.clone();
        for _ in 0..k - 1 {
            results.append(&mut result.clone());
        }
        results
    }
}

fn factorize(n: &UBig, b1: u64, b2: u64, d: u64, thread_num: usize, seed: Option<u64>) -> Vec<UBig> {
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or_else(|| thread_rng().gen::<u64>()));

    let (mut result, n) = trial_division(n, 10000);
    let result_pollard = pollard_rho(&n, 100000, 12, thread_num, &mut rng);

    let prime_less_than_b1 = eratosthenes(b1);
    let step1_mul = prime_less_than_b1.into_iter().map(|p| pow_less_than(p, b1)).collect::<Vec<_>>();
    let step1_hint = step1_mul.iter().map(|x| compute_optimal_hint(x.clone())).collect::<Vec<_>>();

    for i in result_pollard.into_iter() {
        result.append(&mut factorize_sub(&i, b1, b2, d, &step1_mul, &step1_hint, thread_num, &mut rng));
    }
    result.sort();
    result
}


fn main() {
    // println!("Hello, world!");
    // println!("{:?}", eratosthenes(100));
    let start_time = Instant::now();
    println!("result: {:?}", factorize(&UBig::from_str("283598282799012588354313727318318100165490374946550831678436461954855068456871761675152071482710347887068874127489").unwrap(),
                                       200000, 10000000, 2310, 12, Some(12345678)));
    println!("time: {:?}", start_time.elapsed());

    // println!("{:?}", modinv(&BigInt::from(3456757u64), &BigInt::from(5567544567843u64)));
}
