use std::mem::swap;
use std::str::FromStr;
use std::time::Instant;

use ibig::{ubig, UBig};
use ibig::modular::{Modulo, ModuloRing};
use rand::{Rng, SeedableRng, thread_rng};
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::elliptic_curve::{EllipticCurve, EllipticPoint};
use crate::poly::MultipointEvaluation;
use crate::utils::{gcd, mod_inv};

mod elliptic_curve;
mod poly;
mod utils;

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

fn ecm_sub(n: &UBig, b1: u64, b2: u64, d: u64, sigma: usize) -> Option<UBig> {
    let ring = ModuloRing::new(&n);
    let primes = eratosthenes(b1);

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

    for p in primes.iter() {
        let p_pow = pow_less_than(*p, b1);
        q = curve.mul(&q, p_pow, &ring);
    }

    let g = n.gcd(&q.z.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
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
    let mut s = curve.mul(&q2, i, &ring);
    let mut t = curve.mul(&q2, i - 2, &ring);

    let mut h = ring.from(1);

    let multi_eval = MultipointEvaluation::new(&a, &ring);

    while i * d < b2 {
        let mut b = Vec::new();
        for _ in 0..k {
            let g = n.gcd(&s.z.residue());
            if &ubig!(1) == &g {
                b.push(s.affine_x(&ring));
            } else if &g != n {
                return Some(g);
            }
            t = curve.add_h(&s, &stride, &t);
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

fn ecm(n: &UBig, b1: u64, b2: u64, d: u64, thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    loop {
        let sigma_vec = (0..thread_num).into_iter().map(|_| rng.gen_range(6..1usize << 63)).collect::<Vec<_>>();

        let result_tmp = sigma_vec.into_par_iter().map(|sigma| {
            ecm_sub(n, b1, b2, d, sigma)
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

fn pollard_rho(n: &UBig, k: u64, trial_num: u64, rng: &mut StdRng) -> Vec<UBig> {
    let ring = ModuloRing::new(&n);

    let mut result = Vec::new();

    fn f<'a>(x: &Modulo<'a>, c: &Modulo<'a>) -> Modulo<'a> {
        x * x + c
    }

    for _ in 0..(trial_num / 50).max(1) {
        let c = ring.from(rng.gen_range(ubig!(1)..n - 1));

        let mut a = ring.from(rng.gen_range(ubig!(1)..n - 1));
        let mut b = a.clone();

        for _ in 0..k {
            a = f(&a, &c);
            b = f(&f(&b, &c), &c);

            let g = n.gcd(&(&a - &b).residue());
            if &ubig!(1) < &g && &g < n {
                result.push(g);
            }
        }
    }

    for _ in 0..trial_num {
        let c = ring.from(rng.gen_range(ubig!(1)..n - 1));

        let mut a = ring.from(rng.gen_range(ubig!(1)..n - 1));
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
    }


    println!("{:?}", result);
    clean_divisor(&result, n)
}

fn miller_rabin(n: &UBig, k: u64, rng: &mut StdRng) -> bool {
    let ring = ModuloRing::new(&n);
    let d = (n - &ubig!(1)).trailing_zeros().unwrap();
    let s = (n - &ubig!(1)) >> d;

    'outer: for _ in 0..k {
        let a = ring.from(rng.gen_range(ubig!(2)..n.clone()));
        let mut x = a.pow(&s);

        if &x.residue() == &ubig!(1) || &(&x.residue() + &ubig!(1)) == n {
            continue;
        }

        for _ in 0..d - 1 {
            x = &x * &x;
            if &(&x.residue() + &ubig!(1)) == n {
                continue 'outer;
            }
        }
        return false;
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

fn factorize_sub(n: &UBig, b1: u64, b2: u64, d: u64, rng: &mut StdRng) -> Vec<UBig> {
    if n == &ubig!(1) {
        return vec![];
    }

    let (p, k) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100, rng) {
        vec![p; k as usize]
    } else {
        let q = ecm(&p, b1, b2, d, 12, rng);

        let mut result = Vec::new();
        for p in q {
            result.append(&mut factorize_sub(&p, b1, b2, d, rng));
            // result.push(p);
        }

        let mut results = result.clone();
        for _ in 0..k - 1 {
            results.append(&mut result.clone());
        }
        results
    }
}

fn factorize(n: &UBig, b1: u64, b2: u64, d: u64, seed: Option<u64>) -> Vec<UBig> {
    let mut rng = StdRng::seed_from_u64(seed.unwrap_or_else(|| thread_rng().gen::<u64>()));

    let (mut result, n) = trial_division(n, 10000);
    let result_pollard = pollard_rho(&n, 100000, 10, &mut rng);

    for i in result_pollard.into_iter() {
        result.append(&mut factorize_sub(&i, b1, b2, d, &mut rng));
    }
    result.sort();
    result
}


fn main() {
    // println!("Hello, world!");
    // println!("{:?}", eratosthenes(100));
    let start_time = Instant::now();
    println!("result: {:?}", factorize(&UBig::from_str("283598282799012588354313727318318100165490374946550831678436461954855068456871761675152071482710347887068874127489").unwrap(),
                                       200000, 10000000, 2310, Some(123456)));
    println!("time: {:?}", start_time.elapsed());

    // println!("{:?}", modinv(&BigInt::from(3456757u64), &BigInt::from(5567544567843u64)));
}
