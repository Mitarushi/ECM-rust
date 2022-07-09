use std::str::FromStr;
use std::time::Instant;
use ibig::modular::{ModuloRing, Modulo};
use ibig::{ubig, UBig};
use rand::{thread_rng, Rng};
use rayon::prelude::*;

mod elliptic_curve;

use crate::elliptic_curve::{EllipticPoint, EllipticCurve};

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

fn segment_sieve(n: u64, d: u64) -> Vec<u64> {
    let mut sieve = vec![true; d as usize];

    let mut i = 2;
    while i * i <= n + d {
        let from = (n + i - 1) / i * i;
        for j in (from..n + d).step_by(i as usize) {
            sieve[(j - n) as usize] = false;
        }
        i += 1;
    }
    sieve.iter()
        .enumerate()
        .filter(|&(_, &x)| x)
        .map(|(i, _)| i as u64 + n)
        .collect()
}

fn mod_inv(a: &UBig, n: &UBig) -> UBig {
    if a == &ubig!(1) {
        ubig!(1)
    } else {
        n - (n * mod_inv(&(n % a), a)) / a
    }
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

fn ecm_sub(n: &UBig, b1: u64, b2: u64, d: u64) -> Option<UBig> {
    let ring = ModuloRing::new(&n);
    let primes = eratosthenes(b1);

    let sigma = ring.from(thread_rng().gen_range(ubig!(6)..n.clone()));
    // let sigma = ring.from(UBig::from_str("8689346476060549").unwrap());
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
        }
    }
    let c = &c_a * ring.from(&mod_inv(&c_b.residue(), &n));

    let mut q = EllipticPoint::new(cube(&u), cube(&v));
    let curve = EllipticCurve::new(c);

    for p in primes.iter() {
        let p_pow = pow_less_than(*p, b1);
        q = curve.mul(&q, p_pow, &ring);
    }

    let g = n.gcd(&q.z.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
    }

    let mut s = vec![curve.double_h(&q)];
    s.push(curve.double_h(&s[0]));
    let mut beta = Vec::new();
    for i in 0..d as usize {
        if i > 1 {
            s.push(curve.add_h(&s[i - 1], &s[0], &s[i - 2]));
        }
        beta.push(&s[i].x * &s[i].z);
    }

    let mut h = ring.from(1);
    let b = b1 - 1;

    let mut t = curve.mul(&q, b - 2 * d, &ring);
    let mut r = curve.mul(&q, b, &ring);

    let mut u = b;
    while u < b2 {
        let alpha = &r.x * &r.z;

        for p2 in segment_sieve(u + 2, 2 * d - 1) {
            let delta = ((p2 - u) / 2 - 1) as usize;
            h *= (&r.x - &s[delta].x) * (&r.z + &s[delta].z) - &alpha + &beta[delta];
        }

        let tmp = r.clone();
        r = curve.add_h(&r, &s[d as usize - 1], &t);
        t = tmp;
        u += 2 * d;
    }

    let g = n.gcd(&h.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
    }
    None
}

fn ecm(n: &UBig, b1: u64, b2: u64, d: u64, k: usize) -> Vec<UBig> {
    loop {
        let result_tmp = (0..k).into_par_iter().map(|_| {
            ecm_sub(n, b1, b2, d)
        }).collect::<Vec<_>>();

        let is_end = result_tmp.iter().any(|x| x.is_some());
        if is_end {
            let mut result = Vec::new();
            for r in result_tmp {
                if let Some(mut r) = r {
                    println!("divisor: {:?}", r);
                    if result.is_empty() {
                        result.push(r.clone());
                        result.push(n / &r)
                    } else {
                        let mut next_result = Vec::new();
                        for  s in result.iter() {
                            let g = s.gcd(&r);
                            if &g != &ubig!(1) {
                                if &g != s {
                                    next_result.push(s / &g);
                                    next_result.push(g.clone());
                                } else {
                                    next_result.push(s.clone());
                                }
                                r /= &g;
                            } else {
                                next_result.push(s.clone());
                            }
                        }
                        result.clear();
                        result.append(&mut next_result);
                    }
                }
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

fn miller_rabin(n: &UBig, k: u64) -> bool {
    let ring = ModuloRing::new(&n);
    let d = (n - &ubig!(1)).trailing_zeros().unwrap();
    let s = (n - &ubig!(1)) >> d;

    'outer: for _ in 0..k {
        let a = ring.from(thread_rng().gen_range(ubig!(2)..n.clone()));
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

fn factorize_sub(n: &UBig, b1: u64, b2: u64, d: u64) -> Vec<UBig> {
    if n == &ubig!(1) {
        return vec![];
    }

    let (p, k) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100) {
        vec![p; k as usize]
    } else {
        let q = ecm(&p, b1, b2, d, 48);
        let mut result = Vec::new();
        for p in q {
            result.append(&mut factorize_sub(&p, b1, b2, d));
            // result.push(p);
        }

        let mut results = result.clone();
        for _ in 0..k - 1 {
            results.append(&mut result.clone());
        }
        results
    }
}

fn factorize(n: &UBig, b1: u64, b2: u64, d: u64) -> Vec<UBig> {
    let (mut result, n) = trial_division(n, 10000);
    result.append(&mut factorize_sub(&n, b1, b2, d));
    result.sort();
    result
}


fn main() {
    // println!("Hello, world!");
    // println!("{:?}", eratosthenes(100));
    let start_time = Instant::now();
    println!("result: {:?}", factorize(&UBig::from_str("627057063764139831929324851379409869378845668175598843037877190478889006888518431438644711527536922839520331484815861906173161536477065546885468336421475511783984145060592245840032548652210559519683510271").unwrap(),
                                       10000000, 1000000000, 100000));
    println!("time: {:?}", start_time.elapsed());

    // println!("{:?}", modinv(&BigInt::from(3456757u64), &BigInt::from(5567544567843u64)));
}
