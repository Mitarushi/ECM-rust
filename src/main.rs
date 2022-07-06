mod modint;
mod elliptic_curve;

use std::str::FromStr;
use num_bigint;
use num_bigint::{BigInt, UniformBigInt};
use num_traits::{One, Zero};
use rand::distributions::uniform::UniformSampler;
use crate::elliptic_curve::{EllipticPoint, EllipticCurve};
use crate::modint::Montgomery;

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

fn gcd(a: &BigInt, b: &BigInt) -> BigInt {
    if a.is_zero() {
        return b.clone();
    }
    if b.is_zero() {
        return a.clone();
    }
    let v1 = a.trailing_zeros().unwrap();
    let v2 = b.trailing_zeros().unwrap();
    let beta = std::cmp::min(v1, v2);
    let mut a = a >> v1;
    let mut b = b >> v2;
    while a != b {
        if &a < &b {
            std::mem::swap(&mut a, &mut b);
        }
        let v = &a - &b;
        a = b;
        b = &v >> &v.trailing_zeros().unwrap();
    }
    a << beta
}

fn mod_inv(a: &BigInt, n: &BigInt) -> BigInt {
    if a.is_one() {
        BigInt::one()
    } else {
        n + (-n * mod_inv(&(n % a), a) + 1) / a
    }
}

fn pow_less_than(p: u64, n: u64) -> u64 {
    let mut r = 1;
    while r * p <= n {
        r *= p;
    }
    r
}

fn ecm(n: &BigInt, b1: u64, b2: u64, d: u64) -> BigInt {
    let mr = Montgomery::new(n.clone());

    let mut rng = rand::thread_rng();
    let sampler = UniformBigInt::new(BigInt::from(6u64), BigInt::min(n.clone(), BigInt::from(100000000000000000u64)));

    let primes = eratosthenes(b1);

    loop {
        let sigma = mr.convert_bigint(&sampler.sample(&mut rng));
        // let sigma = mr.convert_bigint(&BigInt::from_str("8689346476060549").unwrap());
        println!("{}", mr.val(&sigma));
        let u = mr.sub(&mr.square(&sigma), &mr.convert_u64(5));
        let v = mr.mul(&sigma, &mr.convert_u64(4));
        let c_b = mr.mul(&mr.pow(&u, 3), &mr.mul(&v, &mr.convert_u64(4)));
        let tmp = mr.pow(&mr.sub(&v, &u), 3);
        let tmp = mr.mul(&tmp, &mr.add(&mr.mul(&u, &mr.convert_u64(3)), &v));
        let c_a = mr.sub(&tmp, &mr.add(&c_b, &c_b));
        let g = gcd(&mr.val(&c_b), &n);
        if !&g.is_one() {
            if &g == n {
                continue;
            } else {
                return g;
            }
        }
        let c = mr.mul(&c_a, &mr.convert_bigint(&mod_inv(&mr.val(&c_b), &n)));

        let mut q = EllipticPoint::new(mr.pow(&u, 3), mr.pow(&v, 3));
        let curve = EllipticCurve::new(c);

        for p in primes.iter() {
            let p_pow = pow_less_than(*p, b1);
            q = curve.mul(&q, p_pow, &mr);
        }

        let g = gcd(&mr.val(&q.z), &n);
        if &BigInt::one() < &g && &g < n {
            return g;
        }

        let mut s = vec![curve.double_h(&q, &mr)];
        s.push(curve.double_h(&s[0], &mr));
        let mut beta = Vec::new();
        for i in 0..d as usize {
            if i > 1 {
                s.push(curve.add_h(&s[i - 1], &s[0], &s[i - 2], &mr));
            }
            beta.push(mr.mul(&s[i].x, &s[i].z));
        }

        let mut h = mr.convert_u64(1);
        let b = b1 - 1;

        let mut t = curve.mul(&q, b - 2 * d, &mr);
        let mut r = curve.mul(&q, b, &mr);

        let mut u = b;
        while u < b2 {
            let alpha = mr.mul(&r.x, &r.z);

            for p2 in segment_sieve(u + 2, 2 * d - 1) {
                let delta = ((p2 - u) / 2 - 1) as usize;
                let tmp = mr.mul(&mr.sub(&r.x, &s[delta].x), &mr.add(&r.z, &s[delta].z));
                let tmp = mr.sub(&tmp, &mr.sub(&alpha, &beta[delta]));
                h = mr.mul(&h, &tmp);
            }

            let tmp = r.clone();
            r = curve.add_h(&r, &s[d as usize - 1], &t, &mr);
            t = tmp;
            u += 2 * d;
        }

        let g = gcd(&mr.val(&h), &n);
        if &BigInt::one() < &g && &g < n {
            return g;
        }
    }
}

fn trial_division(n: &BigInt, k: u64) -> (Vec<BigInt>, BigInt) {
    let mut primes = Vec::new();
    let mut n = n.clone();
    for p in 2..k + 1 {
        while &n % &BigInt::from(p) == BigInt::zero() {
            n /= &BigInt::from(p);
            primes.push(BigInt::from(p));
        }
    }
    (primes, n)
}

fn miller_rabin(n: &BigInt, k: u64) -> bool {
    let mr = Montgomery::new(n.clone());
    let d = (n - &BigInt::one()).trailing_zeros().unwrap();
    let s = (n - &BigInt::one()) >> d;

    let mut rng = rand::thread_rng();
    let sampler = UniformBigInt::new(BigInt::from(2u64), n.clone());

    'outer: for _ in 0..k {
        let a = mr.convert_bigint(&sampler.sample(&mut rng));
        let mut x = mr.binary_pow(&a, &s);

        if &mr.val(&x) == &BigInt::one()  || &(&mr.val(&x) + &BigInt::one()) == n{
            continue;
        }

        for _ in 0..d - 1 {
            x = mr.mul(&x, &x);
            if &(&mr.val(&x) + &BigInt::one()) == n {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn pow_check(n: &BigInt) -> Option<(BigInt, u64)> {
    let max_k = n.bits();

    for k in (1..max_k).rev() {
        let root = n.nth_root(k as u32);
        if &root.pow(k as u32) == n {
            return Some((root, k));
        }
    }
    None
}

fn factorize_sub(n: &BigInt, b1:u64, b2: u64, d: u64) -> Vec<BigInt> {
    if n.is_one() {
        return vec![];
    }

    let (p, k) = pow_check(&n).unwrap();
    if miller_rabin(&p, 100) {
        vec![p; k as usize]
    } else {
        let q = ecm(&p, b1, b2, d);
        println!("{}", q);
        let mut result = factorize_sub(&q, b1, b2, d);
        result.append(&mut factorize_sub(&(n / &q), b1, b2, d));

        let mut results = result.clone();
        for _ in 0..k - 1 {
            results.append(&mut result.clone());
        }
        results
    }
}

fn factorize(n :&BigInt, b1: u64, b2: u64, d: u64) -> Vec<BigInt> {
    let (mut result, mut n) = trial_division(n, 10000);
    result.append(&mut factorize_sub(&n, b1, b2, d));
    result.sort();
    result
}


fn main() {
    println!("Hello, world!");
    println!("{:?}", eratosthenes(100));
    println!("{:?}", factorize(&BigInt::from_str("627057063764139831929324851379409869378845668175598843037877190478889006888518431438644711527536922839520331484815861906173161536477065546885468336421475511783984145060592245840032548652210559519683510271").unwrap(),
                               100000000, 1000000000, 10000000));
    // println!("{:?}", modinv(&BigInt::from(3456757u64), &BigInt::from(5567544567843u64)));
}
