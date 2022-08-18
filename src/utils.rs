use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};

pub fn mod_inv<'a>(a: &Modulo<'a>, ring: &'a ModuloRing) -> Modulo<'a> {
    ring.from(mod_inv_ubig(&a.residue(), &ring.modulus()))
}

fn mod_inv_ubig(a: &UBig, n: &UBig) -> UBig {
    if a == &ubig!(1) {
        ubig!(1)
    } else if a == &ubig!(0) {
        ubig!(0)
    } else {
        n - (n * mod_inv_ubig(&(n % a), a)) / a
    }
}

pub fn bit_length(n: usize) -> u32 {
    64 - (n - 1).leading_zeros()
}

pub fn eratosthenes(n: u64) -> Vec<u64> {
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

pub fn pow_less_than(p: u64, n: u64) -> u64 {
    let mut r = 1;
    while r * p <= n {
        r *= p;
    }
    r
}

pub fn clean_divisor(a: &Vec<UBig>, n: &UBig) -> Vec<UBig> {
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

pub fn pow_check(n: &UBig) -> Option<(UBig, usize)> {
    let max_k = n.bit_len();

    for k in (1..max_k).rev() {
        let root = nth_root(n, k);
        if &root.pow(k) == n {
            return Some((root, k));
        }
    }
    None
}