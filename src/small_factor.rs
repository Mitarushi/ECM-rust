use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};
use rand::Rng;
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::clean_divisor;

pub fn trial_division(n: &UBig, k: u64) -> (Vec<UBig>, UBig) {
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

pub fn pollard_rho(n: &UBig, k: u64, trial_num: u64, thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
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

    clean_divisor(&result, n)
}