use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};
use rand::Rng;
use rand::rngs::StdRng;
use rayon::prelude::*;

use crate::{bit_length, clean_divisor, EllipticCurve, EllipticPoint, mod_inv, MultipointEvaluation};

fn cube<'a>(x: &'a Modulo<'a>) -> Modulo<'a> {
    x * x * x
}

fn ecm_sub(n: &UBig, _b1: u64, b2: u64, d: u64, k: u64, sigma: usize, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>) -> Option<UBig> {
    let ring = ModuloRing::new(&n);

    assert!(sigma >= 6);
    let sigma = ring.from(sigma);
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

    let mut a = Vec::new();
    let mut q2 = q.clone();

    for i in 0..d {
        let g = n.gcd(&q.z.residue());
        if &ubig!(1) == &g {
            a.push(-q.affine_x(&ring));
        } else if &g != n {
            return Some(g);
        }

        if i == d / 2 {
            q2 = q.clone();
        }

        for _ in 0..k {
            q = curve.triple_h(&q);
        }
    }

    let a_len = 1 << bit_length(a.len());
    a.resize(a_len, ring.from(1));

    let mut i = 0;

    let mut h = ring.from(1);

    let multi_eval = MultipointEvaluation::new(&a, &ring);

    while i * d < b2 {
        let mut b = Vec::new();
        for _ in 0..a_len {
            for _ in 0..k / 2 {
                q2 = curve.double_h(&q2);
            }

            let g = n.gcd(&q2.z.residue());
            if &ubig!(1) == &g {
                b.push(q2.affine_x(&ring));
            } else if &g != n {
                return Some(g);
            }

            i += 1;
        }
        b.resize(a_len, ring.from(1));

        h *= multi_eval.eval(&b);
    }

    let g = n.gcd(&h.residue());
    if &ubig!(1) < &g && &g < n {
        return Some(g);
    }
    None
}

pub fn ecm(n: &UBig, b1: u64, b2: u64, d: u64, k: u64, step1_mul: &Vec<u64>, step1_hint: &Vec<u64>,
       thread_num: usize, rng: &mut StdRng) -> Vec<UBig> {
    loop {
        let sigma_vec = (0..thread_num).into_iter().map(|_| rng.gen_range(6..1usize << 63)).collect::<Vec<_>>();

        let result_tmp = sigma_vec.into_par_iter().map(|sigma| {
            ecm_sub(n, b1, b2, d, k, sigma, step1_mul, step1_hint)
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
