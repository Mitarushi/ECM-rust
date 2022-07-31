use std::collections::VecDeque;

use ibig::modular::{Modulo, ModuloRing};
use ibig::UBig;
use rug::{Integer, integer::Order};

#[derive(Debug, Clone)]
pub struct Poly<'a> {
    pub a: Vec<Modulo<'a>>,
    pub ring: &'a ModuloRing,
    mod_log: usize,
}

impl<'a> Poly<'a> {
    pub fn new(a: Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Self {
        let mod_log = ring.modulus().bit_len() / 8 + 1;
        Poly { a, ring, mod_log }
    }

    fn to_rug(&self, padding: usize) -> Integer {
        let mut bytes = vec![0; padding * self.a.len()];
        for (i, a) in self.a.iter().enumerate() {
            let a_byte = a.residue().to_le_bytes();
            bytes[i * padding..i * padding + a_byte.len()].copy_from_slice(&a_byte);
        }
        Integer::from_digits(&bytes, Order::Lsf)
    }

    fn from_rug(a: &Integer, padding: usize, ring: &'a ModuloRing) -> Self {
        let a = a.to_digits(Order::Lsf);
        let n = (a.len() + padding - 1) / padding;
        let mut result = Vec::with_capacity(n);

        for i in 0..n {
            let idx = i * padding;
            let x = UBig::from_le_bytes(&a[idx..(idx + padding).min(a.len())]);
            result.push(ring.from(x));
        }

        Poly::new(result, ring)
    }

    pub fn set_len(&mut self, n: usize, ring: &'a ModuloRing) {
        self.a.resize(n, ring.from(0));
    }

    pub fn truncate(&mut self, n: usize) {
        self.a.truncate(n);
    }

    pub fn div_set_len(&mut self, n: usize) {
        if self.a.len() > n {
            self.a = self.a[self.a.len() - n..].to_vec();
        }
    }

    pub fn monic_linear(a: Modulo<'a>, ring: &'a ModuloRing) -> Self {
        Poly::new(vec![a, ring.from(1)], ring)
    }

    pub fn zero(ring: &'a ModuloRing) -> Self {
        Poly::new(vec![], ring)
    }

    fn large_mul(&self, rhs: &Self) -> Self {
        let n = self.a.len() + rhs.a.len() - 1;
        let padding = self.mod_log * 2 + (64 - n.leading_zeros() as usize) / 8 + 1;
        let a = self.to_rug(padding);
        let b = rhs.to_rug(padding);
        let mut result = Poly::from_rug(&(a * b), padding, self.ring);
        result.set_len(n, self.ring);
        result
    }

    fn small_mul(&self, rhs: &Self) -> Self {
        let mut result_a = vec![self.ring.from(0); self.a.len() + rhs.a.len() - 1];
        for (i, a) in self.a.iter().enumerate() {
            for (j, b) in rhs.a.iter().enumerate() {
                result_a[i + j] += a * b;
            }
        }
        Poly::new(result_a, self.ring)
    }

    pub fn reverse(&self) -> Self {
        let mut result = self.clone();
        result.a.reverse();
        result
    }

    pub fn inv(&self, degree: usize) -> Self {
        assert_eq!(self.a[0].residue(), ibig::ubig!(1));

        let mut k = 1;
        let mut g = Poly::new(vec![self.ring.from(1)], self.ring);
        while k < degree {
            k *= 2;
            let mut f = self.clone();
            f.truncate(k);
            let mut fg = &f * &g;
            fg.truncate(k);
            fg = -fg;
            fg.a[0] += self.ring.from(2);
            g = &fg * &g;
            g.truncate(k);
        }
        g
    }
}

impl<'a> std::ops::Mul for &Poly<'a> {
    type Output = Poly<'a>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.a.len().min(rhs.a.len()) <= 16 {
            self.small_mul(rhs)
        } else {
            self.large_mul(rhs)
        }
    }
}

impl<'a> std::ops::Sub for &Poly<'a> {
    type Output = Poly<'a>;
    fn sub(self, rhs: Self) -> Self::Output {
        let n = self.a.len().max(rhs.a.len());
        let mut c = vec![self.ring.from(0); n];
        for i in 0..self.a.len() {
            c[i] += &self.a[i];
        }
        for i in 0..rhs.a.len() {
            c[i] -= &rhs.a[i];
        }
        Poly::new(c, self.ring)
    }
}

impl<'a> std::ops::Neg for Poly<'a> {
    type Output = Poly<'a>;
    fn neg(self) -> Self::Output {
        let mut result = self.clone();
        for a in &mut result.a {
            *a = -a.clone();
        }
        result
    }
}

// pub struct multipoint_evaluation<'a> {
//     pub n: usize,
//     pub mul_table: Vec<Poly<'a>>,
//     ring: &'a ModuloRing,
// }
//
// impl<'a> multipoint_evaluation<'a> {
//     pub fn new(point1: &Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Self {
//         let n = point1.len();
//         assert!(n.is_power_of_two());
//
//         let mut mul_table = vec![Poly::zero(ring); n * 2];
//
//         for i in 0..n {
//             mul_table[i + n] = Poly::monic_linear(point1[i].clone(), ring);
//         }
//
//         for i in (1..n).rev() {
//             mul_table[i] = &mul_table[i * 2] * &mul_table[i * 2 + 1];
//         }
//
//         multipoint_evaluation { n, mul_table, ring }
//     }
//
//     fn prod(&self, a: Vec<Poly<'a>>) -> Poly<'a> {
//         let mut a = a;
//         let mut k = self.n / 2;
//         while k > 0 {
//             for i in 0..k {
//                 a[i] = &a[i * 2] * &a[i * 2 + 1];
//             }
//             k /= 2;
//         }
//         a.into_iter().nth(0).unwrap()
//     }
//
//     pub fn eval(&self, point2: &Vec<Modulo<'a>>) -> Modulo<'a> {
//         let mut up_tree_t = vec![Poly::zero(self.ring); self.n * 2];
//         up_tree_t[1] = self.mul_table[1].inv(self.n);
//
//         None
//     }
// }

pub fn multipoint_evaluation_prod<'a>(point1: &Vec<Modulo<'a>>, point2: &Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Modulo<'a> {
    let n = point1.len();
    let mut mul1 = vec![Poly::zero(ring); n * 2];
    let mut mul2 = vec![Poly::zero(ring); n * 2];

    for i in 0..n {
        mul1[i + n] = Poly::monic_linear(point1[i].clone(), ring);
        mul2[i + n] = Poly::monic_linear(-point2[i].clone(), ring);
    }

    for i in (1..n).rev() {
        mul1[i] = &mul1[2 * i] * &mul1[2 * i + 1];
        mul2[i] = &mul2[2 * i] * &mul2[2 * i + 1];
    }

    let mut inv = vec![Poly::zero(ring); n * 2];

    for i in 0..n {
        inv[i + n] = Poly::monic_linear(point2[i].clone(), ring);
    }

    for i in (1..n).rev() {
        let k = mul2[i].a.len();
        let mut g = &inv[2 * i] * &inv[2 * i + 1];
        g.div_set_len(inv[2 * i].a.len());
        let mut t = &mul2[i] * &g;
        t.div_set_len(k);
        t = -t;
        t.a[k - 1] += ring.from(2);
        t = &t * &g;
        t.div_set_len(k);
        inv[i] = t;
    }

    let mut rem = vec![Poly::zero(ring); n * 2];
    rem[1] = &mul1[1] - &mul2[1];
    rem[1].truncate(n);
    for i in 2..2 * n {
        let k = inv[i].a.len();
        let mut q = rem[i / 2].clone();
        q.div_set_len(k - 1);
        q = &q * &inv[i];
        q.div_set_len(k - 1);
        let mut r = &rem[i / 2] - &(&mul2[i] * &q);
        r.truncate(k - 1);
        rem[i] = r;
    }

    let mut prod = ring.from(1);
    for i in 0..n {
        prod *= &rem[i + n].a[0];
    }
    prod
}