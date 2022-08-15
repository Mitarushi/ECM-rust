use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};
use rug::{Integer, integer::Order};

use crate::mod_inv;
use crate::utils::bit_length;

#[derive(Debug, Clone)]
pub struct Poly<'a> {
    pub a: Vec<Modulo<'a>>,
    pub ring: &'a ModuloRing,
    mod_log: usize,
}

impl<'a> Poly<'a> {
    pub fn new(a: Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Self {
        let mod_log = ring.modulus().bit_len();
        Poly { a, ring, mod_log }
    }

    fn to_rug(&self, padding: usize) -> Integer {
        let mut bytes = vec![0; padding * self.len()];
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

    pub fn div_set_len(&mut self, n: usize) {
        if self.len() > n {
            self.a = self.a[self.len() - n..].to_vec();
        }
    }

    pub fn monic_linear(a: Modulo<'a>, ring: &'a ModuloRing) -> Self {
        Poly::new(vec![a, ring.from(1)], ring)
    }

    pub fn zero(ring: &'a ModuloRing) -> Self {
        Poly::new(vec![], ring)
    }

    pub fn one(ring: &'a ModuloRing) -> Self {
        Poly::new(vec![ring.from(1)], ring)
    }

    fn large_mul(&self, rhs: &Self) -> Self {
        let n = self.len() + rhs.len() - 1;
        let max_plus = self.len().min(rhs.len());
        let padding = (self.mod_log * 2 + bit_length(max_plus) as usize) / 8 + 1;
        let a = self.to_rug(padding);
        let b = rhs.to_rug(padding);
        let mut result = Poly::from_rug(&(a * b), padding, self.ring);
        result.set_len(n, self.ring);
        result
    }

    fn small_mul(&self, rhs: &Self) -> Self {
        let mut result_a = vec![ubig!(0); self.len() + rhs.len() - 1];
        for (i, a) in self.a.iter().enumerate() {
            for (j, b) in rhs.a.iter().enumerate() {
                result_a[i + j] += a.residue() * b.residue();
            }
        }
        Poly::new(result_a.into_iter().map(|x| self.ring.from(x)).collect(), self.ring)
    }

    pub fn reverse(&self) -> Self {
        let mut result = self.clone();
        result.a.reverse();
        result
    }

    pub fn inv(&self, degree: usize) -> Self {
        let mut k = 1;
        let mut g = Poly::new(vec![mod_inv(self.a.last().unwrap(), self.ring)], self.ring);
        while k < degree {
            let mut f = self.clone();
            f.div_set_len(k * 2);
            let mut fg = &f * &g;
            fg.div_set_len(k * 2);
            fg.set_len(k, self.ring);
            let mut delta_g = &fg * &g;
            delta_g.div_set_len(k);
            g.shift_mul(k);
            g -= &delta_g;
            k *= 2;
        }
        g.div_set_len(degree);
        g
    }

    pub fn small_mul_t(&self, x: &Self) -> Self {
        let c = self.len();
        let a = x.len();
        let b = c + 1 - a;

        let mut result = vec![self.ring.from(0); b];
        for i in 0..b {
            for (j, p) in x.a.iter().enumerate() {
                result[i] += &self.a[i + j] * p;
            }
        }

        Poly::new(result, self.ring)
    }

    pub fn dual_mul_t(&self, x1: &Self, x2: &Self) -> (Self, Self) {
        assert_eq!(x1.len(), x2.len());
        let c = self.len();
        let a = x1.len();
        let b = c + 1 - a;

        if a < 32 {
            return (self.small_mul_t(x1), self.small_mul_t(x2));
        }

        let mut x = x2.clone();
        x.set_len(c, self.ring);
        x.a.append(&mut x1.a.clone());
        let rev_x = x.reverse();

        let result = self * &rev_x;

        let result1 = Poly::new(result.a[a - 1..a - 1 + b].to_vec(), self.ring);
        let result2 = Poly::new(result.a[c + a - 1..c + a - 1 + b].to_vec(), self.ring);

        (result1, result2)
    }

    pub fn len(&self) -> usize {
        self.a.len()
    }

    pub fn shift_mul(&mut self, k: usize) {
        self.a = vec![vec![self.ring.from(0); k], self.a.clone()].concat();
    }

    pub fn inplace_mul_monic_linear(&mut self, a: &Modulo<'a>) {
        self.a.push(self.ring.from(0));

        for i in (0..self.len() - 1).rev() {
            let p = &mut self.a[i + 1] as *mut Modulo<'a>;
            unsafe {
                *p += &self.a[i];
            }
            self.a[i] *= a;
        }
    }

    pub fn inplace_mul_t_monic_linear(&mut self, a: &Modulo<'a>) {
        for i in 0..self.len() - 1 {
            self.a[i] *= a;
            let p = &mut self.a[i] as *mut Modulo<'a>;
            unsafe {
                *p += &self.a[i + 1];
            }
        }
        self.a.pop();
    }

    pub fn mul_t_inv_monic_linear(&self, rhs: &Self) -> Modulo<'a> {
        let mut result = self.ring.from(0);
        for (i, x) in rhs.a.iter().enumerate() {
            result += &self.a[i + 1] * x;
        }
        result
    }
}

impl<'a> std::ops::Mul for &Poly<'a> {
    type Output = Poly<'a>;
    fn mul(self, rhs: Self) -> Self::Output {
        if self.len().min(rhs.len()) <= 4 {
            self.small_mul(rhs)
        } else {
            self.large_mul(rhs)
        }
    }
}

impl<'a> std::ops::SubAssign<&Poly<'a>> for Poly<'a> {
    fn sub_assign(&mut self, rhs: &Self) {
        self.set_len(self.len().max(rhs.len()), self.ring);
        for (i, x) in rhs.a.iter().enumerate() {
            self.a[i] -= x;
        }
    }
}

impl<'a> std::ops::Sub for &Poly<'a> {
    type Output = Poly<'a>;
    fn sub(self, rhs: Self) -> Self::Output {
        let mut result = self.clone();
        result -= rhs;
        result
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

pub struct MultipointEvaluation<'a> {
    pub n: usize,
    pub a: Vec<Modulo<'a>>,
    pub lower_mul_table: Vec<Poly<'a>>,
    pub upper_mul_table: Vec<Poly<'a>>,
    pub all_inv: Poly<'a>,
    ring: &'a ModuloRing,

}

impl<'a> MultipointEvaluation<'a> {
    const B1: usize = 8;
    const B2: usize = 16;

    pub fn new(point1: &Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Self {
        let n = point1.len();
        assert!(n.is_power_of_two());
        assert!(n >= Self::B1 && n >= Self::B2);

        let mut lower_mul_table = vec![Poly::zero(ring); n];
        let mut upper_mul_table = vec![Poly::zero(ring); n / Self::B2 * 2];

        for i in (0..n).step_by(Self::B2) {
            lower_mul_table[i] = Poly::one(ring);
            for j in 0..Self::B2 - 1 {
                lower_mul_table[i + j + 1] = lower_mul_table[i + j].clone();
                lower_mul_table[i + j + 1].inplace_mul_monic_linear(&(-&(point1[i + j])));
            }

            upper_mul_table[(i + n) / Self::B2] = lower_mul_table[i + Self::B2 - 1].clone();
            upper_mul_table[(i + n) / Self::B2].inplace_mul_monic_linear(&(-&point1[i + Self::B2 - 1]));
        }

        for i in (1..n / Self::B2).rev() {
            upper_mul_table[i] = &upper_mul_table[i * 2] * &upper_mul_table[i * 2 + 1];
        }
        let all_inv = Poly::inv(&upper_mul_table[1], n);

        MultipointEvaluation { n, a: point1.clone(), lower_mul_table, upper_mul_table, all_inv, ring }
    }

    fn monic_linear_prod(&self, a: &Vec<Modulo<'a>>) -> Poly<'a> {
        let mut b = vec![Poly::zero(self.ring); self.n / Self::B1];

        for i in (0..self.n).step_by(Self::B1) {
            let mut t = Poly::one(self.ring);
            for j in 0..Self::B1 {
                t.inplace_mul_monic_linear(&a[i + j]);
            }
            b[i / Self::B1] = t;
        }

        let mut k = b.len() / 2;
        while k > 0 {
            for i in 0..k {
                b[i] = &b[i * 2] * &b[i * 2 + 1];
            }
            k /= 2;
        }
        b.into_iter().nth(0).unwrap()
    }

    pub fn eval(&self, point2: &Vec<Modulo<'a>>) -> Modulo<'a> {
        let prod = self.monic_linear_prod(point2);
        let prod = &prod - &self.upper_mul_table[1];

        let mut up_tree_t = vec![Poly::zero(self.ring); self.n / Self::B2];
        let mut t = &self.all_inv * &prod;
        t.div_set_len(self.n + 1);
        up_tree_t[0] = t.reverse();

        let mut k = 1;
        while k <= self.n / Self::B2 / 2 {
            for i in (0..k).rev() {
                let tmp = up_tree_t[i].dual_mul_t(&self.upper_mul_table[(k + i) * 2], &self.upper_mul_table[(k + i) * 2 + 1]);
                (up_tree_t[i * 2 + 1], up_tree_t[i * 2]) = tmp;
            }
            k *= 2;
        }

        let mut result = self.ring.from(1);
        for i in (0..self.n).step_by(Self::B2) {
            let i2 = i / Self::B2;
            for j in (1..Self::B2).rev() {
                result *= up_tree_t[i2].mul_t_inv_monic_linear(&self.lower_mul_table[i + j]);
                up_tree_t[i2].inplace_mul_t_monic_linear(&(-&self.a[i + j]));
            }
            result *= &up_tree_t[i2].mul_t_inv_monic_linear(&self.lower_mul_table[i]);
        }
        result
    }
}
