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

    pub fn truncate(&mut self, n: usize) {
        self.a.truncate(n);
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

    fn large_mul(&self, rhs: &Self) -> Self {
        let n = self.len() + rhs.len() - 1;
        let padding = self.mod_log * 2 + (64 - n.leading_zeros() as usize) / 8 + 1;
        let a = self.to_rug(padding);
        let b = rhs.to_rug(padding);
        let mut result = Poly::from_rug(&(a * b), padding, self.ring);
        result.set_len(n, self.ring);
        result
    }

    fn small_mul(&self, rhs: &Self) -> Self {
        let mut result_a = vec![self.ring.from(0); self.len() + rhs.len() - 1];
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
        assert_eq!(self.a.last().unwrap().residue(), ibig::ubig!(1));

        let mut k = 1;
        let mut g = Poly::new(vec![self.ring.from(1)], self.ring);
        while k < degree {
            k *= 2;
            let mut f = self.clone();
            f.div_set_len(k);
            let mut fg = &f * &g;
            fg.div_set_len(k);
            fg = -fg;
            fg.a[k - 1] += self.ring.from(2);
            g = &fg * &g;
            g.div_set_len(k);
        }
        g.div_set_len(degree);
        g
    }

    pub fn mul_t(&self, rhs: &Self) -> Self {
        let c = self.len();
        let a = rhs.len();
        let b = c + 1 - a;

        let rev_rhs = rhs.reverse();
        let mut result = self * &rev_rhs;
        result.set_len(a + b - 1, self.ring);
        result.a = result.a[a - 1..a - 1 + b].to_vec();

        result
    }

    pub fn len(&self) -> usize {
        self.a.len()
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

impl<'a> std::ops::Sub for &Poly<'a> {
    type Output = Poly<'a>;
    fn sub(self, rhs: Self) -> Self::Output {
        let n = self.len().max(rhs.len());
        let mut c = vec![self.ring.from(0); n];
        for i in 0..self.len() {
            c[i] += &self.a[i];
        }
        for i in 0..rhs.len() {
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

pub struct MultipointEvaluation<'a> {
    pub n: usize,
    pub mul_table: Vec<Poly<'a>>,
    pub all_inv: Poly<'a>,
    ring: &'a ModuloRing,
}

impl<'a> MultipointEvaluation<'a> {
    pub fn new(point1: &Vec<Modulo<'a>>, ring: &'a ModuloRing) -> Self {
        let n = point1.len();
        assert!(n.is_power_of_two());

        let mut mul_table = vec![Poly::zero(ring); n * 2];

        for i in 0..n {
            mul_table[i + n] = Poly::monic_linear(-point1[i].clone(), ring);
        }

        for i in (1..n).rev() {
            mul_table[i] = &mul_table[i * 2] * &mul_table[i * 2 + 1];
        }

        let all_inv = Poly::inv(&mul_table[1], n + 1);

        MultipointEvaluation { n, mul_table, all_inv, ring }
    }

    fn prod(&self, a: Vec<Poly<'a>>) -> Poly<'a> {
        let mut a = a;
        let mut k = self.n / 2;
        while k > 0 {
            for i in 0..k {
                a[i] = &a[i * 2] * &a[i * 2 + 1];
            }
            k /= 2;
        }
        a.into_iter().nth(0).unwrap()
    }

    pub fn eval(&self, point2: &Vec<Modulo<'a>>) -> Modulo<'a> {
        let prod = self.prod(point2.iter().map(|x| Poly::monic_linear(x.clone(), self.ring)).collect());
        let prod = &prod - &self.mul_table[1];

        let mut up_tree_t = vec![Poly::zero(self.ring); self.n * 2];
        let mut t = &self.all_inv * &prod;
        t.div_set_len(self.n + 1);
        up_tree_t[1] = t.reverse();

        for i in 1..self.n {
            up_tree_t[i * 2] = up_tree_t[i].mul_t(&self.mul_table[i * 2 + 1]);
            up_tree_t[i * 2 + 1] = up_tree_t[i].mul_t(&self.mul_table[i * 2]);
        }

        let mut result = self.ring.from(1);
        for i in 0..self.n {
            result *= &up_tree_t[i + self.n].a[1];
        }
        result
    }
}
