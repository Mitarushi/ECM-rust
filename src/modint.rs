use num_bigint::{BigInt, ToBigInt};
use num_traits::{One, Zero};
use num_integer::Integer;


#[derive(Debug, Clone)]
pub struct Montgomery
{
    pub n: BigInt,
    pub r_log: u64,
    pub n_inv: BigInt,
    pub r_mask: BigInt,
    pub r2: BigInt,
}


impl Montgomery
{
    pub fn new(n: BigInt) -> Self {
        let r_log = n.bits();
        let r = BigInt::one() << r_log;
        let n_inv = &r - Montgomery::get_n_inv(r_log, &n);
        let r_mask = &r - BigInt::one();
        let r2 = &r * &r % &n;
        Montgomery { n, r_log, n_inv, r_mask, r2 }
    }

    fn get_n_inv(r_log: u64, n: &BigInt) -> BigInt {
        let mut precision = 1;
        let mut res = n.clone();
        let mask = (BigInt::one() << r_log) - BigInt::one();
        while precision < r_log {
            res = &res * ((2.to_bigint().unwrap() - n * &res) & &mask);
            res = &res & &mask;
            precision *= 2;
        }
        res
    }

    pub fn mr(&self, t: &BigInt) -> BigInt {
        let tn = (t & &self.r_mask) * &self.n_inv;
        let tn = tn & &self.r_mask;
        let tn = t + tn * &self.n;
        let tn = tn >> self.r_log;
        if &tn >= &self.n {
            tn - &self.n
        } else {
            tn
        }
    }

    pub fn add(&self, x: &ModInt, y: &ModInt) -> ModInt {
        let mut z = &x.x + &y.x;
        if &z >= &self.n {
            z = z - &self.n;
        }
        ModInt { x: z }
    }

    pub fn sub(&self, x: &ModInt, y: &ModInt) -> ModInt {
        let mut z = &x.x - &y.x;
        if &z < &BigInt::zero() {
            z = z + &self.n;
        }
        ModInt { x: z }
    }

    pub fn mul(&self, x: &ModInt, y: &ModInt) -> ModInt {
        ModInt { x: self.mr(&(&x.x * &y.x)) }
    }

    pub fn square(&self, x: &ModInt) -> ModInt {
        self.mul(x, x)
    }

    pub fn pow(&self, x: &ModInt, n: u64) -> ModInt {
        let mut res = x.clone();
        for _ in 0..n - 1 {
            res = self.mul(&res, x);
        }
        res
    }

    pub fn binary_pow(&self, x: &ModInt, n: &BigInt) -> ModInt {
        let mut res = self.convert_u64(1);
        let mut n = n.clone();
        let mut tmp = x.clone();
        while !n.is_zero() {
            if n.is_odd() {
                res = self.mul(&res, &tmp);
            }
            tmp = self.mul(&tmp, &tmp);
            n >>= 1;
        }
        res
    }

    pub fn convert_bigint(&self, x: &BigInt) -> ModInt {
        ModInt { x: self.mr(&(x * &self.r2)) }
    }

    pub fn convert_u64(&self, x: u64) -> ModInt {
        self.convert_bigint(&x.into())
    }

    pub fn val(&self, x: &ModInt) -> BigInt {
        self.mr(&x.x)
    }
}

#[derive(Debug, Clone)]
pub struct ModInt {
    pub x: BigInt,
}

impl ModInt {
    pub fn new(x: BigInt) -> Self {
        ModInt { x }
    }

    pub fn zero() -> Self {
        ModInt { x: BigInt::zero() }
    }
}
