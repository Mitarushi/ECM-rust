use ibig::modular::{Modulo, ModuloRing};

use crate::{bit_length, mod_inv};

#[derive(Debug, Clone)]
pub struct EllipticPoint<'a> {
    pub x: Modulo<'a>,
    pub z: Modulo<'a>,
}

#[derive(Debug, Clone)]
pub struct EllipticCurve<'a> {
    pub c: Modulo<'a>,
}

impl<'a> EllipticPoint<'a> {
    pub fn new(x: Modulo<'a>, z: Modulo<'a>) -> Self {
        EllipticPoint { x, z }
    }

    pub fn affine_x(&self, ring: &'a ModuloRing) -> Modulo<'a> {
        &self.x * mod_inv(&self.z, ring)
    }

    pub fn zero(ring: &'a ModuloRing) -> Self {
        EllipticPoint::new(ring.from(1), ring.from(0))
    }
}

enum ChainAddInfo {
    AddToLeft,
    AddToRight,
}


fn generate_chain(n: u64, hint: u64) -> Vec<ChainAddInfo> {
    let mut result = Vec::new();
    let mut p = n;
    let mut q = hint;

    while p != 1 || q != 1 {
        if p > q {
            result.push(ChainAddInfo::AddToLeft);
            p -= q;
        } else {
            result.push(ChainAddInfo::AddToRight);
            q -= p;
        }
    }

    result.reverse();
    result
}

impl<'a> EllipticCurve<'a> {
    pub fn new(c: Modulo<'a>, ring: &'a ModuloRing) -> Self {
        let c = (c + ring.from(2)) * mod_inv(&ring.from(4), ring);
        EllipticCurve { c }
    }

    pub fn add_h(&self, p1: &EllipticPoint<'a>, p2: &EllipticPoint<'a>, pm: &EllipticPoint<'a>) -> EllipticPoint {
        let t1 = (&p1.x - &p1.z) * (&p2.x + &p2.z);
        let t2 = (&p1.x + &p1.z) * (&p2.x - &p2.z);
        let t = &t1 + &t2;
        let x = &t * &t * &pm.z;
        let t = &t1 - &t2;
        let z = &t * &t * &pm.x;
        EllipticPoint { x, z }
    }

    pub fn double_h(&'a self, p: &EllipticPoint<'a>) -> EllipticPoint {
        let u = &p.x + &p.z;
        let u = &u * &u;
        let v = &p.x - &p.z;
        let v = &v * &v;
        let k = &u - &v;
        let t = &self.c * &k + &v;
        let x = &u * &v;
        let z = &k * &t;
        EllipticPoint { x, z }
    }

    pub fn triple_h(&'a self, p: &EllipticPoint<'a>) -> EllipticPoint {
        let u1 = &p.x + &p.z;
        let u = &u1 * &u1;
        let v1 = &p.x - &p.z;
        let v = &v1 * &v1;
        let k = &u - &v;
        let x = &u * &v;
        let z = &k * (&self.c * &k + &v);

        let t1 = (&x - &z) * &u1;
        let t2 = (&x + &z) * &v1;
        let t = &t1 + &t2;
        let x = &t * &t * &p.z;
        let t = &t1 - &t2;
        let z = &t * &t * &p.x;

        EllipticPoint { x, z }
    }

    pub fn _mul(&'a self, p: &EllipticPoint<'a>, n: u64, ring: &'a ModuloRing) -> EllipticPoint {
        if n == 0 {
            return EllipticPoint::zero(ring);
        } else if n == 1 {
            return p.clone();
        } else if n == 2 {
            return self.double_h(p);
        }

        let mut u = p.clone();
        let mut t = self.double_h(p);
        let b = bit_length(n as usize);

        for j in (1..b - 1).rev() {
            if (n >> j) & 1 == 1 {
                u = self.add_h(&t, &u, p);
                t = self.double_h(&t);
            } else {
                t = self.add_h(&u, &t, p);
                u = self.double_h(&u);
            }
        }

        if n & 1 == 1 {
            self.add_h(&u, &t, p)
        } else {
            self.double_h(&u)
        }
    }

    pub fn mul_with_hint(&'a self, p: &EllipticPoint<'a>, n: u64, hint: u64, ring: &'a ModuloRing) -> EllipticPoint {
        if n == 0 {
            return EllipticPoint::zero(ring);
        } else if n == 1 {
            return p.clone();
        }

        let chain = generate_chain(n, hint);
        let mut u = self.double_h(p);
        let mut v = p.clone();
        let mut t = p.clone();

        match chain[0] {
            ChainAddInfo::AddToRight => (u, v) = (v, u),
            _ => (),
        }

        for c in chain.into_iter().skip(1) {
            match c {
                ChainAddInfo::AddToLeft => {
                    let tmp = self.add_h(&u, &v, &t);
                    t = u;
                    u = tmp;
                }
                ChainAddInfo::AddToRight => {
                    let tmp = self.add_h(&v, &u, &t);
                    t = v;
                    v = tmp;
                }
            }
        }
        u
    }
}