use ibig::modular::{Modulo, ModuloRing};

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
}

impl<'a> EllipticCurve<'a> {
    pub fn new(c: Modulo<'a>) -> Self {
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
        let x2 = &p.x * &p.x;
        let z2 = &p.z * &p.z;
        let xz = &p.x * &p.z;
        let t = &x2 - &z2;
        let x = &t * &t;
        let t = (&x2 + &xz * &self.c + &z2) * &xz;
        let t = &t + &t;
        let z = &t + &t;
        EllipticPoint { x, z }
    }

    pub fn mul(&'a self, p: &EllipticPoint<'a>, n: u64, ring: &'a ModuloRing) -> EllipticPoint {
        if n == 0 {
            return EllipticPoint::new(ring.from(0), ring.from(0));
        } else if n == 1 {
            return p.clone();
        } else if n == 2 {
            return self.double_h(p);
        }

        let mut u = p.clone();
        let mut t = self.double_h(p);
        let b = 64 - n.leading_zeros();

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
}