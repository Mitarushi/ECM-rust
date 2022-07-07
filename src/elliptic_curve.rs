use ibig::modular::{ModuloRing, Modulo};

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
        let t = &p1.x * &p2.x - &p1.z * &p2.z;
        let x = &t * &t * &pm.z;
        let t = &p1.x * &p2.z - &p1.z * &p2.x;
        let z = &t * &t * &pm.x;
        EllipticPoint { x, z }
    }

    pub fn double_h(&'a self, p: &EllipticPoint<'a>, ring: &'a ModuloRing) -> EllipticPoint {
        let t = &p.x * &p.x - &p.z * &p.z;
        let x = &t * &t;
        let z = (&p.x * (&p.x + &p.z * &self.c) + &p.z * &p.z) * &p.x * &p.z * ring.from(4);
        EllipticPoint { x, z }
    }

    pub fn mul(&'a self, p: &EllipticPoint<'a>, n: u64, ring: &'a ModuloRing) -> EllipticPoint {
        if n == 0 {
            return EllipticPoint::new(ring.from(0), ring.from(0));
        } else if n == 1 {
            return p.clone();
        } else if n == 2 {
            return self.double_h(p, ring);
        }

        let mut u = p.clone();
        let mut t = self.double_h(p, ring);
        let b = 64 - n.leading_zeros();

        for j in (1..b - 1).rev() {
            if (n >> j) & 1 == 1 {
                u = self.add_h(&t, &u, p);
                t = self.double_h(&t, ring);
            } else {
                t = self.add_h(&u, &t, p);
                u = self.double_h(&u, ring);
            }
        }

        if n & 1 == 1 {
            self.add_h(&u, &t, p)
        } else {
            self.double_h(&u, ring)
        }
    }
}