use crate::modint::{ModInt, Montgomery};

#[derive(Debug, Clone)]
pub struct EllipticPoint  {
    pub x: ModInt,
    pub z: ModInt,
}

#[derive(Debug, Clone)]
pub struct EllipticCurve  {
    pub c: ModInt,
}

impl EllipticPoint {
    pub fn new(x: ModInt, z: ModInt) -> Self {
        EllipticPoint { x, z }
    }

    pub fn zero() -> Self {
        EllipticPoint { x: ModInt::zero(), z: ModInt::zero() }
    }
}

impl EllipticCurve {
    pub fn new(c: ModInt) -> Self {
        EllipticCurve { c }
    }

    pub fn add_h(&self, p1: &EllipticPoint, p2: &EllipticPoint, pm: &EllipticPoint,
                 mr: &Montgomery) -> EllipticPoint {
        let t = mr.sub(&mr.mul(&p1.x, &p2.x), &mr.mul(&p1.z, &p2.z));
        let x = mr.mul(&mr.square(&t), &pm.z);
        let t = mr.sub(&mr.mul(&p1.x, &p2.z), &mr.mul(&p1.z, &p2.x));
        let z = mr.mul(&mr.square(&t), &pm.x);
        EllipticPoint { x, z }
    }

    pub fn double_h(&self, p: &EllipticPoint, mr: &Montgomery) -> EllipticPoint {
        let x = mr.square(&mr.sub(&mr.square(&p.x), &mr.square(&p.z)));
        let t = mr.add(&p.x, &mr.mul(&p.z, &self.c));
        let t = mr.add(&mr.mul(&p.x, &t), &mr.square(&p.z));
        let t = mr.mul(&mr.mul(&t, &p.x), &p.z);
        let z = mr.mul(&t, &mr.convert_u64(4));
        EllipticPoint { x, z }
    }

    pub fn mul(&self, p: &EllipticPoint, n: u64, mr: &Montgomery) -> EllipticPoint {
        if n == 0 {
            return EllipticPoint::zero();
        } else if n == 1 {
            return p.clone();
        } else if n == 2 {
            return self.double_h(p, mr);
        }

        let mut u = p.clone();
        let mut t = self.double_h(p, mr);
        let b = 64 - n.leading_zeros();

        for j in (1 .. b-1).rev() {
            if (n >> j) & 1 == 1 {
                u = self.add_h(&t, &u, p, mr);
                t = self.double_h(&t, mr);
            } else {
                t = self.add_h(&u, &t, p, mr);
                u = self.double_h(&u, mr);
            }
        }

        if n & 1 == 1 {
            self.add_h(&u, &t, p, mr)
        } else {
            self.double_h(&u, mr)
        }
    }
}