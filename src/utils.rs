use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};

pub fn mod_inv<'a>(a: &Modulo<'a>, ring: &'a ModuloRing) -> Modulo<'a> {
    ring.from(mod_inv_ubig(&a.residue(), &ring.modulus()))
}

fn mod_inv_ubig(a: &UBig, n: &UBig) -> UBig {
    if a == &ubig!(1) {
        ubig!(1)
    } else if a == &ubig!(0) {
        ubig!(0)
    } else {
        n - (n * mod_inv_ubig(&(n % a), a)) / a
    }
}

pub fn gcd(a: u64, b: u64) -> u64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

pub fn bit_length(n: usize) -> u32 {
    64 - n.leading_zeros()
}