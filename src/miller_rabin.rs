use ibig::{UBig, ubig};
use ibig::modular::{Modulo, ModuloRing};
use rand::Rng;
use rand::rngs::StdRng;

pub fn miller_rabin(n: &UBig, k: u64, thread_num: usize, rng: &mut StdRng) -> bool {
    let ring = ModuloRing::new(&n);
    let d = (n - &ubig!(1)).trailing_zeros().unwrap();
    let s = (n - &ubig!(1)) >> d;

    let main = |a: Modulo| {
        let mut x = a.pow(&s);

        if &x.residue() == &ubig!(1) || &(&x.residue() + &ubig!(1)) == n {
            return true;
        }

        for _ in 0..d - 1 {
            x = &x * &x;
            if &(&x.residue() + &ubig!(1)) == n {
                return true;
            }
        }
        false
    };

    for _ in 0..(k as usize + thread_num - 1) / thread_num {
        let a_vec = (0..thread_num).into_iter().map(|_| ring.from(rng.gen_range(ubig!(2)..n.clone()))).collect::<Vec<_>>();
        let result = a_vec.into_iter().map(|a| main(a)).collect::<Vec<_>>();

        if result.into_iter().any(|x| !x) {
            return false;
        }
    }
    true
}