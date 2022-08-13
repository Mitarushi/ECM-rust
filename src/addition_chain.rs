pub fn continued_fraction_len(x: u64, y: u64) -> u32 {
    let mut x = x;
    let mut y = y;

    if x < y {
        (x, y) = (y, x);
    }

    let mut result = 0;

    while x.max(y) > 1 {
        (x, y) = (y, x - y);
        result += 1;

        if x < y {
            if x == 0 {
                return u32::MAX / 2;
            }

            let d = (y - 1) / x;
            y = y - d * x;
            result += d as u32;
        }
    }

    result
}

const BLOCK_SIZE: usize = 16;

fn continued_fraction_len_block(n: u32, a: [u32; BLOCK_SIZE], max_len: u32) -> u32 {
    let mut x = [n; BLOCK_SIZE];
    let mut y = a;

    for i in 0..max_len - 1 {
        for j in 0..BLOCK_SIZE {
            x[j] -= y[j];
            if x[j] < y[j] {
                (x[j], y[j]) = (y[j], x[j]);
            }
        }

        let is_one = x.iter().find(|&&x| x == 1).is_some();
        if is_one {
            return i + 1;
        }
    }

    max_len
}

const SEARCH_WIDTH: isize = 500;

pub fn compute_optimal_hint(n: u64) -> u64 {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let x = (n as f64 / phi).round() as u64;

    let mut min_len = continued_fraction_len(n, x) + 1;
    let mut min_block = 0;

    for i in -SEARCH_WIDTH..SEARCH_WIDTH {
        let lo = (x as isize + i * BLOCK_SIZE as isize).max(1) as u32;
        let mut a = [0; BLOCK_SIZE];
        for i in 0..BLOCK_SIZE {
            a[i] = (lo + i as u32).min(n as u32 - 1);
        }
        let len = continued_fraction_len_block(n as u32, a, min_len);

        if len < min_len {
            min_len = len;
            min_block = i;
        }
    }

    let lo = (x as isize + min_block * BLOCK_SIZE as isize).max(1) as u32;
    let mut optimal_hint = lo;
    for i in 0..BLOCK_SIZE {
        let a = (lo + i as u32).min(n as u32 - 1);

        if continued_fraction_len(n, a as u64) <= min_len {
            optimal_hint = a;
        }
    }

    optimal_hint as u64
}