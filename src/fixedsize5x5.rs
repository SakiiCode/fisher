use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{cmp::min, ops::SubAssign, simd::Simd};

use crate::math::Quotient;

macro_rules! get {
    ($arr:ident, $r:expr, $c:expr, $cols:expr) => {
        unsafe { *$arr.get_unchecked($r * $cols + $c) }
    };
}

macro_rules! set {
    ($arr:ident, $r:expr, $c:expr, $cols:expr, $val:expr) => {
        unsafe { *$arr.get_unchecked_mut($r * $cols + $c) = $val }
    };
}

fn fill_5x5(mat_new: &mut [u32; 16], r_sum: &[u32; 5], c_sum: &[u32; 5], p_0: f64) -> f64 {
    let r = r_sum.len();
    let c = c_sum.len();

    let mut r_vec: Vec<Simd<u32, 4>> = Vec::with_capacity(r);

    for i in 0..r - 1 {
        let start = i * 4;
        let end = start + 4;
        r_vec.push(Simd::from_slice(&mat_new[start..end]));
    }
    let mut r_vec_red: Simd<u32, 4> = Simd::from_slice(&c_sum[0..4]);

    let mut c_vec: Vec<Simd<u32, 4>> = Vec::with_capacity(c);

    for i in 0..c - 1 {
        let mut arr = [0; 4];
        for j in 0..r - 1 {
            arr[j] = mat_new[j * (c - 1) + i];
        }
        c_vec.push(Simd::from_array(arr));
    }

    let mut c_vec_red: Simd<u32, 4> = Simd::from_slice(&r_sum[0..4]);

    for i in 0..c - 1 {
        c_vec_red.sub_assign(c_vec[i]);
    }

    for i in 0..r - 1 {
        r_vec_red.sub_assign(r_vec[i]);
    }

    let mut reduced = r_sum[r - 1];
    for j in r_vec_red.as_array() {
        if reduced < *j {
            //println!("");
            return 0.0;
        } else {
            reduced -= *j;
        }
    }

    let n = r_sum.iter().sum();

    let mut p_1 = Quotient::new(2 * n as usize, 2 * n as usize);

    p_1.mul_fact(r_sum);
    p_1.mul_fact(c_sum);

    p_1.div_fact(&[n; 1]);
    for i in 0..r - 1 {
        p_1.div_fact(r_vec[i].as_array());
    }
    p_1.div_fact(r_vec_red.as_array());

    /*for i in 0..c - 1 {
        p_1.div_fact(c_vec[i].as_array());
    }*/

    p_1.div_fact(c_vec_red.as_array());

    p_1.div_fact(&[reduced; 1]);

    let p_1_res = p_1.solve();
    if p_1_res <= p_0 {
        //println!(" p={}", p_1_res);
        return p_1_res;
    } else {
        //println!(" p=0.0");
        return 0.0;
    }
}

pub fn dfs_5x5(
    mat_new: &mut [u32; 16],
    xx: usize,
    yy: usize,
    r_sum: &[u32; 5],
    c_sum: &[u32; 5],
    p_0: f64,
) -> f64 {
    let r = r_sum.len();
    let c = c_sum.len();
    let mut max_1 = r_sum[xx];
    let mut max_2 = c_sum[yy];

    for j in 0..c - 1 {
        max_1 -= get!(mat_new, xx, j, c - 1);
    }

    for i in 0..r - 1 {
        max_2 -= get!(mat_new, i, yy, c - 1);
    }

    return (0..=min(max_1, max_2))
        .into_par_iter()
        .map(|k| {
            let mut mat_new2 = mat_new.clone();
            set!(mat_new2, xx, yy, c - 1, k);
            if xx + 2 == r && yy + 2 == c {
                return fill_5x5(&mut mat_new2, r_sum, c_sum, p_0);
            } else if xx + 2 == r {
                return dfs_5x5(&mut mat_new2, 0, yy + 1, r_sum, c_sum, p_0);
            } else {
                return dfs_5x5(&mut mat_new2, xx + 1, yy, r_sum, c_sum, p_0);
            }
        })
        .sum();
}
