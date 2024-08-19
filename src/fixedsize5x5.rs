use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    cmp::min,
    ops::SubAssign,
    simd::{num::SimdUint, LaneCount, Simd, SupportedLaneCount},
};

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

fn fill<const N: usize>(
    mat_new: &mut [u32; N * N],
    r_sum: &[u32; N + 1],
    c_sum: &[u32; N + 1],
    p_0: f64,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
    let mut r_vec: Vec<Simd<u32, N>> = Vec::with_capacity(N);

    for i in 0..N {
        let start = i * N;
        r_vec.push(Simd::from_slice(&mat_new[start..]));
    }
    let mut r_vec_red: Simd<u32, N> = Simd::from_slice(c_sum);

    for i in 0..N {
        r_vec_red.sub_assign(r_vec[i]);
    }

    let r_red_sum = r_vec_red.reduce_sum();
    let mut r_last = r_sum[N];

    if r_last < r_red_sum {
        //println!("");
        return 0.0;
    }
    r_last -= r_red_sum;

    let mut c_vec: Vec<Simd<u32, N>> = Vec::with_capacity(N);

    for i in 0..N {
        let mut arr = [0; N];
        for j in 0..N {
            arr[j] = mat_new[j * N + i];
        }
        c_vec.push(Simd::from_array(arr));
    }

    let mut c_vec_red: Simd<u32, N> = Simd::from_slice(r_sum);

    for i in 0..N {
        c_vec_red.sub_assign(c_vec[i]);
    }

    let n = r_sum.iter().sum();

    let mut p_1 = Quotient::new(2 * n as usize, 2 * n as usize);

    p_1.mul_fact(r_sum);
    p_1.mul_fact(c_sum);

    p_1.div_fact(&[n; 1]);
    for i in 0..N {
        p_1.div_fact(r_vec[i].as_array());
    }
    p_1.div_fact(r_vec_red.as_array());

    p_1.div_fact(c_vec_red.as_array());

    p_1.div_fact(&[r_last; 1]);

    let p_1_res = p_1.solve();
    if p_1_res <= p_0 {
        //println!(" p={}", p_1_res);
        p_1_res
    } else {
        //println!(" p=0.0");
        0.0
    }
}

pub fn dfs<const N: usize>(
    mat_new: &mut [u32; N * N],
    xx: usize,
    yy: usize,
    r_sum: &[u32; N + 1],
    c_sum: &[u32; N + 1],
    p_0: f64,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
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
                return fill::<N>(&mut mat_new2, r_sum, c_sum, p_0);
            } else if xx + 2 == r {
                return dfs::<N>(&mut mat_new2, 0, yy + 1, r_sum, c_sum, p_0);
            } else {
                return dfs::<N>(&mut mat_new2, xx + 1, yy, r_sum, c_sum, p_0);
            }
        })
        .sum();
}
