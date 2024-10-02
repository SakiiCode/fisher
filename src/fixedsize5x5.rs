use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    cell::RefCell,
    cmp::min,
    ops::{IndexMut, SubAssign},
    simd::{num::SimdInt, LaneCount, Simd, SupportedLaneCount},
};
use thread_local::ThreadLocal;

use crate::math::Quotient;

macro_rules! set {
    ($arr:ident, $r:expr, $c:expr, $cols:expr, $val:expr) => {
        unsafe { *$arr.get_unchecked_mut($r * $cols + $c) = $val }
    };
}

#[inline(never)]
fn fill<const N: usize>(
    mat_new: &mut [i32; N * N],
    r_sum: &[i32; N + 1],
    c_sum: &[i32; N + 1],
    p_0: f64,
    tl: &ThreadLocal<Box<RefCell<Quotient>>>,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
    let mut r_vec_red: Simd<i32, N> = Simd::from_slice(c_sum);

    let mut r_vec = [Simd::from_array([0; N]); N];
    for i in 0..N {
        let start = i * N;
        r_vec[i] = Simd::from_slice(&mat_new[start..]);
        r_vec_red.sub_assign(r_vec[i]);
    }

    let r_red_sum = r_vec_red.reduce_sum();
    let mut r_last = r_sum[N];

    if r_last < r_red_sum {
        //println!("");
        return 0.0;
    }
    r_last -= r_red_sum;

    let mut c_vec_red: Simd<i32, N> = Simd::from_slice(r_sum);

    for i in 0..N {
        let mut col_simd = Simd::from_array([0; N]);
        for j in 0..N {
            *col_simd.index_mut(j) = mat_new[j * N + i];
        }
        c_vec_red.sub_assign(col_simd);
    }

    // r_sum is N+1 length, SIMD cannot be used
    let n: i32 = r_sum.iter().sum();

    let p_1_ref = tl.get_or(|| {
        let mut init_n = Vec::with_capacity(2 * (N + 1));
        let init_d = vec![n];
        init_n.extend_from_slice(r_sum);
        init_n.extend_from_slice(c_sum);
        Box::new(RefCell::new(Quotient::new(n as usize, &init_n, &init_d)))
    });

    let mut p_1 = (p_1_ref).borrow_mut();
    p_1.clear();

    p_1.div_fact(mat_new);
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
    mat_new: &mut [i32; N * N],
    xx: usize,
    yy: usize,
    r_sum: &[i32; N + 1],
    c_sum: &[i32; N + 1],
    p_0: f64,
    tl: &ThreadLocal<Box<RefCell<Quotient>>>,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
    let r = r_sum.len();
    let c = c_sum.len();
    // do not put this below max_1
    let mut max_2 = c_sum[yy];

    let max_1 = r_sum[xx] - Simd::from_slice(&mat_new[xx * N..]).reduce_sum();

    for index in (yy..N * N).step_by(N) {
        max_2 -= mat_new[index]; //get!(mat_new, i, yy, c - 1);
    }

    let next_cycle = |k| {
        let mut mat_new2 = mat_new.clone();
        set!(mat_new2, xx, yy, c - 1, k);
        if xx + 2 == r && yy + 2 == c {
            return fill::<N>(&mut mat_new2, r_sum, c_sum, p_0, tl);
        } else if xx + 2 == r {
            return dfs::<N>(&mut mat_new2, 0, yy + 1, r_sum, c_sum, p_0, tl);
        } else {
            return dfs::<N>(&mut mat_new2, xx + 1, yy, r_sum, c_sum, p_0, tl);
        }
    };

    if yy == 0 {
        return (0..=min(max_1, max_2))
            .into_par_iter()
            .map(next_cycle)
            .sum();
    } else {
        return (0..=min(max_1, max_2)).map(next_cycle).sum();
    }
}
