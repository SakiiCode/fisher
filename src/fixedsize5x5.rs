use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    cell::RefCell,
    cmp::min,
    ops::{IndexMut, SubAssign},
    simd::{num::SimdInt, LaneCount, Simd, SupportedLaneCount},
};
use thread_local::ThreadLocal;

use crate::math::Quotient;

/*macro_rules! get {
    ($arr:ident, $r:expr, $c:expr, $cols:expr) => {
        unsafe { *$arr.get_unchecked($r * $cols + $c) }
    };
}*/

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
    /*
    let mut r_vec: Vec<Simd<i32, N>> = Vec::with_capacity(N);

    for i in 0..N {
        let start = i * N;
        r_vec.push(Simd::from_slice(&mat_new[start..]));
    }
    */

    let mut r_vec_red: Simd<i32, N> = Simd::from_slice(c_sum);

    /*
    let r_vec: Box<[Simd<i32, N>]> = (0..N)
    .map(|i| {
        let start = i * N;
        let row_simd = Simd::from_slice(&mat_new[start..]);
        r_vec_red.sub_assign(row_simd);
        row_simd
    })
    .collect();
    */

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

    /*
    let mut c_vec: Vec<Simd<i32, N>> = Vec::with_capacity(N);

    for i in 0..N {
        let mut arr = [0; N];
        for j in 0..N {
            arr[j] = mat_new[j * N + i];
        }
        c_vec.push(Simd::from_array(arr));
    }
    */

    let mut c_vec_red: Simd<i32, N> = Simd::from_slice(r_sum);
    /*
    let mut c_vec = r_vec.clone();
    let max_step_size = usize::ilog2(N / 2);
    let step_size_iter = (1..=max_step_size).map(|x| 2usize.pow(x)).rev();
    for step_size in step_size_iter {
        for index in 0..step_size {
            let (v1, v2) = c_vec[index].interleave(c_vec[index + step_size]);
            c_vec[index] = v1;
            c_vec[index + step_size] = v2;
        }
    }

    for i in (0..N).step_by(2) {
        let (v1, v2) = c_vec[i].interleave(c_vec[i + 1]);
        c_vec[i] = v1;
        c_vec[i + 1] = v2;
    }

    for i in 0..N {
        c_vec_red.sub_assign(c_vec[i]);
    }
    */

    for i in 0..N {
        /*let mut col = [0; N];
        for j in 0..N {
            col[j] = mat_new[j * N + i];
        }*/

        let mut col_simd = Simd::from_array([0; N]);
        for j in 0..N {
            *col_simd.index_mut(j) = mat_new[j * N + i];
        }
        c_vec_red.sub_assign(col_simd);
    }

    // r_sum is N+1 length, SIMD cannot be used
    let n: i32 = r_sum.iter().sum();

    //let mut p_1 = Quotient::new(2 * n as usize, 2 * n as usize);
    let p_1_ref = tl.get_or(|| {
        let mut init_n = Vec::with_capacity(2 * (N + 1));
        let init_d = vec![n];
        init_n.extend_from_slice(r_sum);
        init_n.extend_from_slice(c_sum);
        Box::new(RefCell::new(Quotient::new(n as usize, &init_n, &init_d)))
    });

    let mut p_1 = (p_1_ref).borrow_mut();
    p_1.clear();

    //p_1.mul_fact(r_sum);
    //p_1.mul_fact(c_sum);

    //p_1.div_fact(&[n; 1]);
    //foreach is slower here
    /*for i in 0..N {
        p_1.div_fact(r_vec[i].as_array());
    }*/
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
    //let max_1 = r_sum[xx] - &mat_new[xx * N..(xx + 1) * N].iter().sum();

    for index in (yy..N * N).step_by(N) {
        max_2 -= mat_new[index]; //get!(mat_new, i, yy, c - 1);
    }

    /*
    // last row and col is empty, loops are needed
    for j in 0..c - 1 {
        max_1 -= get!(mat_new, xx, j, c - 1);
    }

    for i in 0..r - 1 {
        max_2 -= get!(mat_new, i, yy, c - 1);
    }
    */

    /*let max_1 = r_sum[xx] - Simd::from_slice(&mat_new[xx * N..]).reduce_sum();

    let mut max_2 = c_sum[yy];
    for i in 0..r - 1 {
        max_2 -= get!(mat_new, i, yy, c - 1);
    }*/

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

    if xx + 2 != r {
        return (0..=min(max_1, max_2))
            .into_par_iter()
            .map(next_cycle)
            .sum();
    } else {
        return (0..=min(max_1, max_2)).map(next_cycle).sum();
    }
}
