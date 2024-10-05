use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    cell::RefCell,
    cmp::min,
    convert::Infallible,
    ops::SubAssign,
    simd::{num::SimdInt, LaneCount, Simd, SupportedLaneCount},
};
use thread_local::ThreadLocal;

use crate::math::Quotient;

macro_rules! set {
    ($arr:ident, $r:expr, $c:expr, $cols:expr, $val:expr) => {
        unsafe { *$arr.get_unchecked_mut($r * $cols + $c) = $val }
    };
}

// faster if not inlined
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

    let col = r_vec.map(|row| row.reduce_sum());
    let col_simd = Simd::from(col);
    c_vec_red.sub_assign(col_simd);

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

pub fn calculate(table: Vec<Vec<i32>>) -> Result<f64, Infallible> {
    let mut row_sum: Vec<i32> = table.iter().map(|row| row.iter().sum()).collect();
    let mut col_sum: Vec<i32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    if table.len() != table[0].len() {
        println!("fisher.calculate() can only be used with square matrices yet!");
        return Ok(-1.0);
    }

    let n: i32 = row_sum.iter().sum();

    let mut p_0 = Quotient::new(n.try_into().unwrap(), &[], &[]);

    p_0.mul_fact(&row_sum);
    p_0.mul_fact(&col_sum);

    p_0.div_fact(&[n; 1]);
    p_0.div_fact(&table.iter().flatten().cloned().collect::<Vec<i32>>());

    let stat = p_0.solve() + f64::EPSILON;

    let lanes: usize = match table.len() {
        1 => {
            println!("Invalid table size!");
            return Ok(-1.0);
        }
        2 => 1,
        3 => 2,
        4 | 5 => 4,
        6..=9 => 8,
        10..=17 => 16,
        _ => {
            println!("fisher.calculate() can only be used with up to 16x16 matrices!");
            return Ok(-1.0);
        }
    };
    let tl = ThreadLocal::new();

    row_sum.resize(lanes + 1, 0);
    col_sum.resize(lanes + 1, 0);

    let p = match lanes {
        1 => dfs::<1>(
            &mut [0; 1],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
        ),
        2 => dfs::<2>(
            &mut [0; 4],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
        ),
        4 => dfs::<4>(
            &mut [0; 16],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
        ),
        8 => dfs::<8>(
            &mut [0; 64],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
        ),
        16 => dfs::<16>(
            &mut [0; 256],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
        ),
        _ => {
            println!("Error in matching lane count. This should never happen");
            return Ok(-1.0);
        }
    };

    return Ok(p);
}

#[test]
fn fixed2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.5920745920745918,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed3x3() {
    let input = vec![vec![32, 10, 20], vec![20, 25, 18], vec![11, 17, 14]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.011074529608901276,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.010961244321907074,
        epsilon = 0.000001
    ));
}

#[test]
#[ignore]
fn fixed5x5_large() {
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
        vec![1, 1, 0, 0, 3],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.24678711559405725,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed5x5_small() {
    let input = vec![
        vec![1, 0, 0, 0, 0],
        vec![1, 1, 0, 1, 0],
        vec![1, 1, 0, 0, 1],
        vec![0, 0, 1, 2, 1],
        vec![1, 1, 2, 1, 1],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.9712771262351092);
}
