use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{
    cell::RefCell,
    cmp::min,
    convert::Infallible,
    num::NonZero,
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
    mat_new: [i32; N * N],
    r_sum: &[i32; N + 1],
    c_sum: &[i32; N + 1],
    p_0: f64,
    tl: &ThreadLocal<Box<RefCell<Quotient>>>,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
    //print!("{:?} -> ", &mat_new);
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
        //println!();
        return 0.0;
    }
    r_last -= r_red_sum;

    let mut c_vec_red: Simd<i32, N> = Simd::from_slice(r_sum);

    let col = r_vec.map(|row| row.reduce_sum());
    let col_simd = Simd::from(col);
    c_vec_red.sub_assign(col_simd);

    let p_1_ref = tl.get_or(|| {
        // r_sum is N+1 length, SIMD cannot be used
        let n: i32 = r_sum.iter().sum();
        let init_d = vec![n];
        let mut init_n = Vec::with_capacity(2 * (N + 1) + 1);
        init_n.extend_from_slice(r_sum);
        init_n.extend_from_slice(c_sum);
        Box::new(RefCell::new(Quotient::new(n as usize, &init_n, &init_d)))
    });

    let mut p_1 = (p_1_ref).borrow_mut();
    p_1.clear();

    p_1.div_fact(&mat_new);
    p_1.div_fact(r_vec_red.as_array());

    p_1.div_fact(c_vec_red.as_array());

    p_1.div_fact(&[r_last]);

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
    mat_new: [i32; N * N],
    xx: usize,
    yy: usize,
    r_sum: &[i32; N + 1],
    c_sum: &[i32; N + 1],
    p_0: f64,
    tl: &ThreadLocal<Box<RefCell<Quotient>>>,
    threads: i32,
    max_threads: i32,
) -> f64
where
    LaneCount<N>: SupportedLaneCount,
{
    let r = r_sum.len();
    let max_1 = r_sum[xx] - Simd::from_slice(&mat_new[xx * N..]).reduce_sum();

    let c = c_sum.len();
    let mut max_2 = c_sum[yy];

    for index in (yy..N * N).step_by(N) {
        max_2 -= mat_new[index]; //get!(mat_new, i, yy, c - 1);
    }

    let max = min(max_1, max_2);

    let next_cycle = |k| {
        let mut mat_new2 = mat_new.clone();
        set!(mat_new2, xx, yy, c - 1, k);
        if xx + 2 == r && yy + 2 == c {
            return fill::<N>(mat_new2, r_sum, c_sum, p_0, tl);
        }

        let (next_x, next_y) = next_cell(xx, yy, r, c);
        return dfs::<N>(
            mat_new2,
            next_x,
            next_y,
            r_sum,
            c_sum,
            p_0,
            tl,
            threads + max,
            max_threads,
        );
    };

    if threads < max_threads {
        return (0..=max).into_par_iter().map(next_cycle).sum();
    } else {
        return (0..=max).map(next_cycle).sum();
    }
}

fn next_cell(xx: usize, yy: usize, r: usize, c: usize) -> (usize, usize) {
    let mut next_x = xx;
    let mut next_y = yy;
    let yellow = ((xx + yy) % 2) == 0;
    if (yellow && xx == 0) || (!yellow && xx + 2 == c) {
        //print!("red ");
        next_y += 1;
    } else if (!yellow && yy == 0) || (yellow && yy + 2 == r) {
        //print!("blue ");
        next_x += 1;
    } else if yellow {
        //print!("yellow ");
        next_y += 1;
        next_x -= 1;
    } else {
        //print!("green ");
        next_x += 1;
        next_y -= 1
    }
    //println!("{},{}", next_x, next_y);
    return (next_x, next_y);
}

pub fn calculate(table: Vec<Vec<i32>>) -> Result<f64, Infallible> {
    let mut row_sum: Vec<i32> = table.iter().map(|row| row.iter().sum()).collect();
    let mut col_sum: Vec<i32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    let n: i32 = row_sum.iter().sum();
    let seq = &table.iter().flatten().cloned().collect::<Vec<i32>>();

    if seq.iter().any(|x| *x < 0) {
        println!("ERROR: Negative value in the table!");
        return Ok(-1.0);
    }

    if seq.iter().all(|x| *x == 0) {
        println!("ERROR: All elements in the table are zero!");
        return Ok(-1.0);
    }

    let mut p_0 = Quotient::new(n.try_into().unwrap(), &[], &[]);

    p_0.mul_fact(&row_sum);
    p_0.mul_fact(&col_sum);

    p_0.div_fact(&[n; 1]);
    p_0.div_fact(seq);

    let stat = p_0.solve() + f64::EPSILON;
    let needed_lanes = usize::max(row_sum.len(), col_sum.len());

    let lanes = match needed_lanes {
        1 => {
            println!("ERROR: Invalid table size!");
            return Ok(-1.0);
        }
        2 => 1,
        3 => 2,
        4 | 5 => 4,
        6..=9 => 8,
        _ => {
            println!("ERROR: fisher.recursive() can only be used with up to 9x9 matrices!");
            return Ok(-1.0);
        }
    };
    let tl = ThreadLocal::new();

    row_sum.resize(lanes + 1, 0);
    col_sum.resize(lanes + 1, 0);

    row_sum.sort();
    col_sum.sort();

    let max_threads = std::thread::available_parallelism()
        .unwrap_or(NonZero::new(12).unwrap())
        .get() as i32;

    let p = match lanes {
        1 => dfs::<1>(
            [0; 1],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
            0,
            max_threads,
        ),
        2 => dfs::<2>(
            [0; 4],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
            0,
            max_threads,
        ),
        4 => dfs::<4>(
            [0; 16],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
            0,
            max_threads,
        ),
        8 => dfs::<8>(
            [0; 64],
            0,
            0,
            &row_sum.try_into().unwrap(),
            &col_sum.try_into().unwrap(),
            stat,
            &tl,
            0,
            max_threads,
        ),
        _ => {
            println!("Error in matching lane count. This should never happen");
            return Ok(-1.0);
        }
    };

    return Ok(p);
}

#[test]
fn fixed1x1_error() {
    let input = vec![vec![5]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert_eq!(output, -1.0);
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
fn fixed2x2_error() {
    let input = vec![vec![3, 4], vec![4, -2]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert_eq!(output, -1.0);
}

#[test]
fn fixed3x2() {
    let input = vec![vec![1000, 626, 782], vec![976, 814, 892]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, 0.0, epsilon = 0.000001));
}

#[test]
fn fixed3x3() {
    let input = vec![vec![32, 10, 20], vec![20, 25, 18], vec![11, 17, 14]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.010967949934049852,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed3x3_unit() {
    let input = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert_eq!(output, 1.0);
}

#[test]
fn fixed3x3_zero() {
    let input = vec![vec![0, 0, 0], vec![0, 0, 0], vec![0, 0, 0]];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, -1.0, epsilon = 0.000001));
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
fn fixed4x4_large() {
    let input = vec![
        vec![28, 28, 28, 0],
        vec![0, 0, 0, 16],
        vec![0, 0, 0, 5],
        vec![0, 0, 0, 7],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, 0.0, epsilon = 0.000001));
}

#[test]
#[ignore]
fn fixed4x5_large() {
    let input = vec![
        vec![8, 3, 5, 5, 6],
        vec![4, 3, 8, 6, 5],
        vec![2, 5, 3, 7, 6],
        vec![4, 8, 2, 3, 6],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.39346963278427133,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed5x5() {
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
        0.22200753799676337,
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
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.9712771262351103,
        epsilon = 0.000001
    ));
}

#[test]
#[ignore]
fn fixed5x5_large() {
    let input = vec![
        vec![8, 8, 3, 5, 2],
        vec![5, 3, 3, 0, 2],
        vec![8, 9, 9, 0, 0],
        vec![9, 4, 5, 3, 2],
        vec![4, 6, 6, 1, 0],
    ];
    let result = calculate(input).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.26314046636138944,
        epsilon = 0.000001
    ));
}

#[test]
#[ignore]
fn fixed9x7() {
    let input = vec![
        vec![0, 0, 2, 0, 0, 0, 1],
        vec![0, 1, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 0, 2],
        vec![1, 2, 0, 0, 1, 0, 0],
        vec![0, 0, 0, 0, 0, 0, 1],
        vec![0, 0, 0, 0, 0, 0, 0],
        vec![1, 2, 0, 0, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 0, 1],
        vec![0, 1, 0, 0, 2, 1, 0],
    ];
    let output = calculate(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.11590654664515711,
        epsilon = 0.000001
    ));
}

#[test]
fn fixed10x8_error() {
    let input = vec![
        vec![0, 4, 1, 0, 0, 0, 1, 0],
        vec![0, 1, 0, 0, 0, 0, 0, 0],
        vec![0, 1, 0, 0, 0, 0, 0, 0],
        vec![1, 8, 1, 0, 1, 0, 0, 0],
        vec![0, 1, 1, 1, 0, 1, 0, 0],
        vec![0, 5, 0, 0, 1, 0, 0, 0],
        vec![1, 3, 0, 1, 2, 2, 1, 0],
        vec![2, 7, 0, 0, 1, 4, 1, 1],
        vec![0, 1, 0, 0, 0, 0, 0, 0],
        vec![0, 1, 1, 0, 0, 1, 0, 0],
    ];
    let result = calculate(input).unwrap();
    dbg!(result);
    assert_eq!(result, -1.0);
}
