#![allow(clippy::needless_return)]
#![allow(clippy::ptr_arg)]
use lazy_static::lazy_static;
use math::Quotient;
use pyo3::prelude::*;
use rand::Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{cmp::min, sync::Mutex, vec};

mod asa159;
mod asa643;
mod math;

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
fn fill(mat_new: &mut Vec<u32>, r_sum: &Vec<u32>, c_sum: &Vec<u32>, p_0: f64) -> f64 {
    let r = r_sum.len();
    let c = c_sum.len();
    //print!("{:?} -> ", &mat_new);

    for j in 0..c - 1 {
        let mut temp = c_sum[j];
        for i in 0..r - 1 {
            temp -= get!(mat_new, i, j, c);
        }
        set!(mat_new, r - 1, j, c, temp);
    }

    let temp = r_sum[r - 1];
    let mut sum = 0;
    for j in 0..c - 1 {
        sum += get!(mat_new, r - 1, j, c);
    }
    if temp < sum {
        //println!();
        return 0.0;
    }

    set!(mat_new, r - 1, c - 1, c, temp - sum);
    //print!("{:?} ", &mat_new);

    for i in 0..r - 1 {
        let mut temp = r_sum[i];
        for j in 0..c - 1 {
            temp -= get!(mat_new, i, j, c);
        }
        set!(mat_new, i, c - 1, c, temp);
    }

    let n = r_sum.iter().sum();

    let mut p_1 = Quotient::new(2 * n as usize, 2 * n as usize);

    p_1.mul_fact(r_sum);
    p_1.mul_fact(c_sum);

    p_1.div_fact(&[n; 1]);
    p_1.div_fact(mat_new);

    let p_1_res = p_1.solve();
    if p_1_res <= p_0 {
        //println!(" p={}", p_1_res);
        return p_1_res;
    } else {
        //println!(" p=0.0");
        return 0.0;
    }
}

fn _dfs(
    mat_new: &mut Vec<u32>,
    xx: usize,
    yy: usize,
    r_sum: &Vec<u32>,
    c_sum: &Vec<u32>,
    p_0: f64,
) -> f64 {
    let r = r_sum.len();
    let c = c_sum.len();
    let mut max_1;
    let mut max_2;
    unsafe {
        max_1 = *r_sum.get_unchecked(xx);
        max_2 = *c_sum.get_unchecked(yy);
    }

    for j in 0..c {
        max_1 -= get!(mat_new, xx, j, c);
    }

    for i in 0..r {
        max_2 -= get!(mat_new, i, yy, c);
    }

    return (0..=min(max_1, max_2))
        .into_par_iter()
        .map(|k| {
            let mut mat_new2 = mat_new.clone();
            set!(mat_new2, xx, yy, c, k);
            if xx + 2 == r && yy + 2 == c {
                return fill(&mut mat_new2, r_sum, c_sum, p_0);
            } else if xx + 2 == r {
                return _dfs(&mut mat_new2, 0, yy + 1, r_sum, c_sum, p_0);
            } else {
                return _dfs(&mut mat_new2, xx + 1, yy, r_sum, c_sum, p_0);
            }
        })
        .sum();
}

#[pyfunction]
pub fn recursive(table: Vec<Vec<u32>>) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let col_sum: Vec<u32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    let mut mat = vec![0; col_sum.len() * row_sum.len()];

    let mut p_0 = Quotient::default();

    p_0.mul_fact(&row_sum);
    p_0.mul_fact(&col_sum);

    p_0.div_fact(&[row_sum.iter().sum(); 1]);
    p_0.div_fact(&table.iter().flatten().map(|x| *x).collect::<Vec<u32>>());

    let p = _dfs(&mut mat, 0, 0, &row_sum, &col_sum, p_0.solve() + 0.00000001);

    return Ok(p);
}

#[pyfunction]
pub fn sim(table: Vec<Vec<u32>>, iterations: u32) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let col_sum: Vec<u32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    let n = row_sum.iter().sum::<u32>().try_into().unwrap();

    let mut fact = vec![0.0; n + 1];
    fact[0] = 0.0;
    let mut x = 1.0;
    for i in 1..=n {
        fact[i] = fact[i - 1] + f64::ln(x);
        x += 1.0;
    }

    let row_sum_i = row_sum.iter().map(|s| (*s).try_into().unwrap()).collect();
    let col_sum_i = col_sum.iter().map(|s| (*s).try_into().unwrap()).collect();

    let nrow = row_sum.len();
    let ncol = col_sum.len();
    let statistic = find_statistic_r(&table, &fact) + 0.00000001;

    let test = generate(&row_sum_i, &col_sum_i, &fact);
    if let Err(error) = test {
        println!("{}", error.1);
        return Ok(-f64::from(error.0));
    }

    // STATISTIC <- -sum(lfactorial(x))
    let sum = (0..iterations)
        .into_par_iter()
        .map(|_| generate(&row_sum_i, &col_sum_i, &fact))
        .map(|table| find_statistic_c(&table.unwrap(), nrow, ncol, &fact))
        .filter(|ans| *ans <= statistic)
        .count() as f64;

    // PVAL <- (1 + sum(tmp <= STATISTIC/almost.1)) / (B + 1)
    let pvalue = (1.0 + sum) / (iterations as f64 + 1.0);

    return Ok(pvalue);
}

lazy_static! {
    static ref FEXACT_LOCK: Mutex<()> = Mutex::new(());
}

#[pyfunction]
#[pyo3(signature = (table, workspace=None))]
pub fn exact(table: Vec<Vec<u32>>, workspace: Option<u32>) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let col_sum: Vec<u32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    // seq needs to be column-major
    let mut seq: Vec<f64> = (0..(table[0].len()))
        .flat_map(|index| table.iter().map(move |row| row[index] as f64))
        .collect();

    let wsize = match workspace {
        Some(size) => size,
        None => {
            let sum: u32 = row_sum.iter().sum();
            let exp = sum / 20;
            200 * 10u32.pow(exp.clamp(3, 6))
        }
    };
    //dbg!(wsize);

    let result;
    let code;
    unsafe {
        let _guard = FEXACT_LOCK.lock();
        let nrow = row_sum.len() as i32;
        let ncol = col_sum.len() as i32;
        let mut expect = 0.0;
        let mut percnt = 0.0;
        let mut emin = 0.0;
        let mut prt = 0.0;
        let mut pre = 0.0;
        let ws = wsize.try_into().unwrap();
        code = asa643::fexact_(
            nrow,
            ncol,
            seq.as_mut_ptr(),
            nrow,
            &mut expect,
            &mut percnt,
            &mut emin,
            &mut prt,
            &mut pre,
            ws,
        );

        result = pre;
    }
    if code < 0 {
        return Ok(f64::from(code));
    } else {
        return Ok(result);
    }
}

fn find_statistic_c(table: &Vec<i32>, nrow: usize, ncol: usize, fact: &Vec<f64>) -> f64 {
    let mut ans = 0.0;
    for i in 0..nrow {
        for j in 0..ncol {
            ans -= fact[table[i * ncol + j] as usize];
        }
    }
    return ans;
}

fn find_statistic_r(table: &Vec<Vec<u32>>, fact: &Vec<f64>) -> f64 {
    let mut ans = 0.0;
    for row in table {
        for cell in row {
            ans -= fact[*cell as usize];
        }
    }
    return ans;
}

fn generate(
    row_sum: &Vec<i32>,
    col_sum: &Vec<i32>,
    fact: &Vec<f64>,
) -> Result<Vec<i32>, (i32, &'static str)> {
    let mut rng = rand::thread_rng();
    let mut seed = rng.gen::<i32>();

    let result = asa159::rcont2(
        i32::try_from(row_sum.len()).unwrap(),
        i32::try_from(col_sum.len()).unwrap(),
        row_sum,
        col_sum,
        &mut 0,
        &mut seed,
        fact,
    );

    return result;
}

/// A Python module implemented in Rust.
#[pymodule]
fn fisher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(recursive, m)?)?;
    m.add_function(wrap_pyfunction!(sim, m)?)?;
    m.add_function(wrap_pyfunction!(exact, m)?)?;
    Ok(())
}

#[test]
fn rec2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.5920745920745921);
}

#[test]
fn rec3x3() {
    let input = vec![vec![4, 1, 0], vec![1, 5, 0], vec![1, 1, 4]];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.005293725881961177);
}

#[test]
fn rec4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.010961244321907074);
}

#[test]
#[ignore]
fn rec4x4big() {
    let input = vec![
        vec![28, 28, 0, 28],
        vec![0, 0, 16, 0],
        vec![0, 0, 5, 0],
        vec![0, 0, 7, 0],
    ];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.00594949486831653); // this should be 0.0
}

#[test]
#[ignore]
fn rec3x4big() {
    // DNF
    let input = vec![
        vec![11, 12, 18, 15],
        vec![15, 13, 13, 15],
        vec![15, 19, 19, 15],
    ];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.8823465060350275,
        epsilon = 0.0001
    ));
}

#[test]
fn rec5x4() {
    // 0.11s
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
    ];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.63888061913001);
}

#[test]
#[ignore]
fn rec5x5() {
    // 12.13s
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
        vec![1, 1, 0, 0, 3],
    ];
    let result = recursive(input).unwrap();
    dbg!(result);
    assert_eq!(result, 0.24678711559405725);
}

#[test]
fn proc2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.5920745920745918,
        epsilon = 0.000001
    ));
}

#[test]
fn proc4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.01096124432190692,
        epsilon = 0.000001
    ));
}

#[test]
fn proc4x4big() {
    let input = vec![
        vec![28, 28, 28, 0],
        vec![0, 0, 0, 16],
        vec![0, 0, 0, 5],
        vec![0, 0, 0, 7],
    ];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, 0.0, epsilon = 0.000001));
}

#[test]
fn proc3x4big() {
    let input = vec![
        vec![11, 12, 18, 15],
        vec![15, 13, 13, 15],
        vec![15, 19, 19, 15],
    ];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.8821660735808727,
        epsilon = 0.000001
    ));
}

#[test]
fn proc4x5big() {
    let input = vec![
        vec![8, 3, 5, 5, 6],
        vec![4, 3, 8, 6, 5],
        vec![2, 5, 3, 7, 6],
        vec![4, 8, 2, 3, 6],
    ];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.39346963278427133,
        epsilon = 0.000001
    ));
}

#[test]
fn proc5x5() {
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
        vec![1, 1, 0, 0, 3],
    ];
    let result = exact(input, None).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.22200753799676337,
        epsilon = 0.000001
    ));
}

#[test]
#[ignore]
fn proc5x5_large() {
    let input = vec![
        vec![8, 8, 3, 5, 2],
        vec![5, 3, 3, 0, 2],
        vec![8, 9, 9, 0, 0],
        vec![9, 4, 5, 3, 2],
        vec![4, 6, 6, 1, 0],
    ];
    let result = exact(input, None).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.26314046636138944,
        epsilon = 0.00001
    ));
}

#[test]
#[ignore]
fn proc7x4_2e8() {
    let input = vec![
        vec![41, 22, 18, 5],
        vec![5, 3, 3, 0],
        vec![20, 9, 9, 0],
        vec![10, 4, 5, 3],
        vec![16, 6, 6, 1],
        vec![13, 8, 5, 2],
        vec![19, 12, 12, 6],
    ];
    let result = exact(input, None).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.9239149531167176,
        epsilon = 0.00001
    ));
}

#[test]
fn proc4x15() {
    let input = vec![
        vec![23, 22, 13, 22, 19, 16, 22, 22, 24, 20, 14, 16, 19, 16, 19],
        vec![26, 20, 6, 20, 13, 12, 21, 18, 19, 14, 14, 14, 18, 11, 14],
        vec![26, 22, 14, 22, 14, 17, 22, 21, 23, 23, 14, 18, 16, 12, 13],
        vec![26, 23, 13, 24, 18, 19, 24, 25, 22, 18, 18, 17, 21, 21, 18],
    ];
    let result = exact(input, Some(200000000)).unwrap();
    dbg!(result);
    assert_eq!(result, -501.0);
}

#[test]
fn proc16x8() {
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
        vec![0, 3, 1, 0, 0, 0, 1, 0],
        vec![0, 0, 0, 0, 1, 0, 0, 0],
        vec![0, 0, 0, 0, 0, 3, 0, 0],
        vec![0, 1, 0, 0, 0, 0, 0, 0],
        vec![1, 2, 1, 1, 0, 1, 0, 1],
        vec![1, 0, 1, 0, 1, 3, 0, 0],
    ];
    let result = exact(input, Some(200000000)).unwrap();
    dbg!(result);
    assert_eq!(result, -502.0);
}

#[test]
fn sim2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let result = sim(input, 10000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.592, epsilon = 0.02));
}

#[test]
fn sim4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let result = sim(input, 10000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.011, epsilon = 0.004));
}

#[test]
fn sim4x4error() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![0, 0, 0, 0],
        vec![1, 1, 0, 3],
    ];
    assert!(sim(input, 10000).unwrap() < 0.0);
}

#[test]
fn sim4x4big() {
    let input = vec![
        vec![28, 28, 28, 0],
        vec![0, 0, 0, 16],
        vec![0, 0, 0, 5],
        vec![0, 0, 0, 7],
    ];
    let result = sim(input, 100000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.0, epsilon = 0.004));
}

#[test]
fn sim5x5() {
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
        vec![1, 1, 0, 0, 3],
    ];
    let result = sim(input, 10000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.222, epsilon = 0.02));
}
