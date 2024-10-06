use core::f64;
use std::convert::Infallible;

use rand::Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};

use crate::asa159;

pub fn calculate(table: Vec<Vec<i32>>, iterations: i32) -> Result<f64, Infallible> {
    let row_sum: Vec<i32> = table.iter().map(|row| row.iter().sum()).collect();
    let col_sum: Vec<i32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    let n = row_sum.iter().sum::<i32>().try_into().unwrap();

    let mut fact = vec![0.0; n + 1];
    fact[0] = 0.0;
    let mut x = 1.0;
    for i in 1..=n {
        fact[i] = fact[i - 1] + f64::ln(x);
        x += 1.0;
    }

    let seq = table.iter().flatten().cloned().collect::<Vec<i32>>();
    if seq.iter().any(|x| *x < 0) {
        println!("ERROR: Negative value in matrix!");
        return Ok(-1.0);
    }

    let nrow = row_sum.len();
    let ncol = col_sum.len();
    let statistic = find_statistic_r(&table, &fact) + f64::EPSILON;

    let row_sum_i = row_sum.iter().map(|s| (*s).try_into().unwrap()).collect();
    let col_sum_i = col_sum.iter().map(|s| (*s).try_into().unwrap()).collect();

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

fn find_statistic_c(table: &Vec<i32>, nrow: usize, ncol: usize, fact: &Vec<f64>) -> f64 {
    let mut ans = 0.0;
    for i in 0..nrow {
        for j in 0..ncol {
            ans -= fact[table[i * ncol + j] as usize];
        }
    }
    return ans;
}

fn find_statistic_r(table: &Vec<Vec<i32>>, fact: &Vec<f64>) -> f64 {
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

#[test]
fn sim1x1_error() {
    let input = vec![vec![5]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert_eq!(output, -1.0);
}

#[test]
fn sim2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.5920745920745918,
        epsilon = 0.002
    ));
}

#[test]
fn sim2x2_error() {
    let input = vec![vec![3, 4], vec![4, -2]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert_eq!(output, -1.0);
}

#[test]
fn sim3x2() {
    let input = vec![vec![32, 10, 20], vec![20, 25, 18]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.009645916798182401,
        epsilon = 0.002
    ));
}

#[test]
fn sim3x3() {
    let input = vec![vec![32, 10, 20], vec![20, 25, 18], vec![11, 17, 14]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.010967949934049852,
        epsilon = 0.002
    ));
}

#[test]
fn sim3x3_unit() {
    let input = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert_eq!(output, 1.0);
}

#[test]
fn sim3x3_zero() {
    let input = vec![vec![0, 0, 0], vec![0, 0, 0], vec![0, 0, 0]];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert_eq!(output, -3.0);
}

#[test]
fn sim3x4_large() {
    let input = vec![
        vec![11, 12, 18, 15],
        vec![15, 13, 13, 15],
        vec![15, 19, 19, 15],
    ];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.8821660735808727,
        epsilon = 0.002
    ));
}

#[test]
fn sim4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.01096124432190692,
        epsilon = 0.002
    ));
}

#[test]
fn sim4x4_error() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![0, 0, 0, 0],
        vec![1, 1, 0, 3],
    ];
    assert!(calculate(input, 1000000).unwrap() < 0.0);
}

#[test]
fn sim4x4_large() {
    let input = vec![
        vec![28, 28, 28, 0],
        vec![0, 0, 0, 16],
        vec![0, 0, 0, 5],
        vec![0, 0, 0, 7],
    ];
    let result = calculate(input, 1000000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.0, epsilon = 0.004));
}

#[test]
fn sim4x5_large() {
    let input = vec![
        vec![8, 3, 5, 5, 6],
        vec![4, 3, 8, 6, 5],
        vec![2, 5, 3, 7, 6],
        vec![4, 8, 2, 3, 6],
    ];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.39346963278427133,
        epsilon = 0.002
    ));
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
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.22200753799676337,
        epsilon = 0.002
    ));
}

#[test]
fn sim5x5_small() {
    let input = vec![
        vec![1, 0, 0, 0, 0],
        vec![1, 1, 0, 1, 0],
        vec![1, 1, 0, 0, 1],
        vec![0, 0, 1, 2, 1],
        vec![1, 1, 2, 1, 1],
    ];
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.9712771262351103,
        epsilon = 0.002
    ));
}

#[test]
fn sim5x5_large() {
    let input = vec![
        vec![8, 8, 3, 5, 2],
        vec![5, 3, 3, 0, 2],
        vec![8, 9, 9, 0, 0],
        vec![9, 4, 5, 3, 2],
        vec![4, 6, 6, 1, 0],
    ];
    let result = calculate(input, 1000000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.26314046636138944,
        epsilon = 0.002
    ));
}

#[test]
fn sim9x7_error() {
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
    let output = calculate(input, 1000000).unwrap();
    dbg!(output);
    assert_eq!(output, -3.0);
}
