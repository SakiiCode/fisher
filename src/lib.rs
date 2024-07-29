use math::Quotient;
use pyo3::prelude::*;
use rand::Rng;
use rayon::iter::{IntoParallelIterator, ParallelIterator};
use std::{cmp::min, vec};

mod asa159;
mod asa643;
mod math;

macro_rules! idx {
    ($arr:expr, $r:expr, $c:expr, $cols:expr) => {
        $arr[$r * $cols + $c]
    };
}

fn fill(mat_new: &mut Vec<u32>, r_sum: &Vec<u32>, c_sum: &Vec<u32>, p_0: f64) -> f64 {
    let r = r_sum.len();
    let c = c_sum.len();
    //print!("{:?} -> ", &mat_new);

    for i in 0..r - 1 {
        let mut temp = r_sum[i];
        for j in 0..c - 1 {
            temp -= idx!(mat_new, i, j, c);
        }
        idx!(mat_new, i, c - 1, c) = temp;
    }
    for j in 0..c - 1 {
        let mut temp = c_sum[j];
        for i in 0..r - 1 {
            temp -= idx!(mat_new, i, j, c);
        }
        idx!(mat_new, r - 1, j, c) = temp;
    }

    let mut temp = r_sum[r - 1];
    for j in 0..c - 1 {
        if temp < idx!(mat_new, r - 1, j, c) {
            //println!();
            return 0.0;
        } else {
            temp -= idx!(mat_new, r - 1, j, c);
        }
    }

    idx!(mat_new, r - 1, c - 1, c) = temp;
    //print!("{:?} ", &mat_new);

    let mut p_1 = Quotient::default();

    r_sum.iter().for_each(|x| p_1.mul_fact(*x));
    c_sum.iter().for_each(|y| p_1.mul_fact(*y));

    p_1.div_fact(r_sum.iter().sum());
    mat_new.iter().for_each(|i| p_1.div_fact(*i));

    let p_1_res = p_1.solve();
    if p_1_res <= p_0 + 0.00000001 {
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
    let mut max_1 = r_sum[xx];
    let mut max_2 = c_sum[yy];

    for j in 0..c {
        max_1 -= idx!(mat_new, xx, j, c);
    }

    for i in 0..r {
        max_2 -= idx!(mat_new, i, yy, c);
    }

    return (0..=min(max_1, max_2))
        .into_par_iter()
        .map(|k| {
            let mut mat_new2 = mat_new.clone();
            idx!(mat_new2, xx, yy, c) = k;
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

/// Formats the sum of two numbers as string.
#[pyfunction]
pub fn recursive(table: Vec<Vec<u32>>) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let mut col_sum: Vec<u32> = vec![];

    for j in 0..table[0].len() {
        let mut temp = 0;
        for i in 0..table.len() {
            temp += table[i][j];
        }
        col_sum.push(temp);
    }

    let mut mat = vec![0; col_sum.len() * row_sum.len()];

    let mut p_0 = Quotient::default();

    row_sum.iter().for_each(|x| p_0.mul_fact(*x));
    col_sum.iter().for_each(|y| p_0.mul_fact(*y));

    p_0.div_fact(row_sum.iter().sum());
    table.iter().flatten().for_each(|i| p_0.div_fact(*i));

    let p = _dfs(&mut mat, 0, 0, &row_sum, &col_sum, p_0.solve());

    return Ok(p);
}

#[pyfunction]
pub fn sim(table: Vec<Vec<u32>>, iterations: u32) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let mut col_sum: Vec<u32> = vec![];

    for j in 0..table[0].len() {
        let mut temp = 0;
        for i in 0..table.len() {
            temp += table[i][j];
        }
        col_sum.push(temp);
    }

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
    let statistic = find_statistic_r(&table, &fact);

    let test = generate(&row_sum_i, &col_sum_i, &fact);
    if test.is_err() {
        println!("{}", test.as_ref().unwrap_err().1);
        return Ok(-f64::from(test.unwrap_err().0));
    }

    // STATISTIC <- -sum(lfactorial(x))
    let sum = (0..iterations)
        .into_par_iter()
        .map(|_| generate(&row_sum_i, &col_sum_i, &fact))
        .map(|table| find_statistic_c(&table.unwrap(), nrow, ncol, &fact))
        .filter(|ans| *ans <= statistic + 0.00000001)
        .count() as f64;

    // PVAL <- (1 + sum(tmp <= STATISTIC/almost.1)) / (B + 1)
    let pvalue = (1.0 + sum) / (iterations as f64 + 1.0);

    return Ok(pvalue);
}

#[pyfunction]
#[pyo3(signature = (table, workspace=None))]
pub fn exact(table: Vec<Vec<u32>>, workspace: Option<u32>) -> PyResult<f64> {
    let row_sum: Vec<u32> = table.iter().map(|row| row.iter().sum()).collect();
    let mut col_sum: Vec<u32> = vec![];

    let mut seq = vec![];
    for j in 0..table[0].len() {
        let mut temp = 0;
        for i in 0..table.len() {
            temp += table[i][j];
            seq.push(table[i][j] as f64);
        }
        col_sum.push(temp);
    }

    let wsize;
    if workspace.is_some() {
        wsize = workspace.unwrap();
    } else {
        let sum: u32 = row_sum.iter().sum();
        let exp = sum / 20;
        wsize = 200 * 10u32.pow(exp.clamp(3, 6));
    }
    //dbg!(wsize);

    let result;
    let code;
    unsafe {
        let mut nrow = [row_sum.len() as i32; 1];
        let mut ncol = [col_sum.len() as i32; 1];
        let mut expect = [0.0];
        let mut percnt = [0.0];
        let mut emin = [0.0];
        let mut prt = [0.0];
        let mut pre = [0.0];
        let mut ws = [wsize.try_into().unwrap()];
        code = asa643::fexact_(
            nrow.as_mut_ptr(),
            ncol.as_mut_ptr(),
            seq.as_mut_ptr(),
            nrow.as_mut_ptr(),
            expect.as_mut_ptr(),
            percnt.as_mut_ptr(),
            emin.as_mut_ptr(),
            prt.as_mut_ptr(),
            pre.as_mut_ptr(),
            ws.as_mut_ptr(),
        );

        result = pre[0];
    }
    if code < 0 {
        return Ok(f64::try_from(code).unwrap());
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
        &fact,
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
    assert_eq!(output, 0.5920745920745922);
}

#[test]
fn rec3x3() {
    let input = vec![vec![4, 1, 0], vec![1, 5, 0], vec![1, 1, 4]];
    let output = recursive(input).unwrap();
    dbg!(output);
    assert_eq!(output, 0.005293725881961175);
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
    assert_eq!(output, 0.01096124432190708);
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
    assert_eq!(output, 0.6388806191300103);
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
    assert_eq!(result, 0.24678711559405733);
}

#[test]
fn proc2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = exact(input, None).unwrap();
    dbg!(output);
    assert_eq!(output, 0.5920745920745918);
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
        0.01096,
        epsilon = 0.0001
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
    assert!(float_cmp::approx_eq!(f64, output, 0.0, epsilon = 0.0001));
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
        0.8821660735808745,
        epsilon = 0.0001
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
        0.39346963278449454,
        epsilon = 0.0001
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
        0.22200753799676035,
        epsilon = 0.0001
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
    assert!(float_cmp::approx_eq!(f64, result, 0.9239, epsilon = 0.0001));
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
    assert!(sim(input, 10000).is_err());
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
