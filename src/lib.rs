use math::Quotient;
use pyo3::prelude::*;
use rand::Rng;
use std::{cmp::min, vec};

mod asa159;
mod math;

struct Pos(i64, i64);

fn _dfs(
    mat: &mut Vec<Vec<u64>>,
    pos: Pos,
    r_sum: &Vec<u64>,
    c_sum: &Vec<u64>,
    p_0: f64,
    p: &mut Vec<f64>,
) {
    let Pos { 0: xx, 1: yy } = pos;
    let r = r_sum.len();
    let c = c_sum.len();

    let mut mat_new = vec![];

    for i in 0..mat.len() {
        mat_new.push(mat[i].clone());
    }

    if xx == -1 && yy == -1 {
        for i in 0..r - 1 {
            let mut temp = r_sum[i];
            for j in 0..c - 1 {
                temp -= mat_new[i][j];
            }
            mat_new[i][c - 1] = temp;
        }
        for j in 0..c - 1 {
            let mut temp = c_sum[j];
            for i in 0..r - 1 {
                temp -= mat_new[i][j];
            }
            mat_new[r - 1][j] = temp;
        }

        let mut temp = r_sum[r - 1];
        for j in 0..c - 1 {
            if temp < mat_new[r - 1][j] {
                return;
            } else {
                temp -= mat_new[r - 1][j];
            }
        }

        mat_new[r - 1][c - 1] = temp.try_into().unwrap();

        let mut p_1 = Quotient::default();

        let mut n = 0;

        for x in r_sum.clone() {
            p_1.mul_fact(x);
            n += x;
        }
        for y in c_sum.clone() {
            p_1.mul_fact(y);
        }
        p_1.div_fact(n);

        for i in 0..mat_new.len() {
            for j in 0..mat_new[0].len() {
                p_1.div_fact(mat_new[i][j]);
            }
        }

        let p_1_res = p_1.solve();
        if p_1_res <= p_0 + 0.00000001 {
            p[0] += p_1_res;
        }
    } else {
        let x_idx: usize = xx.try_into().unwrap();
        let y_idx: usize = yy.try_into().unwrap();

        let mut max_1 = r_sum[x_idx];
        let mut max_2 = c_sum[y_idx];

        for j in 0..c {
            max_1 -= mat_new[x_idx][j];
        }
        for i in 0..r {
            max_2 -= mat_new[i][y_idx];
        }
        for k in 0..=min(max_1, max_2) {
            mat_new[x_idx][y_idx] = k;
            let pos_new: Pos;
            if x_idx + 2 == r && y_idx + 2 == c {
                pos_new = Pos(-1, -1);
            } else if x_idx + 2 == r {
                pos_new = Pos(0, yy + 1);
            } else {
                pos_new = Pos(xx + 1, yy);
            }

            _dfs(&mut mat_new, pos_new, r_sum, c_sum, p_0, p);
        }
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
pub fn exact(table: Vec<Vec<u64>>) -> PyResult<f64> {
    let mut row_sum: Vec<u64> = vec![];
    let mut col_sum: Vec<u64> = vec![];

    for i in 0..table.len() {
        let mut temp = 0u64;
        for j in 0..table[0].len() {
            temp += table[i][j];
        }
        row_sum.push(temp);
    }

    for j in 0..table[0].len() {
        let mut temp = 0;
        for i in 0..table.len() {
            temp += table[i][j];
        }
        col_sum.push(temp);
    }

    let mut mat = vec![vec![0; col_sum.len()]; row_sum.len()];
    let pos = Pos(0, 0);

    let mut p_0 = Quotient::default();
    let mut n = 0;

    for x in row_sum.clone() {
        p_0.mul_fact(x);
        n += x;
    }
    for y in col_sum.clone() {
        p_0.mul_fact(y);
    }

    p_0.div_fact(n);

    for i in 0..table.len() {
        for j in 0..table[0].len() {
            p_0.div_fact(table[i][j]);
        }
    }

    let mut p = vec![0.0];
    _dfs(&mut mat, pos, &row_sum, &col_sum, p_0.solve(), &mut p);

    return Ok(p[0]);
}

#[pyfunction]
pub fn sim(table: Vec<Vec<u32>>, iterations: u32) -> PyResult<f64> {
    let mut row_sum: Vec<u32> = vec![];
    let mut col_sum: Vec<u32> = vec![];

    for i in 0..table.len() {
        let mut temp = 0u32;
        for j in 0..table[0].len() {
            temp += table[i][j];
        }
        row_sum.push(temp);
    }

    for j in 0..table[0].len() {
        let mut temp = 0;
        for i in 0..table.len() {
            temp += table[i][j];
        }
        col_sum.push(temp);
    }

    let mut rng = rand::thread_rng();

    let mut results = vec![];

    let mut n: usize = 0;
    for sum in row_sum.clone() {
        n += usize::try_from(sum).unwrap();
    }

    let mut fact = vec![0.0; n + 1];
    fact[0] = 0.0;
    fact[1] = 0.0;
    let mut x = 2.0;
    for i in 2..=n {
        fact[i] = fact[i - 1] + f64::ln(x);
        x += 1.0;
    }

    //dbg!(fact.clone());

    let mut row_sum_i = vec![];
    let mut col_sum_i = vec![];
    for s in row_sum.clone() {
        row_sum_i.push(s.try_into().unwrap());
    }
    for s in col_sum.clone() {
        col_sum_i.push(s.try_into().unwrap());
    }

    for _ in 0..iterations {
        let nrow = row_sum.len();
        let ncol = col_sum.len();
        unsafe {
            let table = generate(&mut row_sum_i, &mut col_sum_i, rng.gen::<i32>());
            //dbg!(table.clone());
            // STATISTIC <- -sum(lfactorial(x))
            let ans = find_statistic_c(table, nrow, ncol, &fact);
            results.push(ans);
        }
    }

    //dbg!(results.clone());

    // PVAL <- (1 + sum(tmp <= STATISTIC/almost.1)) / (B + 1)

    let mut sum = 0.0;
    let statistic = find_statistic_r(table, &fact);

    for i in 0..iterations {
        if results[i as usize] <= statistic + 0.00000001 {
            sum += 1.0;
        }
    }

    let pvalue = (1.0 + sum) / (iterations as f64 + 1.0);

    return Ok(pvalue);
}

fn find_statistic_c(table: Vec<i32>, nrow: usize, ncol: usize, fact: &Vec<f64>) -> f64 {
    let mut ans = 0.0;
    for i in 0..nrow {
        for j in 0..ncol {
            ans -= fact[table[i * ncol + j] as usize];
        }
    }
    return ans;
}

fn find_statistic_r(table: Vec<Vec<u32>>, fact: &Vec<f64>) -> f64 {
    let mut ans = 0.0;
    for row in table {
        for cell in row {
            ans -= fact[cell as usize];
        }
    }
    return ans;
}

unsafe fn generate(row_sum: &Vec<i32>, col_sum: &Vec<i32>, seed: i32) -> Vec<i32> {
    let mut seed = [seed; 1];
    let mut key = [0; 1];
    let size = row_sum.len() * col_sum.len();
    let mut matrix = vec![0; size];
    let mut ierror = [0; 1];
    asa159::rcont2(
        i32::try_from(row_sum.len()).unwrap(),
        i32::try_from(col_sum.len()).unwrap(),
        row_sum.clone().as_mut_ptr(),
        col_sum.clone().as_mut_ptr(),
        key.as_mut_ptr(),
        seed.as_mut_ptr(),
        matrix.as_mut_ptr(),
        ierror.as_mut_ptr(),
    );
    //dbg!(ierror);
    //dbg!(matrix);
    if ierror[0] != 0 {
        panic!("Error generating contingency table ({})", ierror[0]);
    }
    return matrix;
}

/// A Python module implemented in Rust.
#[pymodule]
fn fisher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(exact, m)?)?;
    m.add_function(wrap_pyfunction!(sim, m)?)?;
    Ok(())
}

#[test]
fn size2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = exact(input).unwrap();
    assert_eq!(output, 0.5920745920745921);
}

#[test]
fn size4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let output = exact(input).unwrap();
    assert_eq!(output, 0.01096124432190799);
}

/*
#[test]
fn size5x4() {
    let input = vec![
        vec![3, 1, 1, 1, 0],
        vec![1, 4, 1, 0, 0],
        vec![2, 1, 3, 2, 0],
        vec![1, 1, 1, 2, 0],
    ];
    let output = exact(input).unwrap();
    assert_eq!(output, 0.6388806191300838);
}*/

#[test]
fn sim2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let result = sim(input.clone(), 10000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.592, epsilon = 0.02));
    //sim(input.clone(), 12345, 1).unwrap();
    //assert_eq!(output, 0.5920745920745921);
}

#[test]
fn sim4x4() {
    let input = vec![
        vec![4, 1, 0, 1],
        vec![1, 5, 0, 0],
        vec![1, 1, 4, 2],
        vec![1, 1, 0, 3],
    ];
    let result = sim(input.clone(), 10000).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(f64, result, 0.011, epsilon = 0.004));
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
