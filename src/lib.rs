use factorial::Factorial;
use pyo3::prelude::*;
use std::cmp::min;

struct Pos(i32, i32);

macro_rules! fac {
    ($x:expr) => {
        f64::try_from($x.checked_factorial().unwrap()).unwrap()
    };
}

fn _dfs(
    mat: Vec<Vec<u32>>,
    pos: Pos,
    r_sum: &Vec<u32>,
    c_sum: &Vec<u32>,
    p_0: f64,
    p: &mut Vec<f64>,
) {
    let Pos { 0: xx, 1: yy } = pos;
    let r = r_sum.len();
    let c = c_sum.len();

    let mut mat_new = vec![];

    for i in 0..mat.len() {
        let mut temp = vec![];
        for j in 0..mat[0].len() {
            temp.push(mat[i][j]);
        }
        mat_new.push(temp);
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
        let mut temp: i32 = r_sum[r - 1].try_into().unwrap();
        for j in 0..c - 1 {
            temp -= i32::try_from(mat_new[r - 1][j]).unwrap();
        }
        if temp < 0 {
            return;
        }
        mat_new[r - 1][c - 1] = temp.try_into().unwrap();

        let mut p_1 = 1.0;

        for x in r_sum.clone() {
            p_1 *= fac!(x);
        }
        for y in c_sum {
            p_1 *= fac!(y);
        }

        let mut n = 0;
        for x in r_sum {
            n += x;
        }
        p_1 /= fac!(n);

        for i in 0..mat_new.len() {
            for j in 0..mat_new[0].len() {
                p_1 /= fac!(mat_new[i][j]);
            }
        }

        if p_1 <= p_0 + 0.00000001 {
            p[0] += p_1;
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
            max_2 = mat_new[i][y_idx];
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

            _dfs(mat_new.clone(), pos_new, r_sum, c_sum, p_0, p);
        }
    }
}

/// Formats the sum of two numbers as string.
#[pyfunction]
fn exact(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.
#[pymodule]
fn fisher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(exact, m)?)?;
    Ok(())
}
