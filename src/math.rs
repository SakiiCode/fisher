use std::{
    ops::{AddAssign, DivAssign},
    simd::{num::SimdFloat, Simd},
};

const N_U: usize = 4;
const N_F: f64 = 4.0;
const DEFAULT_SIMD: [f64; 4] = [1.0, 2.0, 3.0, 4.0];

pub struct Quotient {
    container: Vec<f64>,
    initial_sln: f64,
    initial_idx: usize,
    idx: usize,
    solution: f64,
    offset_simd: Simd<f64, N_U>,
}

impl Quotient {
    pub fn new(n: usize, init_n: &[i32], init_d: &[i32]) -> Quotient {
        let size = 2 * n;
        let container = Vec::with_capacity(size);
        let mut result: Quotient = Quotient {
            container,
            initial_sln: 1.0,
            initial_idx: 0,
            idx: 0,
            solution: 1.0,
            offset_simd: Simd::splat(N_F),
        };
        result.mul_fact(init_n);
        result.div_fact(init_d);
        result.initial_idx = init_d.iter().map(|i| usize::try_from(*i).unwrap()).sum();
        result.initial_sln = result.solution;
        result
    }

    pub fn mul_fact(&mut self, arr: &[i32]) {
        for x in arr {
            self.container.extend((1..=*x).map(|x| x as f64));
        }
    }

    pub fn div_fact(&mut self, arr: &[i32]) {
        for x in arr {
            let mut i = 0usize;
            let max = *x as usize;
            let mut d_simd: Simd<f64, N_U> = Simd::from_array(DEFAULT_SIMD);
            while i < max {
                let slice = &self.container[(self.idx + i)..(self.idx + max)];
                let mut num = Simd::load_or(slice, d_simd);
                num.div_assign(d_simd);
                self.solution *= num.reduce_product();
                d_simd.add_assign(self.offset_simd);
                i += N_U;
            }
            self.idx += max;
        }
    }

    pub fn solve(&self) -> f64 {
        self.solution
    }

    pub fn clear(&mut self) {
        self.idx = self.initial_idx;
        self.solution = self.initial_sln;
    }
}

#[test]
fn test1() {
    let mut row_sum = vec![4, 5, 3, 3, 5];
    let col_sum = vec![3, 7, 2, 3, 5];

    let n: i32 = row_sum.iter().sum();

    row_sum.extend_from_slice(&col_sum);
    let mut q = Quotient::new(n as usize, &row_sum, &[n; 1]);

    let table = [
        1, 1, 1, 0, 0, 2, 1, 0, 1, 3, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1,
    ];

    q.div_fact(&table);

    let result = q.solve();

    assert!(float_cmp::approx_eq!(f64, result, 1.85572e-5, epsilon = 1e-9));

    let mut q = Quotient::new(19, &[7, 6, 6, 6], &[13]);

    q.div_fact(&[6]);
    q.div_fact(&[6]);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 1.0 / 1716.0));
}

#[test]
fn test2() {
    let q = Quotient::new(3, &vec![3], &vec![2, 1, 0]);
    let s = q.solve();
    assert_eq!(s, 3.0);
}

#[test]
fn test3() {
    let mut q = Quotient::new(3, &[], &vec![]);
    q.mul_fact(&[3]);
    q.div_fact(&[2, 1, 0]);
    let s = q.solve();
    assert_eq!(s, 3.0);
}

/*
#[test]
fn test3() {
    let mut q = Quotient::default();
    q.mul_fact(5);
    q.mul_fact(6);

    q.div_fact(3);
    q.div_fact(4);
    q.div_fact(6);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 5.0 / 6.0));
}

#[test]
fn test4() {
    let mut q = Quotient::default();
    q.mul_fact(6);
    q.mul_fact(6);
    q.mul_fact(8);
    q.mul_fact(5);
    q.mul_fact(7);
    q.mul_fact(8);
    q.mul_fact(4);
    q.mul_fact(6);

    q.div_fact(25);
    q.div_fact(6);
    q.div_fact(5);
    q.div_fact(3);
    q.div_fact(4);
    q.div_fact(5);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 1.0 / 2629308825.0));
}
*/
