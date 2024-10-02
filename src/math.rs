#[cfg(test)]
use crate::asa159;
#[cfg(test)]
use rand::Rng;

pub struct Quotient {
    container: Vec<f64>,
    size: usize,
    initial_sln: f64,
    initial_idx: usize,
    //n_idx: usize,
    //d_idx: usize,
    idx: usize,
    solution: f64,
}

// const MAX_FACTORIAL: usize = 13;
// const FACTORIALS: [i32; MAX_FACTORIAL] = [
//     1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
// ];

impl Quotient {
    pub fn new(n: usize, init_n: &[i32], init_d: &[i32]) -> Quotient {
        let size = 2 * n;
        let container_size = size;
        let container = Vec::with_capacity(container_size);
        let mut result: Quotient = Quotient {
            container,
            size,
            initial_sln: 1.0,
            initial_idx: 0,
            //n_idx: 0,
            //d_idx: size,
            idx: 0,
            solution: 1.0,
        };
        result.mul_fact(init_n);
        result.div_fact(init_d);
        //result.initial_n = init_n.iter().map(|i| usize::try_from(*i).unwrap()).sum();
        result.initial_idx = init_d.iter().map(|i| usize::try_from(*i).unwrap()).sum();
        result.initial_sln = result.solution;
        result
    }

    #[inline(never)]
    pub fn mul_fact(&mut self, arr: &[i32]) {
        for x in arr {
            // for i in 1..=*x {
            //     unsafe {
            //         *self.container.get_unchecked_mut(self.idx) = i as f64;
            //     }
            //     //self.idx += 1;
            // }
            self.container.extend((1..=*x).map(|x| x as f64));
            // let idx = *x as usize;
            // if idx < MAX_FACTORIAL {
            //     self.numerator.push(FACTORIALS[idx]);
            // } else {
            //     self.numerator.push(FACTORIALS[MAX_FACTORIAL - 1]);
            //     self.numerator.extend((MAX_FACTORIAL as i32)..=*x);
            // }
        }
    }

    #[inline(never)]
    pub fn div_fact(&mut self, arr: &[i32]) {
        for x in arr {
            for i in 1..=*x {
                unsafe {
                    //*self.container.get_unchecked_mut(self.d_idx) = i as f64;
                    let num = *self.container.get_unchecked(self.idx);
                    self.solution *= num / i as f64;
                }

                self.idx += 1;
            }
            //self.denominator.extend((1..=*x).map(|x| x as f64));
            // let idx = *x as usize;
            // if idx < MAX_FACTORIAL {
            //     self.denominator.push(FACTORIALS[idx]);
            // } else {
            //     self.denominator.push(FACTORIALS[MAX_FACTORIAL - 1]);
            //     self.denominator.extend((MAX_FACTORIAL as i32)..=*x);
            // }
        }
    }

    #[inline(never)]
    pub fn solve(&mut self) -> f64 {
        //let mut result = 1.0;

        //let n = self.numerator.len();
        //let d = self.denominator.len();

        //let len = usize::min(n, d);
        //let len = if n < d { n } else { d };
        //assert!(n == d);
        //let n = self.n_idx;
        //let d = self.d_idx;

        //let len = if n < d { n } else { d };
        // TODO without unsafe
        /*unsafe {
            //for i in 0..d {
            //     result *=
            //         self.container.get_unchecked(i) / self.container.get_unchecked(self.size + i);
            // }
            //if n > d {
            for i in self.idx..self.size {
                self.solution *= self.container.get_unchecked(i);
            }
            // } else {
            //     for i in d..n {
            //         result /= self.container.get_unchecked(self.size + i);
            //     }
            // }
        }*/

        return self.solution;
    }

    #[inline(never)]
    pub fn clear(&mut self) {
        //self.n_idx = self.initial_n;
        //self.d_idx = self.size; // + self.initial_d;
        self.idx = self.initial_idx;
        self.solution = self.initial_sln;
        //self.numerator.truncate(self.initial_n);
        //self.denominator.truncate(self.initial_d);
    }
}

#[test]
fn test1() {
    let row_sum = vec![4, 5, 3, 3, 5];
    let col_sum = vec![3, 7, 2, 3, 5];

    let n: i32 = row_sum.iter().sum();

    let mut q = Quotient::new(n as usize, &[], &[]);

    let mut fact = vec![0.0; (n + 1) as usize];
    fact[0] = 0.0;
    let mut x = 1.0;
    for i in 1..=(n as usize) {
        fact[i] = fact[i - 1] + f64::ln(x);
        x += 1.0;
    }

    let mut rng = rand::thread_rng();
    let mut seed = rng.gen::<i32>();

    let result = asa159::rcont2(
        i32::try_from(row_sum.len()).unwrap(),
        i32::try_from(col_sum.len()).unwrap(),
        &row_sum,
        &col_sum,
        &mut 0,
        &mut seed,
        &fact,
    );

    q.mul_fact(&[4, 5, 3, 3, 3, 7, 2, 3]);
    q.div_fact(&[n; 1]);
    println!("{:?}", &result);
    q.div_fact(&result.unwrap());

    let result = q.solve();
    dbg!(&result);

    //assert!(float_cmp::approx_eq!(f64, result, 6.0 / 5005.0));
}

#[test]
fn test2() {
    let mut q = Quotient::new(19, &[7, 6, 6, 6], &[13]);

    q.div_fact(&[6]);
    q.div_fact(&[6]);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 1.0 / 1716.0));
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
