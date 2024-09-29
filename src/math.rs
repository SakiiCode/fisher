#[derive(Default)]
pub struct Quotient {
    numerator: Vec<i32>,
    denominator: Vec<i32>,
    initial_n: usize,
    initial_d: usize,
}

// const MAX_FACTORIAL: usize = 13;
// const FACTORIALS: [i32; MAX_FACTORIAL] = [
//     1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600,
// ];

impl Quotient {
    pub fn new(size_n: usize, size_d: usize, init_n: &[i32], init_d: &[i32]) -> Quotient {
        let mut result = Quotient {
            numerator: Vec::with_capacity(size_n),
            denominator: Vec::with_capacity(size_d),
            initial_n: 0,
            initial_d: 0,
        };
        result.mul_fact(init_n);
        result.div_fact(init_d);
        result.initial_n = result.numerator.len();
        result.initial_d = result.denominator.len();
        result
    }

    #[inline(never)]
    pub fn mul_fact(&mut self, arr: &[i32]) {
        for x in arr {
            self.numerator.extend(2..=*x);
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
            // for i in 2..=*x {
            //     self.denominator.push(i);
            // }
            self.denominator.extend(2..=*x);
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
        let mut result = 1.0;

        let n = self.numerator.len();
        let d = self.denominator.len();

        //let len = usize::min(n, d);
        let len = if n < d { n } else { d };

        for i in 0..len {
            result *= self.numerator[i] as f64 / self.denominator[i] as f64;
        }

        if n > d {
            for i in d..n {
                result *= self.numerator[i] as f64;
            }
        } else if n < d {
            for i in n..d {
                result /= self.denominator[i] as f64;
            }
        }
        return result;
    }

    #[inline(never)]
    pub fn clear(&mut self) {
        self.numerator.truncate(self.initial_n);
        self.denominator.truncate(self.initial_d);
    }
}

#[test]
fn test1() {
    let mut q = Quotient::default();
    q.mul_fact(&[15, 7, 7]);
    q.div_fact(&[4, 16, 7, 6]);

    assert!(float_cmp::approx_eq!(
        f64,
        q.solve(),
        7.0 / (4.0 * 3.0 * 2.0 * 16.0)
    ));
}

/*
#[test]
fn test2() {
    let mut q = Quotient::default();
    q.mul_fact(7);
    q.mul_fact(6);
    q.mul_fact(6);
    q.mul_fact(6);

    q.div_fact(13);
    q.div_fact(6);
    q.div_fact(6);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 1.0 / 1716.0));
}

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
