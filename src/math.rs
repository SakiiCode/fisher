#![allow(clippy::comparison_chain)]
#[derive(Default)]
pub struct Quotient {
    numerator: Vec<f64>,
    denominator: Vec<f64>,
}

impl Quotient {
    pub fn new(nsize: usize, dsize: usize) -> Quotient {
        Quotient {
            numerator: Vec::with_capacity(nsize),
            denominator: Vec::with_capacity(dsize),
        }
    }

    pub fn mul_fact(&mut self, arr: &[u32]) {
        for x in arr {
            for i in 2..=*x {
                self.numerator.push(i.into());
            }
        }
    }

    pub fn div_fact(&mut self, arr: &[u32]) {
        for x in arr {
            for i in 2..=*x {
                self.denominator.push(i.into());
            }
        }
    }

    pub fn solve(&mut self) -> f64 {
        let mut result = 1.0;

        let n = self.numerator.len();
        let d = self.denominator.len();

        let len = usize::min(n, d);

        for i in 0..len {
            result *= self.numerator[i] / self.denominator[i];
        }

        if n > d {
            for i in d..n {
                result *= self.numerator[i];
            }
        } else if n < d {
            for i in n..d {
                result /= self.denominator[i];
            }
        }
        return result;
    }
}
/*
#[test]
fn test1() {
    let mut q = Quotient::default();
    q.mul_fact(7);
    q.mul_fact(7);
    q.mul_fact(7);

    q.div_fact(5);
    q.div_fact(7);
    q.div_fact(6);

    assert!(float_cmp::approx_eq!(f64, q.solve(), 6.0 * 7.0 * 7.0));
}

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
