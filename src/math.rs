#[derive(Default)]
pub struct Quotient {
    numerator: Vec<u64>,
    denominator: Vec<u64>,
}

impl Quotient {
    pub fn mul_fact(&mut self, x: u64) {
        for i in 2..=x {
            self.numerator.push(i);
        }
    }

    pub fn div_fact(&mut self, x: u64) {
        for i in 2..=x {
            self.denominator.push(i);
        }
    }

    pub fn solve(&self) -> f64 {
        let mut result = 1.0;
        let max_index = std::cmp::min(self.numerator.len(), self.denominator.len());
        for i in 0..max_index {
            result *= self.numerator[i] as f64 / self.denominator[i] as f64;
        }

        if self.numerator.len() == max_index {
            for i in max_index..self.denominator.len() {
                result /= self.denominator[i] as f64;
            }
        } else {
            for i in max_index..self.numerator.len() {
                result *= self.numerator[i] as f64;
            }
        }
        result
    }
}

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
