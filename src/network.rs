use std::{convert::Infallible, sync::Mutex};

use lazy_static::lazy_static;

use crate::asa643;

lazy_static! {
    static ref FEXACT_LOCK: Mutex<()> = Mutex::new(());
}

pub fn calculate(table: Vec<Vec<i32>>, workspace: Option<i32>) -> Result<f64, Infallible> {
    let row_sum: Vec<i32> = table.iter().map(|row| row.iter().sum()).collect();
    let col_sum: Vec<i32> = (0..(table[0].len()))
        .map(|index| table.iter().map(|row| row[index]).sum())
        .collect();

    // seq needs to be column-major
    let mut seq: Vec<f64> = (0..(table[0].len()))
        .flat_map(|index| table.iter().map(move |row| row[index] as f64))
        .collect();

    let wsize = match workspace {
        Some(size) => size,
        None => {
            let sum: u32 = row_sum.iter().sum::<i32>() as u32;
            let exp = sum / 20;
            (200 * 10i32.pow(exp.clamp(3, 6))).into()
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
            nrow.into(),
            ncol.into(),
            seq.as_mut_ptr(),
            nrow.into(),
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

#[test]
fn proc1x1_error() {
    let input = vec![vec![5]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert_eq!(output, -4.0);
}

#[test]
fn proc2x2() {
    let input = vec![vec![3, 4], vec![4, 2]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.5920745920745918,
        epsilon = 0.000001
    ));
}

#[test]
fn proc2x2_error() {
    let input = vec![vec![3, 4], vec![4, -2]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert_eq!(output, -2.0);
}

#[test]
fn proc3x2() {
    let input = vec![vec![1000, 626, 782], vec![976, 814, 892]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, 0.0001679, epsilon = 0.000001));
}

#[test]
fn proc3x3() {
    let input = vec![vec![32, 10, 20], vec![20, 25, 18], vec![11, 17, 14]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.010967949934049852,
        epsilon = 0.000001
    ));
}

#[test]
fn proc3x3_unit() {
    let input = vec![vec![1, 0, 0], vec![0, 1, 0], vec![0, 0, 1]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert_eq!(output, 1.0);
}

#[test]
fn proc3x3_zero() {
    let input = vec![vec![0, 0, 0], vec![0, 0, 0], vec![0, 0, 0]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert_eq!(output, -3.0);
}

#[test]
fn proc3x4_large() {
    let input = vec![vec![11, 12, 18, 15], vec![15, 13, 13, 15], vec![15, 19, 19, 15]];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.8821660735808727,
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
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.01096124432190692,
        epsilon = 0.000001
    ));
}

#[test]
fn proc4x4_large() {
    let input = vec![
        vec![28, 28, 28, 0],
        vec![0, 0, 0, 16],
        vec![0, 0, 0, 5],
        vec![0, 0, 0, 7],
    ];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(f64, output, 0.0, epsilon = 0.000001));
}

#[test]
fn proc4x5_large() {
    let input = vec![
        vec![8, 3, 5, 5, 6],
        vec![4, 3, 8, 6, 5],
        vec![2, 5, 3, 7, 6],
        vec![4, 8, 2, 3, 6],
    ];
    let output = calculate(input, None).unwrap();
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
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.22200753799676337,
        epsilon = 0.000001
    ));
}

#[test]
fn proc5x5_small() {
    let input = vec![
        vec![1, 0, 0, 0, 0],
        vec![1, 1, 0, 1, 0],
        vec![1, 1, 0, 0, 1],
        vec![0, 0, 1, 2, 1],
        vec![1, 1, 2, 1, 1],
    ];
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.9712771262351103,
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
    let result = calculate(input, None).unwrap();
    dbg!(result);
    assert!(float_cmp::approx_eq!(
        f64,
        result,
        0.26314046636138944,
        epsilon = 0.000001
    ));
}

#[test]
fn proc9x7() {
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
    let output = calculate(input, None).unwrap();
    dbg!(output);
    assert!(float_cmp::approx_eq!(
        f64,
        output,
        0.11590654664515711,
        epsilon = 0.000001
    ));
}

#[test]
fn proc4x15_error() {
    let input = vec![
        vec![23, 22, 13, 22, 19, 16, 22, 22, 24, 20, 14, 16, 19, 16, 19],
        vec![26, 20, 6, 20, 13, 12, 21, 18, 19, 14, 14, 14, 18, 11, 14],
        vec![26, 22, 14, 22, 14, 17, 22, 21, 23, 23, 14, 18, 16, 12, 13],
        vec![26, 23, 13, 24, 18, 19, 24, 25, 22, 18, 18, 17, 21, 21, 18],
    ];
    let result = calculate(input, Some(200000000)).unwrap();
    dbg!(result);
    assert_eq!(result, -501.0);
}

#[test]
fn proc16x8_error() {
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
    let result = calculate(input, Some(200000000)).unwrap();
    dbg!(result);
    assert_eq!(result, -502.0);
}
