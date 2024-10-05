#![allow(incomplete_features)]
#![feature(generic_const_exprs)]
#![feature(portable_simd)]
#![allow(clippy::needless_return)]
#![allow(clippy::ptr_arg)]

use core::f64;
use pyo3::prelude::*;

mod asa159;
mod asa643;
mod fixedsize;
mod math;
mod network;
mod simulation;

#[pyfunction]
pub fn recursive(table: Vec<Vec<i32>>) -> PyResult<f64> {
    Ok(fixedsize::calculate(table)?)
}

#[pyfunction]
pub fn sim(table: Vec<Vec<i32>>, iterations: i32) -> PyResult<f64> {
    Ok(simulation::calculate(table, iterations)?)
}

#[pyfunction]
#[pyo3(signature = (table, workspace=None))]
pub fn exact(table: Vec<Vec<i32>>, workspace: Option<i32>) -> PyResult<f64> {
    Ok(network::calculate(table, workspace)?)
}

/// A Python module implemented in Rust.
#[pymodule]
fn fisher(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(recursive, m)?)?;
    m.add_function(wrap_pyfunction!(sim, m)?)?;
    m.add_function(wrap_pyfunction!(exact, m)?)?;
    Ok(())
}
