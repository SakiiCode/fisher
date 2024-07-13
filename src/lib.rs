use pyo3::prelude::*;

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
