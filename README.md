# Fisher RxC

Fast multithreaded implementation for calculating Fisher's exact test for any RxC size tables written in Rust using [Maturin](https://github.com/PyO3/maturin).

## Installation

```bash
pip install fisher-rxc
```

```python
import fisher
```

## Usage

`fisher.exact(table, workspace=None)`

Calculate Fisher's exact test for 2D list according to Mehta & Patel's Network Algorithm. If workspace size is not provided, it will be "guessed" dynamically.

Workspace size of 2e8 takes ~800MB of RAM.

`fisher.sim(table, iterations)`

Calculate Fisher's exact test for 2D list by Monte Carlo simulation. This multithreaded implementation can be more than 100x faster than using R through _rpy2_.

A modern CPU can quickly do 10^9 iterations and get very accurate results.

`fisher.recursive(table)`

Calculate Fisher's exact test by a simple multithreaded recursive algorithm. This is generally **much slower** than the _fisher.exact_ function. Only use for small tables with low numbers.

### Return values

`0 <= x <= 1`: p-value

`x < 0`: error code number, message printed to stdout

## References

Contingency table generator (ASA159): https://people.sc.fsu.edu/~jburkardt/c_src/asa159/asa159.html

Fisher' exact test - network algorithm (ASA643): https://netlib.org/toms/643.gz

Fisher's exact test - recursive: https://stackoverflow.com/questions/25368284/fishers-exact-test-for-bigger-than-2-by-2-contingency-table

Fortran to C transpiler: https://www.netlib.org/toms/

C to Rust transpiler: https://github.com/immunant/c2rust
