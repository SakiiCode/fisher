# Fisher RxC

Fast multithreaded implementation of calculating Fisher's exact test for any RxC size table. Written in Rust using [Maturin](https://github.com/PyO3/maturin).

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

Workspace size of 2e8 takes at most 800MB of RAM.

`fisher.sim(table, iterations)`

Calculate Fisher's exact test for 2D list by multithreaded Monte Carlo simulation. A modern CPU can quickly do 10^7 iterations and get accurate results.

`fisher.recursive(table)`

> [!WARNING]
> This is experimental and may not work on every CPU platform and generation. Any feedback is appreciated

Calculate Fisher's exact test by a multithreaded SIMD recursive algorithm. Despite extensive optimization efforts, this is still generally **much slower** than the _fisher.exact_ function. Only use for small tables with low numbers.

### Return values

`0 <= x <= 1`: p-value

`x < 0`: error code number, message printed to stdout

## Performance

See [benchmark.py](https://github.com/SakiiCode/fisher/blob/main/benchmark.py)

_AMD Ryzen 5600X running LMDE 7_

```
-- EXACT TEST --
fisher-rxc	0.2631	in 20.04s
rpy2		0.2631	in 20.64s
-----------
fisher-rxc	0.9981	in 4.71s
rpy2		0.9981	in 5.28s
-----------
-- MONTE-CARLO SIMULATION --
fisher-rxc	0.2631	in 0.58s
rpy2		0.2630	in 31.53s
-----------
fisher-rxc	0.9981	in 0.67s
rpy2		0.9981	in 36.75s
-----------
```

## Building

**Windows**

```sh
python3 -m venv .env
.env/Scripts/Activate.ps1
maturin develop --release
```

**Linux**

```sh
python3 -m venv .env
source .env/bin/activate
pip3 install maturin[patchelf]
maturin develop --release
# For benchmarks:
sudo apt install r-base
pip install rpy2
```

## References

Contingency table generator (ASA159): https://people.sc.fsu.edu/~jburkardt/c_src/asa159/asa159.html

Fisher' exact test - network algorithm (ASA643): https://netlib.org/toms/643.gz

Fisher's exact test - recursive: https://stackoverflow.com/questions/25368284/fishers-exact-test-for-bigger-than-2-by-2-contingency-table

Fortran to C transpiler: https://www.netlib.org/f2c/

C to Rust transpiler: https://github.com/immunant/c2rust
