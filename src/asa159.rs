#![allow(clippy::needless_return)]

/*
fn i4_max(i1: i32, i2: i32) -> i32 {
    if i2 < i1 {
        i1
    } else {
        i2
    }
}

fn i4_min(i1: i32, i2: i32) -> i32 {
    if i1 < i2 {
        i1
    } else {
        i2
    }
}

fn i4mat_print(mut m: i32, mut n: i32, a: &Vec<i32>, mut title: &str) {
    i4mat_print_some(m, n, a, 1, 1, m, n, title);
}

fn i4mat_print_some(
    mut m: i32,
    mut n: i32,
    mut a: &Vec<i32>,
    mut ilo: i32,
    mut jlo: i32,
    mut ihi: i32,
    mut jhi: i32,
    mut title: &str,
) {
    let mut i2hi: i32 = 0;
    let mut i2lo: i32 = 0;
    let mut j2hi: i32 = 0;
    println!("");
    println!("{}", title);
    if m <= 0 || n <= 0 {
        println!("");
        println!("  (None)");
        return;
    }
    for j2lo in (jlo..=jhi).step_by(10) {
        j2hi = j2lo + 10 - 1;
        j2hi = i4_min(j2hi, n);
        j2hi = i4_min(j2hi, jhi);
        println!("");
        print!("  Col:");
        for j in j2lo..=j2hi {
            print!("  {:>6}", j - 1);
        }
        println!("");
        println!("  Row");
        println!("");
        i2lo = i4_max(ilo, 1);
        i2hi = i4_min(ihi, m);
        for i in i2lo..=i2hi {
            print!("{:>5}:", i - 1);
            for j in j2lo..=j2hi {
                print!("  {:>6}", a[(i - 1 + (j - 1) * m)],);
            }
            println!("");
        }
    }
}

fn i4vec_print(mut n: i32, a: &Vec<i32>, mut title: &str) {
    println!("");
    println!("{}", title);
    println!("");
    for i in 0..n {
        println!("  {:>6}: {:>8}", i, a[i],);
    }
}*/

fn i4vec_sum(a: &Vec<u32>) -> u32 {
    a.iter().sum()
}

/*fn r8_uniform_01(seed: &mut i32) -> f64 {
    let k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if *seed < 0 {
        *seed += i32::MAX;
    }
    let r = *seed as f64 * 4.656612875E-10f64;
    dbg!(r);
    return r;
}*/

pub fn rcont2(
    nrow: usize,
    ncol: usize,
    nrowt: &Vec<u32>,
    ncolt: &Vec<u32>,
    key: &mut i32,
    fact: &Vec<f64>,
) -> Result<Vec<u32>, (i32, &'static str)> {
    let mut done1: i32;
    let mut done2: i32 = 0;
    let mut ia: u32;
    let mut iap: u32;
    let mut ib: u32 = 0;
    let mut ic: u32;
    let mut id: u32;
    let mut idp: u32;
    let mut ie: u32;
    let mut igp: u32;
    let mut ihp: u32;
    let mut ii: u32;
    let mut iip: u32;
    let mut j: u32;
    let mut jc: u32;
    let mut jwork: Vec<u32>;
    let mut lsm: i32;
    let mut lsp: i32;
    let mut nll: u32;
    let mut nlm: u32;
    let mut nlmp: u32;
    let mut nrowtl: u32;
    let mut ntotal: u32 = 0;
    let mut r: f64;
    let mut sumprb: f64;
    let mut x: f64;
    let mut y: f64;
    if *key == 0 {
        *key = 1;
        if nrow <= 1 {
            return Err((
                1,
                "RCONT - Fatal error!\nInput number of rows is less than 2.",
            ));
        }
        if ncol <= 1 {
            return Err((
                2,
                "RCONT - Fatal error!\n\tThe number of columns is less than 2.",
            ));
        }
        for i in 0..nrow {
            if nrowt[i] <= 0 {
                return Err((
                    3,
                    "RCONT - Fatal error!\n\tAn entry in the row sum vector is not positive.",
                ));
            }
        }
        for j in 0..ncol {
            if ncolt[j] <= 0 {
                return Err((
                    4,
                    "RCONT - Fatal error!\n\tAn entry in the column sum vector is not positive.",
                ));
            }
        }
        if i4vec_sum(ncolt) != i4vec_sum(nrowt) {
            return Err((
                5,
                "RCONT - Fatal error!\n\tThe row and column sum vectors do not have the same sum.",
            ));
        }
        ntotal = i4vec_sum(ncolt);
    }
    jwork = vec![0; ncol];
    for i in 0..(ncol - 1) {
        jwork[i] = ncolt[i];
    }
    jc = ntotal;
    let mut matrix = vec![0; nrow * ncol];
    let mut rng = fastrand::Rng::new();
    for l in 0..(nrow - 1) {
        nrowtl = nrowt[l];
        ia = nrowtl;
        ic = jc;
        jc -= nrowtl;
        for m in 0..(ncol - 1) {
            id = jwork[m];
            ie = ic;
            ic -= id;
            ib = ie - ia;
            ii = ib - id;
            if ie == 0 {
                ia = 0;
                for j in m..ncol {
                    matrix[l + j * nrow] = 0;
                }
                break;
            } else {
                r = rng.f64();
                done1 = 0;
                loop {
                    nlm = ((ia * id) as f64 / ie as f64 + 0.5f64) as u32;
                    iap = ia + 1;
                    idp = id + 1;
                    igp = idp - nlm;
                    ihp = iap - nlm;
                    nlmp = nlm + 1;
                    iip = ii + nlmp;
                    x = (fact[(iap - 1) as usize]
                        + fact[(ib) as usize]
                        + fact[(ic) as usize]
                        + fact[(idp - 1) as usize]
                        - fact[(ie) as usize]
                        - fact[(nlmp - 1) as usize]
                        - fact[(igp - 1) as usize]
                        - fact[(ihp - 1) as usize]
                        - fact[(iip - 1) as usize])
                        .exp();
                    if r <= x {
                        break;
                    }
                    sumprb = x;
                    y = x;
                    nll = nlm;
                    lsp = 0;
                    lsm = 0;
                    while lsp == 0 {
                        j = (id - nlm) * (ia - nlm);
                        if j == 0 {
                            lsp = 1;
                        } else {
                            nlm += 1;
                            x = x * j as f64 / (nlm * (ii + nlm)) as f64;
                            sumprb += x;
                            if r <= sumprb {
                                done1 = 1;
                                break;
                            }
                        }
                        done2 = 0;
                        while lsm == 0 {
                            j = nll * (ii + nll);
                            if j == 0 {
                                lsm = 1;
                                break;
                            } else {
                                nll -= 1;
                                y = y * j as f64 / ((id - nll) * (ia - nll)) as f64;
                                sumprb += y;
                                if r <= sumprb {
                                    nlm = nll;
                                    done2 = 1;
                                    break;
                                } else if lsp == 0 {
                                    break;
                                }
                            }
                        }
                        if done2 != 0 {
                            break;
                        }
                    }
                    if done1 != 0 {
                        break;
                    }
                    if done2 != 0 {
                        break;
                    }
                    r = rng.f64();
                    r *= sumprb;
                }
                matrix[l + m * nrow] = nlm;
                ia -= nlm;
                jwork[m] -= nlm;
            }
        }
        matrix[l + (ncol - 1) * nrow] = ia;
    }
    for j in 0..(ncol - 1) {
        matrix[nrow - 1 + j * nrow] = jwork[j];
    }
    matrix[nrow - 1 + (ncol - 1) * nrow] = ib - matrix[nrow - 1 + (ncol - 2) * nrow];
    return Ok(matrix);
}
