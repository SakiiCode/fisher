#![allow(
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut
)]
extern "C" {
    fn malloc(_: u64) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
}
#[no_mangle]
pub unsafe extern "C" fn i4_max(i1: i32, i2: i32) -> i32 {
    if i2 < i1 {
        i1
    } else {
        i2
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4_min(i1: i32, i2: i32) -> i32 {
    if i1 < i2 {
        i1
    } else {
        i2
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4mat_print(mut m: i32, mut n: i32, mut a: *mut i32, mut title: &str) {
    i4mat_print_some(m, n, a, 1, 1, m, n, title);
}
#[no_mangle]
pub unsafe extern "C" fn i4mat_print_some(
    mut m: i32,
    mut n: i32,
    mut a: *mut i32,
    mut ilo: i32,
    mut jlo: i32,
    mut ihi: i32,
    mut jhi: i32,
    mut title: &str,
) {
    let mut i: i32 = 0;
    let mut i2hi: i32 = 0;
    let mut i2lo: i32 = 0;
    let mut j: i32 = 0;
    let mut j2hi: i32 = 0;
    let mut j2lo: i32 = 0;
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
                print!("  {:>6}", *a.offset((i - 1 + (j - 1) * m) as isize),);
            }
            println!("");
        }
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4vec_print(mut n: i32, mut a: *mut i32, mut title: &str) {
    let mut i: i32 = 0;
    println!("");
    println!("{}", title);
    println!("");
    for i in 0..n {
        println!("  {:>6}: {:>8}", i, *a.offset(i as isize),);
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4vec_sum(mut n: i32, mut a: *mut i32) -> i32 {
    let mut sum: i32 = 0;
    sum = 0;
    for i in 0..n {
        sum = sum + *a.offset(i as isize);
    }
    return sum;
}
#[no_mangle]
pub unsafe extern "C" fn r8_uniform_01(mut seed: *mut i32) -> f64 {
    let mut k: i32 = 0;
    let mut r: f64 = 0.;
    k = *seed / 127773;
    *seed = 16807 * (*seed - k * 127773) - k * 2836;
    if *seed < 0 {
        *seed = *seed + 2147483647;
    }
    r = *seed as f64 * 4.656612875E-10f64;
    return r;
}
#[no_mangle]
pub unsafe extern "C" fn rcont2(
    mut nrow: i32,
    mut ncol: i32,
    mut nrowt: *mut i32,
    mut ncolt: *mut i32,
    mut key: *mut i32,
    mut seed: *mut i32,
    mut matrix: *mut i32,
    mut ierror: *mut i32,
) {
    let mut done1: i32 = 0;
    let mut done2: i32 = 0;
    let mut fact: *mut f64 = 0 as *const f64 as *mut f64;
    let mut i: i32 = 0;
    let mut ia: i32 = 0;
    let mut iap: i32 = 0;
    let mut ib: i32 = 0;
    let mut ic: i32 = 0;
    let mut id: i32 = 0;
    let mut idp: i32 = 0;
    let mut ie: i32 = 0;
    let mut igp: i32 = 0;
    let mut ihp: i32 = 0;
    let mut ii: i32 = 0;
    let mut iip: i32 = 0;
    let mut j: i32 = 0;
    let mut jc: i32 = 0;
    let mut jwork: *mut i32 = 0 as *mut i32;
    let mut l: i32 = 0;
    let mut lsm: i32 = 0;
    let mut lsp: i32 = 0;
    let mut m: i32 = 0;
    let mut nll: i32 = 0;
    let mut nlm: i32 = 0;
    let mut nlmp: i32 = 0;
    let mut nrowtl: i32 = 0;
    let mut ntotal: i32 = 0;
    let mut r: f64 = 0.;
    let mut sumprb: f64 = 0.;
    let mut x: f64 = 0.;
    let mut y: f64 = 0.;
    *ierror = 0;
    if *key == 0 {
        *key = 1;
        if nrow <= 1 {
            println!("");
            println!("RCONT - Fatal error!");
            println!("  Input number of rows is less than 2.");
            *ierror = 1;
            return;
        }
        if ncol <= 1 {
            println!("");
            println!("RCONT - Fatal error!");
            println!("  The number of columns is less than 2.");
            *ierror = 2;
            return;
        }
        for i in 0..nrow {
            if *nrowt.offset(i as isize) <= 0 {
                println!("");
                println!("RCONT - Fatal error!");
                println!("  An entry in the row sum vector is not positive.");
                *ierror = 3;
                return;
            }
        }
        for j in 0..ncol {
            if *ncolt.offset(j as isize) <= 0 {
                println!("");
                println!("RCONT - Fatal error!");
                println!("  An entry in the column sum vector is not positive.");
                *ierror = 4;
                return;
            }
        }
        if i4vec_sum(ncol, ncolt) != i4vec_sum(nrow, nrowt) {
            println!("");
            println!("RCONT - Fatal error!");
            println!("  The row and column sum vectors do not have the same sum.");
            *ierror = 6;
            return;
        }
        ntotal = i4vec_sum(ncol, ncolt);
        if !fact.is_null() {
            free(fact as *mut libc::c_void);
        }
        fact = malloc(((ntotal + 1) as u64).wrapping_mul(::core::mem::size_of::<f64>() as u64))
            as *mut f64;
        x = 0.0f64;
        *fact.offset(0 as isize) = 0.0f64;
        for i in 1..=ntotal {
            x = x + (i as f64).ln();
            *fact.offset(i as isize) = x;
        }
    }
    jwork = malloc((ncol as u64).wrapping_mul(::core::mem::size_of::<i32>() as u64)) as *mut i32;
    for i in 0..(ncol - 1) {
        *jwork.offset(i as isize) = *ncolt.offset(i as isize);
    }
    jc = ntotal;
    for l in 0..(nrow - 1) {
        nrowtl = *nrowt.offset(l as isize);
        ia = nrowtl;
        ic = jc;
        jc = jc - nrowtl;
        for m in 0..(ncol - 1) {
            id = *jwork.offset(m as isize);
            ie = ic;
            ic = ic - id;
            ib = ie - ia;
            ii = ib - id;
            if ie == 0 {
                ia = 0;
                for j in m..ncol {
                    *matrix.offset((l + j * nrow) as isize) = 0;
                }
                break;
            } else {
                r = r8_uniform_01(seed);
                done1 = 0;
                loop {
                    nlm = ((ia * id) as f64 / ie as f64 + 0.5f64) as i32;
                    iap = ia + 1;
                    idp = id + 1;
                    igp = idp - nlm;
                    ihp = iap - nlm;
                    nlmp = nlm + 1;
                    iip = ii + nlmp;
                    x = (*fact.offset((iap - 1) as isize)
                        + *fact.offset(ib as isize)
                        + *fact.offset(ic as isize)
                        + *fact.offset((idp - 1) as isize)
                        - *fact.offset(ie as isize)
                        - *fact.offset((nlmp - 1) as isize)
                        - *fact.offset((igp - 1) as isize)
                        - *fact.offset((ihp - 1) as isize)
                        - *fact.offset((iip - 1) as isize))
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
                            nlm = nlm + 1;
                            x = x * j as f64 / (nlm * (ii + nlm)) as f64;
                            sumprb = sumprb + x;
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
                                nll = nll - 1;
                                y = y * j as f64 / ((id - nll) * (ia - nll)) as f64;
                                sumprb = sumprb + y;
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
                    r = r8_uniform_01(seed);
                    r = sumprb * r;
                }
                *matrix.offset((l + m * nrow) as isize) = nlm;
                ia = ia - nlm;
                *jwork.offset(m as isize) = *jwork.offset(m as isize) - nlm;
            }
        }
        *matrix.offset((l + (ncol - 1) * nrow) as isize) = ia;
    }
    for j in 0..(ncol - 1) {
        *matrix.offset((nrow - 1 + j * nrow) as isize) = *jwork.offset(j as isize);
    }
    *matrix.offset((nrow - 1 + (ncol - 1) * nrow) as isize) =
        ib - *matrix.offset((nrow - 1 + (ncol - 2) * nrow) as isize);
    free(jwork as *mut libc::c_void);
}
