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
    fn malloc(_: libc::c_ulong) -> *mut libc::c_void;
    fn free(_: *mut libc::c_void);
    fn printf(_: *const libc::c_char, _: ...) -> libc::c_int;
    fn exp(_: libc::c_double) -> libc::c_double;
    fn log(_: libc::c_double) -> libc::c_double;
}
#[no_mangle]
pub unsafe extern "C" fn i4_max(i1: libc::c_int, i2: libc::c_int) -> libc::c_int {
    if i2 < i1 {
        i1
    } else {
        i2
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4_min(i1: libc::c_int, i2: libc::c_int) -> libc::c_int {
    if i1 < i2 {
        i1
    } else {
        i2
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4mat_print(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut a: *mut libc::c_int,
    mut title: *mut libc::c_char,
) {
    i4mat_print_some(m, n, a, 1 as libc::c_int, 1 as libc::c_int, m, n, title);
}
#[no_mangle]
pub unsafe extern "C" fn i4mat_print_some(
    mut m: libc::c_int,
    mut n: libc::c_int,
    mut a: *mut libc::c_int,
    mut ilo: libc::c_int,
    mut jlo: libc::c_int,
    mut ihi: libc::c_int,
    mut jhi: libc::c_int,
    mut title: *mut libc::c_char,
) {
    let mut i: libc::c_int = 0;
    let mut i2hi: libc::c_int = 0;
    let mut i2lo: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut j2hi: libc::c_int = 0;
    let mut j2lo: libc::c_int = 0;
    printf(b"\n\0" as *const u8 as *const libc::c_char);
    printf(b"%s\n\0" as *const u8 as *const libc::c_char, title);
    if m <= 0 as libc::c_int || n <= 0 as libc::c_int {
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        printf(b"  (None)\n\0" as *const u8 as *const libc::c_char);
        return;
    }
    for j2lo in (jlo..=jhi).step_by(10) {
        j2hi = j2lo + 10 as libc::c_int - 1 as libc::c_int;
        j2hi = i4_min(j2hi, n);
        j2hi = i4_min(j2hi, jhi);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        printf(b"  Col:\0" as *const u8 as *const libc::c_char);
        for j in j2lo..=j2hi {
            printf(
                b"  %6d\0" as *const u8 as *const libc::c_char,
                j - 1 as libc::c_int,
            );
        }
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        printf(b"  Row\n\0" as *const u8 as *const libc::c_char);
        printf(b"\n\0" as *const u8 as *const libc::c_char);
        i2lo = i4_max(ilo, 1 as libc::c_int);
        i2hi = i4_min(ihi, m);
        for i in i2lo..=i2hi {
            printf(
                b"%5d:\0" as *const u8 as *const libc::c_char,
                i - 1 as libc::c_int,
            );
            for j in j2lo..=j2hi {
                printf(
                    b"  %6d\0" as *const u8 as *const libc::c_char,
                    *a.offset((i - 1 as libc::c_int + (j - 1 as libc::c_int) * m) as isize),
                );
            }
            printf(b"\n\0" as *const u8 as *const libc::c_char);
        }
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4vec_print(
    mut n: libc::c_int,
    mut a: *mut libc::c_int,
    mut title: *mut libc::c_char,
) {
    let mut i: libc::c_int = 0;
    printf(b"\n\0" as *const u8 as *const libc::c_char);
    printf(b"%s\n\0" as *const u8 as *const libc::c_char, title);
    printf(b"\n\0" as *const u8 as *const libc::c_char);
    for i in 0..n {
        printf(
            b"  %6d: %8d\n\0" as *const u8 as *const libc::c_char,
            i,
            *a.offset(i as isize),
        );
    }
}
#[no_mangle]
pub unsafe extern "C" fn i4vec_sum(mut n: libc::c_int, mut a: *mut libc::c_int) -> libc::c_int {
    let mut sum: libc::c_int = 0;
    sum = 0 as libc::c_int;
    for i in 0..n {
        sum = sum + *a.offset(i as isize);
    }
    return sum;
}
#[no_mangle]
pub unsafe extern "C" fn r8_uniform_01(mut seed: *mut libc::c_int) -> libc::c_double {
    let mut k: libc::c_int = 0;
    let mut r: libc::c_double = 0.;
    k = *seed / 127773 as libc::c_int;
    *seed = 16807 as libc::c_int * (*seed - k * 127773 as libc::c_int) - k * 2836 as libc::c_int;
    if *seed < 0 as libc::c_int {
        *seed = *seed + 2147483647 as libc::c_int;
    }
    r = *seed as libc::c_double * 4.656612875E-10f64;
    return r;
}
#[no_mangle]
pub unsafe extern "C" fn rcont2(
    mut nrow: libc::c_int,
    mut ncol: libc::c_int,
    mut nrowt: *mut libc::c_int,
    mut ncolt: *mut libc::c_int,
    mut key: *mut libc::c_int,
    mut seed: *mut libc::c_int,
    mut matrix: *mut libc::c_int,
    mut ierror: *mut libc::c_int,
) {
    let mut done1: libc::c_int = 0;
    let mut done2: libc::c_int = 0;
    let mut fact: *mut libc::c_double = 0 as *const libc::c_double as *mut libc::c_double;
    let mut i: libc::c_int = 0;
    let mut ia: libc::c_int = 0;
    let mut iap: libc::c_int = 0;
    let mut ib: libc::c_int = 0;
    let mut ic: libc::c_int = 0;
    let mut id: libc::c_int = 0;
    let mut idp: libc::c_int = 0;
    let mut ie: libc::c_int = 0;
    let mut igp: libc::c_int = 0;
    let mut ihp: libc::c_int = 0;
    let mut ii: libc::c_int = 0;
    let mut iip: libc::c_int = 0;
    let mut j: libc::c_int = 0;
    let mut jc: libc::c_int = 0;
    let mut jwork: *mut libc::c_int = 0 as *mut libc::c_int;
    let mut l: libc::c_int = 0;
    let mut lsm: libc::c_int = 0;
    let mut lsp: libc::c_int = 0;
    let mut m: libc::c_int = 0;
    let mut nll: libc::c_int = 0;
    let mut nlm: libc::c_int = 0;
    let mut nlmp: libc::c_int = 0;
    let mut nrowtl: libc::c_int = 0;
    let mut ntotal: libc::c_int = 0 as libc::c_int;
    let mut r: libc::c_double = 0.;
    let mut sumprb: libc::c_double = 0.;
    let mut x: libc::c_double = 0.;
    let mut y: libc::c_double = 0.;
    *ierror = 0 as libc::c_int;
    if *key == 0 {
        *key = 1 as libc::c_int;
        if nrow <= 1 as libc::c_int {
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            printf(b"RCONT - Fatal error!\n\0" as *const u8 as *const libc::c_char);
            printf(
                b"  Input number of rows is less than 2.\n\0" as *const u8 as *const libc::c_char,
            );
            *ierror = 1 as libc::c_int;
            return;
        }
        if ncol <= 1 as libc::c_int {
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            printf(b"RCONT - Fatal error!\n\0" as *const u8 as *const libc::c_char);
            printf(
                b"  The number of columns is less than 2.\n\0" as *const u8 as *const libc::c_char,
            );
            *ierror = 2 as libc::c_int;
            return;
        }
        i = 0 as libc::c_int;
        while i < nrow {
            if *nrowt.offset(i as isize) <= 0 as libc::c_int {
                printf(b"\n\0" as *const u8 as *const libc::c_char);
                printf(b"RCONT - Fatal error!\n\0" as *const u8 as *const libc::c_char);
                printf(
                    b"  An entry in the row sum vector is not positive.\n\0" as *const u8
                        as *const libc::c_char,
                );
                *ierror = 3 as libc::c_int;
                return;
            }
            i += 1;
        }
        j = 0 as libc::c_int;
        while j < ncol {
            if *ncolt.offset(j as isize) <= 0 as libc::c_int {
                printf(b"\n\0" as *const u8 as *const libc::c_char);
                printf(b"RCONT - Fatal error!\n\0" as *const u8 as *const libc::c_char);
                printf(
                    b"  An entry in the column sum vector is not positive.\n\0" as *const u8
                        as *const libc::c_char,
                );
                *ierror = 4 as libc::c_int;
                return;
            }
            j += 1;
        }
        if i4vec_sum(ncol, ncolt) != i4vec_sum(nrow, nrowt) {
            printf(b"\n\0" as *const u8 as *const libc::c_char);
            printf(b"RCONT - Fatal error!\n\0" as *const u8 as *const libc::c_char);
            printf(
                b"  The row and column sum vectors do not have the same sum.\n\0" as *const u8
                    as *const libc::c_char,
            );
            *ierror = 6 as libc::c_int;
            return;
        }
        ntotal = i4vec_sum(ncol, ncolt);
        if !fact.is_null() {
            free(fact as *mut libc::c_void);
        }
        fact = malloc(
            ((ntotal + 1 as libc::c_int) as libc::c_ulong)
                .wrapping_mul(::core::mem::size_of::<libc::c_double>() as libc::c_ulong),
        ) as *mut libc::c_double;
        x = 0.0f64;
        *fact.offset(0 as libc::c_int as isize) = 0.0f64;
        i = 1 as libc::c_int;
        while i <= ntotal {
            x = x + log(i as libc::c_double);
            *fact.offset(i as isize) = x;
            i += 1;
        }
    }
    jwork = malloc(
        (ncol as libc::c_ulong)
            .wrapping_mul(::core::mem::size_of::<libc::c_int>() as libc::c_ulong),
    ) as *mut libc::c_int;
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
            if ie == 0 as libc::c_int {
                ia = 0 as libc::c_int;
                for j in m..ncol {
                    *matrix.offset((l + j * nrow) as isize) = 0 as libc::c_int;
                }
                break;
            } else {
                r = r8_uniform_01(seed);
                done1 = 0 as libc::c_int;
                loop {
                    nlm = ((ia * id) as libc::c_double / ie as libc::c_double + 0.5f64)
                        as libc::c_int;
                    iap = ia + 1 as libc::c_int;
                    idp = id + 1 as libc::c_int;
                    igp = idp - nlm;
                    ihp = iap - nlm;
                    nlmp = nlm + 1 as libc::c_int;
                    iip = ii + nlmp;
                    x = exp(*fact.offset((iap - 1 as libc::c_int) as isize)
                        + *fact.offset(ib as isize)
                        + *fact.offset(ic as isize)
                        + *fact.offset((idp - 1 as libc::c_int) as isize)
                        - *fact.offset(ie as isize)
                        - *fact.offset((nlmp - 1 as libc::c_int) as isize)
                        - *fact.offset((igp - 1 as libc::c_int) as isize)
                        - *fact.offset((ihp - 1 as libc::c_int) as isize)
                        - *fact.offset((iip - 1 as libc::c_int) as isize));
                    if r <= x {
                        break;
                    }
                    sumprb = x;
                    y = x;
                    nll = nlm;
                    lsp = 0 as libc::c_int;
                    lsm = 0 as libc::c_int;
                    while lsp == 0 {
                        j = (id - nlm) * (ia - nlm);
                        if j == 0 as libc::c_int {
                            lsp = 1 as libc::c_int;
                        } else {
                            nlm = nlm + 1 as libc::c_int;
                            x = x * j as libc::c_double / (nlm * (ii + nlm)) as libc::c_double;
                            sumprb = sumprb + x;
                            if r <= sumprb {
                                done1 = 1 as libc::c_int;
                                break;
                            }
                        }
                        done2 = 0 as libc::c_int;
                        while lsm == 0 {
                            j = nll * (ii + nll);
                            if j == 0 as libc::c_int {
                                lsm = 1 as libc::c_int;
                                break;
                            } else {
                                nll = nll - 1 as libc::c_int;
                                y = y * j as libc::c_double
                                    / ((id - nll) * (ia - nll)) as libc::c_double;
                                sumprb = sumprb + y;
                                if r <= sumprb {
                                    nlm = nll;
                                    done2 = 1 as libc::c_int;
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
        *matrix.offset((l + (ncol - 1 as libc::c_int) * nrow) as isize) = ia;
    }
    for j in 0..(ncol - 1) {
        *matrix.offset((nrow - 1 as libc::c_int + j * nrow) as isize) = *jwork.offset(j as isize);
    }
    *matrix.offset((nrow - 1 as libc::c_int + (ncol - 1 as libc::c_int) * nrow) as isize) =
        ib - *matrix.offset((nrow - 1 as libc::c_int + (ncol - 2 as libc::c_int) * nrow) as isize);
    free(jwork as *mut libc::c_void);
}
