use std::ffi::{c_double, c_int};



extern "C" {
    pub fn fexact_(
        nrow: c_int,
        ncol: c_int,
        table: *mut c_double,
        ldtabl: c_int,
        expect: *mut c_double,
        percnt: *mut c_double,
        emin: *mut c_double,
        prt: *mut c_double,
        pre: *mut c_double,
        ws: c_int,
    ) -> c_int;
}
