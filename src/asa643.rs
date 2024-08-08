#![allow(non_camel_case_types, deprecated)]
use libc::{c_int, int32_t};

type integer = int32_t;
type doublereal = libc::c_double;

extern "C" {
    pub fn fexact_(
        nrow: integer,
        ncol: integer,
        table: *mut doublereal,
        ldtabl: integer,
        expect: *mut doublereal,
        percnt: *mut doublereal,
        emin: *mut doublereal,
        prt: *mut doublereal,
        pre: *mut doublereal,
        ws: integer,
    ) -> c_int;
}
