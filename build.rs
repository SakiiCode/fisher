extern crate cc;

fn main() {
    //println!(r"cargo:rustc-link-lib=static=f2c");

    let mut build_common = cc::Build::new();
    build_common.file("src/asa643.c");

    if cfg!(windows) {
        build_common.flag("/fp:fast").compile("fexact");
    } else {
        build_common
            .flag("-ffast-math")
            .flag("-Wno-unused-result")
            .flag("-Wno-clobbered")
            .compile("fexact");
    }
}
