extern crate cc;

fn main() {
	println!(r"cargo:rustc-link-lib=static=f2c");
    

    cc::Build::new()
        .file("src/asa643.c")
        .compile("fexact");
}
