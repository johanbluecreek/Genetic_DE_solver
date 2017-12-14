/*********************
****** libmevac ******
**********************
Rust library to speed up fitness evaluation, and circumvent Julia bugs.
*********************/
extern crate meval;

use std::ffi::CStr;
use std::os::raw::c_char;

#[no_mangle]
pub extern fn evalpt(exprin: *const c_char, vars: &std::vec::Vec<*const c_char>, pnt: &std::vec::Vec<f64>, len: u32) -> f64 {

    // Create the expression
    let cexpr = unsafe { CStr::from_ptr(exprin) };
    let sexpr = cexpr.to_str().unwrap();
    let expr: meval::Expr = sexpr.parse().unwrap();

    // Create variables/Populate context
    let mut ctx = meval::Context::new();
    // `vars.len()` is not the length as in Julia... use `len`.
    for i in 0 .. len as usize {
        let cvar = unsafe { CStr::from_ptr(vars[i]) };
        let svar = cvar.to_str().unwrap();
        ctx.var(svar, pnt[i]);
    }

    // Evaluate expression in point
    let ret = expr.eval_with_context(ctx).unwrap();
    // Return
    return ret;
}
