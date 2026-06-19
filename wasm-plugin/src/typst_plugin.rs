use crate::{
    convert_to_mtl, convert_to_obj, convert_to_obj_bundle, convert_to_ply,
    convert_to_render_object_bundle, convert_to_stl, maquette_material_map, molecule_info,
};

#[link(wasm_import_module = "typst_env")]
extern "C" {
    fn wasm_minimal_protocol_write_args_to_buffer(ptr: *mut u8);
    fn wasm_minimal_protocol_send_result_to_host(ptr: *const u8, len: usize);
}

fn call2(a_len: usize, b_len: usize, f: fn(&[u8], &[u8]) -> Result<Vec<u8>, String>) -> i32 {
    let mut args = vec![0u8; a_len + b_len];
    unsafe {
        wasm_minimal_protocol_write_args_to_buffer(args.as_mut_ptr());
    }
    let (a, b) = args.split_at(a_len);
    send_result(f(a, b))
}

fn send_result(result: Result<Vec<u8>, String>) -> i32 {
    match result {
        Ok(bytes) => {
            unsafe {
                wasm_minimal_protocol_send_result_to_host(bytes.as_ptr(), bytes.len());
            }
            0
        }
        Err(err) => {
            let bytes = err.into_bytes();
            unsafe {
                wasm_minimal_protocol_send_result_to_host(bytes.as_ptr(), bytes.len());
            }
            1
        }
    }
}

fn call1(a_len: usize, f: fn(&[u8]) -> Result<Vec<u8>, String>) -> i32 {
    let mut args = vec![0u8; a_len];
    unsafe {
        wasm_minimal_protocol_write_args_to_buffer(args.as_mut_ptr());
    }
    send_result(f(&args))
}

#[no_mangle]
pub extern "C" fn to_obj(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_obj)
}

#[no_mangle]
pub extern "C" fn to_obj_bundle(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_obj_bundle)
}

#[no_mangle]
pub extern "C" fn render_object_bundle(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_render_object_bundle)
}

#[no_mangle]
pub extern "C" fn to_mtl(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_mtl)
}

#[no_mangle]
pub extern "C" fn material_map(obj_len: usize) -> i32 {
    call1(obj_len, maquette_material_map)
}

#[no_mangle]
pub extern "C" fn to_stl(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_stl)
}

#[no_mangle]
pub extern "C" fn to_ply(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, convert_to_ply)
}

#[no_mangle]
pub extern "C" fn info(data_len: usize, options_len: usize) -> i32 {
    call2(data_len, options_len, molecule_info)
}
