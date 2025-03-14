/*******************************************************************************
 Copyright (c) 2025
 The Prosperity Public License 3.0.0

 Contributor: [shitohana](https://github.com/shitohana)

 Source Code: https://github.com/shitohana/BSXplorer
 ******************************************************************************/

/// ****************************************************************************
/// * Copyright (c) 2025
/// ***************************************************************************

#![feature(btree_cursors)]
#![feature(ascii_char)]
#![feature(vec_pop_if)]
#![feature(mpmc_channel)]
#![feature(assert_matches)]
#![feature(path_add_extension)]
#![feature(portable_simd)]
#![allow(unused_assignments)]
extern crate core;

pub mod data_structs;
pub mod exports;
pub mod io;
#[cfg(feature = "plots")]
pub mod plots;
pub mod tools;
pub mod utils;
