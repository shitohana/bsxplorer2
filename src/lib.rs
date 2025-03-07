#![feature(btree_cursors)]
#![feature(ascii_char)]
#![feature(vec_pop_if)]
#![feature(mpmc_channel)]
#![feature(assert_matches)]
#![feature(path_add_extension)]
#![allow(unused_assignments)]
extern crate core;

pub mod data_structs;
pub mod exports;
pub mod io;
#[cfg(feature = "plots")]
pub mod plots;
pub mod tools;
pub mod utils;

#[cfg(feature = "python")]
use crate::io::bsx;
#[cfg(feature = "python")]
use crate::io::report::schema;
#[cfg(feature = "python")]
use pyo3::prelude::*;
