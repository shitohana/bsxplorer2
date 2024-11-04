// use polars::prelude::*;
// use std::collections::{HashMap, HashSet};
// use std::hash::Hash;
// use std::io::{Read, Write};
// use polars::error::PolarsError::ColumnNotFound;
// 
// struct BsxWriter<K: Hash + Eq, W: Write + Send, R: Read + Send + Sync> {
//     writer: W,
//     reader: R,
//     hashmap: HashMap<K, u64>,
//     key_cols: Vec<String>,
//     key_bits: Vec<usize>,
//     key_values: Vec<HashSet<K>>
// }
// 
// impl<K, W, R> BsxWriter<K, W, R> {
//     fn new(writer: W, reader: R, key_cols: Vec<String>, key_bits: Vec<usize>) -> Self {
//         let hashmap: HashMap<K, u64> = HashMap::new();
//         let key_values: Vec<HashSet<K>> = key_cols.iter().map(|_| HashSet::new()).collect();
//         Self {writer, reader, hashmap, key_cols, key_bits, key_values}
//     }
//     
//     fn write_batch(&mut self, batch: &DataFrame) -> Result<(), PolarsError> {
//         // Check key columns exist
//         {
//             let batch_schema = batch.schema().iter_names().collect::<Vec<String>>();
//             for col in self.key_cols.iter() {
//                 if !batch_schema.contains(&col) { return Err(ColumnNotFound(col.into())) }
//             }
//         }
//         Ok(())
//     }
// }