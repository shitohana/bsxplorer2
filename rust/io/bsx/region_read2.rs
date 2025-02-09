use crate::io::bsx::read::BsxFileReader;
use bio_types::annot::contig::Contig;
use bio_types::annot::loc::Loc;
use bio_types::annot::pos::{Pos, SeqPosUnstranded};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::io::{Read, Seek};

#[derive(Eq, PartialEq, Serialize, Deserialize, Clone, Debug)]
struct BsxFileIndex {
    index: HashMap<String, BTreeMap<u32, usize>>,
}

impl BsxFileIndex {
    fn new() -> BsxFileIndex {
        let index = HashMap::new();
        Self { index }
    }

    fn from_reader<R: Read + Seek>(reader: BsxFileReader<R>) -> Self {
        let mut new = BsxFileIndex::new();
        for (idx, batch) in reader.enumerate() {
            let batch_pos = batch
                .expect("Failed to deserialize batch")
                .first_seq_pos(&mut None)
                .expect("Failed to get first seq pos");
            new.insert_node(idx, batch_pos);
        }
        new
    }

    fn insert_node<R: Display, S>(&mut self, idx: usize, position: Pos<R, S>) {
        let entry = self
            .index
            .entry(position.refid().to_string())
            .or_insert(BTreeMap::new());
        entry.insert(position.pos() as u32, idx);
    }

    fn query_contig<R, S>(&self, contig: Contig<R, S>) -> Option<Vec<usize>>
    where
        R: Display,
    {
        if let Some(btree) = self.index.get(&contig.refid().to_string()) {
            let start_pos = contig.start() as u32;
            let end_pos = start_pos + contig.length() as u32;
            let batches = btree
                .range(start_pos..end_pos)
                .map(|x| x.1)
                .cloned()
                .collect::<Vec<_>>();
            Some(batches)
        } else {
            None
        }
    }

    fn query_pos(&self, position: SeqPosUnstranded) -> Option<usize> {
        if let Some(btree) = self.index.get(&position.refid().to_string()) {
            let start_pos = position.start() as u32;
            btree.range(start_pos..).next().map(|x| (x.1.clone()))
        } else {
            None
        }
    }

    fn save_to_bin(&self, filename: &str) -> std::io::Result<()> {
        let encoded = bincode::serialize(self).expect("Failed to serialize");
        std::fs::write(filename, encoded)?;
        Ok(())
    }

    fn load_from_bin(filename: &str) -> std::io::Result<Self> {
        let data = std::fs::read(filename)?;
        Ok(bincode::deserialize(&data).expect("Failed to deserialize"))
    }

    fn save_to_json(&self, filename: &str) -> std::io::Result<()> {
        let json = serde_json::to_string(self).expect("Failed to serialize");
        std::fs::write(filename, json)?;
        Ok(())
    }

    fn load_from_json(filename: &str) -> std::io::Result<Self> {
        let data = std::fs::read_to_string(filename)?;
        Ok(serde_json::from_str(&data).expect("Failed to deserialize"))
    }
}
