use polars::frame::DataFrame;
use std::cmp::Ordering;
use crate::region::RegionCoordinates;
use crate::ubatch2::BSXBatch;

pub struct ReadFinalBatch {
    data: DataFrame,
    region: RegionCoordinates,
}

impl PartialEq for ReadFinalBatch {
    fn eq(&self, other: &Self) -> bool {
        self.region == other.region
    }
}

impl Eq for ReadFinalBatch {}

impl PartialOrd for ReadFinalBatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> { Some(self.cmp(other)) }
}

impl Ord for ReadFinalBatch {
    fn cmp(&self, other: &Self) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

impl BSXBatch for ReadFinalBatch {
    fn from_df(data_frame: DataFrame) -> Self {
        let coordinates = RegionCoordinates::new(
            {
                let first_idx = data_frame
                    .column("chr")
                    .unwrap()
                    .categorical()
                    .unwrap()
                    .physical()
                    .first()
                    .unwrap();
                data_frame
                    .column("chr")
                    .unwrap()
                    .categorical()
                    .unwrap()
                    .get_rev_map()
                    .get(first_idx)
                    .to_string()
            },
            data_frame
                .column("position")
                .unwrap()
                .u32()
                .unwrap()
                .first()
                .unwrap(),
            data_frame
                .column("position")
                .unwrap()
                .u32()
                .unwrap()
                .last()
                .unwrap(),
        );
        Self {
            data: data_frame,
            region: coordinates,
        }
    }

    fn get_data(&self) -> &DataFrame {
        &self.data
    }

    fn get_data_mut(&mut self) -> &mut DataFrame {
        &mut self.data
    }

    fn _set_data(&mut self, data_frame: DataFrame) {
        self.data = data_frame;
    }
}