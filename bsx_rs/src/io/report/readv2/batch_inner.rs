use uuid::Uuid;
use std::cmp::Ordering;
use crate::io::report::readv2::{BatchStatus, ReadQueueData};
use crate::io::report::types::ReportType;

pub struct WorkBatch {
    data: ReadQueueData,
    status: BatchStatus,
    report_type: ReportType,
    uuid: Uuid,
}

impl Eq for WorkBatch {}

impl PartialEq<Self> for WorkBatch {
    fn eq(&self, other: &Self) -> bool {
        self.status == other.status
    }
}

impl Ord for WorkBatch {
    fn cmp(&self, other: &Self) -> Ordering {
        self.status.cmp(&other.status)
    }
}

impl PartialOrd<Self> for WorkBatch {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl WorkBatch {
    pub(crate) fn new(data: ReadQueueData, status: BatchStatus, report_type: ReportType) -> WorkBatch {
        let uuid = Uuid::new_v4();
        WorkBatch {
            data,
            status,
            report_type,
            uuid,
        }
    }
    pub(crate) fn into_data(self) -> ReadQueueData {
        self.data
    }

    pub fn data(&self) -> &ReadQueueData {
        &self.data
    }
    pub fn data_mut(&mut self) -> &mut ReadQueueData {
        &mut self.data
    }

    pub fn status(&self) -> BatchStatus {
        self.status
    }

    pub fn report_type(&self) -> ReportType {
        self.report_type
    }

    pub fn uuid(&self) -> Uuid {
        self.uuid
    }

    pub fn set_status(&mut self, status: BatchStatus) {
        self.status = status;
    }

    pub fn set_data(&mut self, data: ReadQueueData) {
        self.data = data;
    }
}