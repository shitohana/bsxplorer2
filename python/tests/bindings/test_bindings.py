from src.bsx2 import _bsx2 as x


def test_module():
    assert hasattr(x, "__all__")


def test_classes():
    for item in [
        "Strand",
        "Context",
        "BsxColumns",
        "AggMethod",
        "BsxBatch",
        "ContextData",
        "ReportTypeSchema",
        "LazyBsxBatch",
        "GenomicPosition",
        "Contig",
        "GffEntryAttributes",
        "GffEntry",
        "HcAnnotStore",
        "HcAnnotStoreIterator"
    ]:
        assert hasattr(x, item)


def test_strand_methods():
    strand_class = x.Strand
    assert hasattr(strand_class, 'Forward')
    assert hasattr(strand_class, 'Reverse')
    assert hasattr(strand_class, 'Null')
    # Check instance methods on enum values
    strand_instance = strand_class.Forward
    assert hasattr(strand_instance, 'name')


def test_context_methods():
    context_class = x.Context
    assert hasattr(context_class, 'CG')
    assert hasattr(context_class, 'CHG')
    assert hasattr(context_class, 'CHH')
    # Check instance methods on enum values
    context_instance = context_class.CG
    assert hasattr(context_instance, 'name')


def test_bsx_columns_methods():
    bsx_columns_class = x.BsxColumns
    assert hasattr(bsx_columns_class, 'Chr')
    assert hasattr(bsx_columns_class, 'Position')
    assert hasattr(bsx_columns_class, 'Strand')
    assert hasattr(bsx_columns_class, 'Context')
    assert hasattr(bsx_columns_class, 'CountM')
    assert hasattr(bsx_columns_class, 'CountTotal')
    assert hasattr(bsx_columns_class, 'Density')
    assert hasattr(bsx_columns_class, 'schema')
    assert hasattr(bsx_columns_class, 'colnames')
    assert hasattr(bsx_columns_class, 'has_name')
    # Check instance methods
    bsx_columns_instance = bsx_columns_class.Chr
    assert hasattr(bsx_columns_instance, 'as_str')
    assert hasattr(bsx_columns_instance, '__str__')
    assert hasattr(bsx_columns_instance, '__repr__')


def test_agg_method_methods():
    agg_method_class = x.AggMethod
    assert hasattr(agg_method_class, 'Mean')
    assert hasattr(agg_method_class, 'GeometricMean')
    assert hasattr(agg_method_class, 'Median')
    assert hasattr(agg_method_class, 'Max')
    assert hasattr(agg_method_class, 'Min')


def test_bsx_batch_methods():
    bsx_batch_class = x.BsxBatch
    assert hasattr(bsx_batch_class, '__init__')
    assert hasattr(bsx_batch_class, 'from_dataframe')
    assert hasattr(bsx_batch_class, 'schema')
    assert hasattr(bsx_batch_class, 'empty')
    assert hasattr(bsx_batch_class, 'concat')
    # Create an empty instance to test instance methods
    bsx_batch_instance = bsx_batch_class.empty()
    assert hasattr(bsx_batch_instance, 'chr')
    assert hasattr(bsx_batch_instance, 'position')
    assert hasattr(bsx_batch_instance, 'strand')
    assert hasattr(bsx_batch_instance, 'context')
    assert hasattr(bsx_batch_instance, 'count_m')
    assert hasattr(bsx_batch_instance, 'count_total')
    assert hasattr(bsx_batch_instance, 'density')
    assert hasattr(bsx_batch_instance, 'column')
    assert hasattr(bsx_batch_instance, 'is_empty')
    assert hasattr(bsx_batch_instance, 'split_at')
    assert hasattr(bsx_batch_instance, 'rechunk')
    assert hasattr(bsx_batch_instance, 'data')
    assert hasattr(bsx_batch_instance, 'into_dataframe')
    assert hasattr(bsx_batch_instance, 'slice')
    assert hasattr(bsx_batch_instance, 'add_context_data')
    assert hasattr(bsx_batch_instance, 'extend')
    assert hasattr(bsx_batch_instance, 'extend_unchecked')
    assert hasattr(bsx_batch_instance, 'discretise')
    assert hasattr(bsx_batch_instance, 'partition')
    assert hasattr(bsx_batch_instance, 'seqname')
    assert hasattr(bsx_batch_instance, 'first_pos')
    assert hasattr(bsx_batch_instance, 'last_pos')
    assert hasattr(bsx_batch_instance, 'first_genomic_pos')
    assert hasattr(bsx_batch_instance, 'last_genomic_pos')
    assert hasattr(bsx_batch_instance, 'as_contig')
    assert hasattr(bsx_batch_instance, 'as_binom')
    assert hasattr(bsx_batch_instance, 'into_report')
    assert hasattr(bsx_batch_instance, 'lazy')
    assert hasattr(bsx_batch_instance, 'height')
    assert hasattr(bsx_batch_instance, '__len__')
    assert hasattr(bsx_batch_instance, '__repr__')


def test_context_data_methods():
    context_data_class = x.ContextData
    assert hasattr(context_data_class, '__init__')
    assert hasattr(context_data_class, 'empty')
    # Create an empty instance to test instance methods
    context_data_instance = context_data_class.empty()
    assert hasattr(context_data_instance, 'is_empty')
    assert hasattr(context_data_instance, 'take')
    assert hasattr(context_data_instance, 'to_decoded_df')
    assert hasattr(context_data_instance, '__len__')


def test_report_type_schema_methods():
    report_type_schema_class = x.ReportTypeSchema
    assert hasattr(report_type_schema_class, 'Bismark')
    assert hasattr(report_type_schema_class, 'CgMap')
    assert hasattr(report_type_schema_class, 'BedGraph')
    assert hasattr(report_type_schema_class, 'Coverage')
    # Check instance methods
    report_type_instance = report_type_schema_class.Bismark
    assert hasattr(report_type_instance, 'col_names')
    assert hasattr(report_type_instance, 'schema')
    assert hasattr(report_type_instance, 'chr_col')
    assert hasattr(report_type_instance, 'position_col')
    assert hasattr(report_type_instance, 'context_col')
    assert hasattr(report_type_instance, 'strand_col')
    assert hasattr(report_type_instance, 'need_align')


def test_lazy_bsx_batch_methods():
    lazy_bsx_batch_class = x.LazyBsxBatch
    assert hasattr(lazy_bsx_batch_class, '__init__')
    # Create an instance to test instance methods
    bsx_batch = x.BsxBatch.empty()
    lazy_bsx_batch_instance = lazy_bsx_batch_class(bsx_batch)
    assert hasattr(lazy_bsx_batch_instance, 'collect')
    assert hasattr(lazy_bsx_batch_instance, 'filter_pos_lt')
    assert hasattr(lazy_bsx_batch_instance, 'filter_pos_gt')
    assert hasattr(lazy_bsx_batch_instance, 'filter_coverage_gt')
    assert hasattr(lazy_bsx_batch_instance, 'filter_strand')
    assert hasattr(lazy_bsx_batch_instance, 'filter_context')


def test_genomic_position_methods():
    genomic_position_class = x.GenomicPosition
    assert hasattr(genomic_position_class, '__init__')
    # Create an instance to test instance methods
    genomic_position_instance = genomic_position_class("chr1", 100)
    assert hasattr(genomic_position_instance, '__str__')
    assert hasattr(genomic_position_instance, '__repr__')
    assert hasattr(genomic_position_instance, '__add__')
    assert hasattr(genomic_position_instance, '__sub__')
    assert hasattr(genomic_position_instance, 'is_zero')
    assert hasattr(genomic_position_instance, 'seqname')
    assert hasattr(genomic_position_instance, 'position')


def test_contig_methods():
    contig_class = x.Contig
    assert hasattr(contig_class, '__init__')
    # Create an instance to test instance methods
    contig_instance = contig_class("chr1", 100, 200, "+")
    assert hasattr(contig_instance, 'seqname')
    assert hasattr(contig_instance, 'start')
    assert hasattr(contig_instance, 'end')
    assert hasattr(contig_instance, 'strand')
    assert hasattr(contig_instance, 'strand_str')
    assert hasattr(contig_instance, 'length')
    assert hasattr(contig_instance, 'start_gpos')
    assert hasattr(contig_instance, 'end_gpos')
    assert hasattr(contig_instance, 'extend_upstream')
    assert hasattr(contig_instance, 'extend_downstream')
    assert hasattr(contig_instance, 'is_in')
    assert hasattr(contig_instance, '__str__')
    assert hasattr(contig_instance, '__repr__')


def test_gff_entry_attributes_methods():
    gff_entry_attributes_class = x.GffEntryAttributes
    assert hasattr(gff_entry_attributes_class, '__init__')
    # Create an instance to test instance methods
    gff_entry_attributes_instance = gff_entry_attributes_class()
    assert hasattr(gff_entry_attributes_instance, 'id')
    assert hasattr(gff_entry_attributes_instance, 'name')
    assert hasattr(gff_entry_attributes_instance, 'alias')
    assert hasattr(gff_entry_attributes_instance, 'parent')
    assert hasattr(gff_entry_attributes_instance, 'other')
    assert hasattr(gff_entry_attributes_instance, '__repr__')


def test_gff_entry_methods():
    gff_entry_class = x.GffEntry
    assert hasattr(gff_entry_class, '__init__')
    # Create an instance to test instance methods
    contig = x.Contig("chr1", 100, 200, "+")
    gff_entry_instance = gff_entry_class(contig, "test", "gene", 1.0, 0, "test_id")
    assert hasattr(gff_entry_instance, 'id')
    assert hasattr(gff_entry_instance, 'contig')
    assert hasattr(gff_entry_instance, 'source')
    assert hasattr(gff_entry_instance, 'feature_type')
    assert hasattr(gff_entry_instance, 'score')
    assert hasattr(gff_entry_instance, 'phase')
    assert hasattr(gff_entry_instance, 'attributes')
    assert hasattr(gff_entry_instance, '__repr__')


def test_hc_annot_store_methods():
    hc_annot_store_class = x.HcAnnotStore
    assert hasattr(hc_annot_store_class, '__init__')
    assert hasattr(hc_annot_store_class, 'from_gff')
    assert hasattr(hc_annot_store_class, 'from_bed')
    # Create an instance to test instance methods
    hc_annot_store_instance = hc_annot_store_class()
    assert hasattr(hc_annot_store_instance, 'len')
    assert hasattr(hc_annot_store_instance, 'is_empty')
    assert hasattr(hc_annot_store_instance, 'get_entry')
    assert hasattr(hc_annot_store_instance, 'get_entries_regex')
    assert hasattr(hc_annot_store_instance, 'get_children')
    assert hasattr(hc_annot_store_instance, 'get_parent')
    assert hasattr(hc_annot_store_instance, 'genomic_query')
    assert hasattr(hc_annot_store_instance, 'get_feature_types')
    assert hasattr(hc_annot_store_instance, 'add_flanks')
    assert hasattr(hc_annot_store_instance, 'iter')
    assert hasattr(hc_annot_store_instance, '__len__')
    assert hasattr(hc_annot_store_instance, '__repr__')
    assert hasattr(hc_annot_store_instance, '__iter__')
    assert hasattr(hc_annot_store_instance, 'init_tree')
    assert hasattr(hc_annot_store_instance, 'init_imap')


def test_hc_annot_store_iterator_methods():
    # Create an iterator instance through HcAnnotStore
    hc_annot_store = x.HcAnnotStore()
    hc_annot_store_iterator = hc_annot_store.iter()
    assert hasattr(hc_annot_store_iterator, '__iter__')
    assert hasattr(hc_annot_store_iterator, '__next__')

def test_bsx_file_reader_methods():
    bsx_file_reader_class = x.BsxFileReader
    assert hasattr(bsx_file_reader_class, '__init__')
    # Note: Can't easily create instance without valid file for testing instance methods
    # but we can test that the class has the expected methods defined


def test_batch_index_methods():
    batch_index_class = x.BatchIndex
    assert hasattr(batch_index_class, '__init__')
    assert hasattr(batch_index_class, 'load')
    # Create an instance to test instance methods
    batch_index_instance = batch_index_class()
    assert hasattr(batch_index_instance, 'insert')
    assert hasattr(batch_index_instance, 'find')
    assert hasattr(batch_index_instance, 'get_chr_order')
    assert hasattr(batch_index_instance, 'get_chr_index')
    assert hasattr(batch_index_instance, 'save')
    assert hasattr(batch_index_instance, 'sort')


def test_ipc_compression_methods():
    ipc_compression_class = x.IpcCompression
    assert hasattr(ipc_compression_class, 'LZ4')
    assert hasattr(ipc_compression_class, 'ZSTD')


def test_bsx_file_writer_methods():
    bsx_file_writer_class = x.BsxFileWriter
    assert hasattr(bsx_file_writer_class, '__init__')
    assert hasattr(bsx_file_writer_class, 'from_sink_and_fai')
    assert hasattr(bsx_file_writer_class, 'from_sink_and_fasta')
    # Note: Can't easily create instance without valid parameters for testing instance methods


def test_compression_methods():
    compression_class = x.Compression
    assert hasattr(compression_class, 'No')
    assert hasattr(compression_class, 'Gz')
    assert hasattr(compression_class, 'Zstd')
    assert hasattr(compression_class, 'Lz4')
    assert hasattr(compression_class, 'Xz2')
    assert hasattr(compression_class, 'Bzip2')
    assert hasattr(compression_class, 'Zip')


def test_report_reader_methods():
    report_reader_class = x.ReportReader
    assert hasattr(report_reader_class, '__init__')
    # Note: Can't easily create instance without valid parameters for testing instance methods


def test_report_writer_methods():
    report_writer_class = x.ReportWriter
    assert hasattr(report_writer_class, '__init__')
    # Note: Can't easily create instance without valid parameters for testing instance methods


def test_filter_operation_methods():
    filter_operation_class = x.FilterOperation
    assert hasattr(filter_operation_class, 'PosLt')
    assert hasattr(filter_operation_class, 'PosGt')
    assert hasattr(filter_operation_class, 'CoverageGt')
    assert hasattr(filter_operation_class, 'Strand')
    assert hasattr(filter_operation_class, 'Context')


def test_region_reader_iterator_methods():
    # Note: This is typically obtained from RegionReader, not instantiated directly
    # Testing the class exists and has expected methods is sufficient
    region_reader_iterator_class = x.RegionReaderIterator
    # These would be tested on an actual instance from RegionReader


def test_region_reader_methods():
    region_reader_class = x.RegionReader
    assert hasattr(region_reader_class, '__init__')
    assert hasattr(region_reader_class, 'from_reader')
    # Note: Can't easily create instance without valid parameters for testing instance methods
