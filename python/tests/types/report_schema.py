from bsx2.types import ReportTypeSchema

def test_report_type_schema_bismark():
    schema_type = ReportTypeSchema.Bismark
    assert schema_type.col_names() == ["chr", "position", "strand", "count_m", "count_um", "context", "trinuc"]
    assert isinstance(schema_type.schema(), object)
    assert schema_type.chr_col() == "chr"
    assert schema_type.position_col() == "position"
    assert schema_type.context_col() == "context"
    assert schema_type.strand_col() == "strand"
    assert not schema_type.need_align()

def test_report_type_schema_cgmap():
    schema_type = ReportTypeSchema.CgMap
    assert schema_type.col_names() == ["chr", "nuc", "position", "context", "dinuc", "density", "count_m", "count_total"]
    assert isinstance(schema_type.schema(), object)
    assert schema_type.chr_col() == "chr"
    assert schema_type.position_col() == "position"
    assert schema_type.context_col() == "context"
    assert schema_type.strand_col() == "strand"
    assert not schema_type.need_align()

def test_report_type_schema_bedgraph():
    schema_type = ReportTypeSchema.BedGraph
    assert schema_type.col_names() == ["chr", "start", "end", "density"]
    assert isinstance(schema_type.schema(), object)
    assert schema_type.chr_col() == "chr"
    assert schema_type.position_col() == "start"
    assert schema_type.context_col() is None
    assert schema_type.strand_col() is None
    assert schema_type.need_align()

def test_report_type_schema_coverage():
    schema_type = ReportTypeSchema.Coverage
    assert schema_type.col_names() == ["chr", "start", "end", "density", "count_m", "count_um"]
    assert isinstance(schema_type.schema(), object)
    assert schema_type.chr_col() == "chr"
    assert schema_type.position_col() == "start"
    assert schema_type.context_col() is None
    assert schema_type.strand_col() is None
    assert schema_type.need_align()

def test_report_type_schema_equality():
    bismark1 = ReportTypeSchema.Bismark
    bismark2 = ReportTypeSchema.Bismark
    cgmap = ReportTypeSchema.CgMap

    assert bismark1 == bismark2
    assert bismark1 is bismark2 # Enums should be singletons

    assert bismark1 != cgmap
    assert bismark1 is not cgmap

    assert ReportTypeSchema.BedGraph != ReportTypeSchema.Coverage
    assert ReportTypeSchema.BedGraph is not ReportTypeSchema.Coverage

# The #[pyclass(eq_int)] makes the enum variants comparable to their
# underlying integer value (starting from 0).
def test_report_type_schema_equality_with_int():
    assert ReportTypeSchema.Bismark == 0
    assert ReportTypeSchema.CgMap == 1
    assert ReportTypeSchema.BedGraph == 2
    assert ReportTypeSchema.Coverage == 3

    assert 0 == ReportTypeSchema.Bismark
    assert 1 == ReportTypeSchema.CgMap
    assert 2 == ReportTypeSchema.BedGraph
    assert 3 == ReportTypeSchema.Coverage

    assert ReportTypeSchema.Bismark != 1
    assert ReportTypeSchema.CgMap != 0
    assert ReportTypeSchema.BedGraph != 3
    assert ReportTypeSchema.Coverage != 2
