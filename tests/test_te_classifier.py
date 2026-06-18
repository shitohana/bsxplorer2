import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parents[1] / "scripts" / "oryza_dmr_functional_downstream"
sys.path.insert(0, str(SCRIPT_DIR))

from pipeline_core import classify_te_like


def status(text):
    return classify_te_like(text)[0]


def test_retrotransposon_ty3_gypsy_is_te_like():
    assert status("retrotransposon protein, Ty3-gypsy subclass") == "TE-like"


def test_transposon_is_te_like():
    assert status("transposon protein") == "TE-like"


def test_centromere_specific_retrotransposon_is_te_like():
    assert status("centromere-specific retrotransposon protein") == "TE-like"


def test_adenylate_kinase_is_non_te():
    assert status("adenylate kinase") == "non-TE protein-coding"


def test_glutathione_s_transferase_is_non_te():
    assert status("glutathione S-transferase") == "non-TE protein-coding"


def test_expressed_protein_is_unknown():
    assert status("expressed protein") == "unknown/hypothetical"


def test_hypothetical_protein_is_unknown():
    assert status("hypothetical protein") == "unknown/hypothetical"


def test_protein_kinase_not_te_like():
    assert status("protein kinase") == "non-TE protein-coding"
