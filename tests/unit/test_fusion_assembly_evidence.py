"""Final assembly membership is auxiliary structure, not signal confirmation."""

import pytest

from fin.fusion.assembly_evidence import load_assembly_junctions
from fin.fusion.annotation import AnnotatedReadPath, FusionArm
from fin.fusion.records import AlignmentBlock
from fin.fusion.correction import junction_locus_key


def _files(tmp_path):
    gtf = tmp_path / "assembly.gtf"
    gtf.write_text(''.join(f'{chrom}\tpyfin\texon\t{start}\t{end}\t.\t{strand}\t.\tgene_id "{chrom}"; transcript_id "{chrom}.tx";\n'
                           for chrom, strand in (("A", "+"), ("B", "-")) for start, end in ((101, 200), (301, 400))))
    support = tmp_path / "junction_support.tsv"
    support.write_text("chrom\tstrand\tinterval\tdonor\tacceptor\traw_reads\tindel_clean_reads\n"
                       "A\t+\tA:0-500\t200\t300\t8\t4\n"
                       "A\t+\tA:1-500\t200\t300\t8\t3\n"
                       "A\t+\tA:0-500\t201\t300\t100\t100\n"
                       "B\t-\tB:0-500\t200\t300\t8\t5\n")
    return gtf, support


def test_requires_final_model_membership_and_does_not_sum_duplicate_rows(tmp_path):
    evidence = load_assembly_junctions(*_files(tmp_path))
    assert evidence.junctions == (("A", "+", 200, 300, 4), ("B", "-", 300, 200, 5))


def test_fusion_endpoint_support_is_not_exact_fusion_support(tmp_path):
    evidence = load_assembly_junctions(*_files(tmp_path))
    left = FusionArm("A", "+", (AlignmentBlock(100, 200, 0, 100, 100),), ("A",))
    right = FusionArm("B", "-", (AlignmentBlock(100, 200, 100, 200, 100),), ("B",))
    path = AnnotatedReadPath("read", (left, right), "interchromosomal_candidate", "resolved")
    audit, = evidence.annotate(path)
    assert audit["donor_clean_reads"] == 4
    assert audit["acceptor_clean_reads"] == 5
    assert audit["exact_intron_clean_reads"] == 0
    assert evidence.registry([path])[junction_locus_key("fusion", left, right)] == ((200, 200),)


def test_clean_threshold_drops_weak_final_junctions(tmp_path):
    evidence = load_assembly_junctions(*_files(tmp_path), min_clean_reads=5)
    assert evidence.junctions == (("B", "-", 300, 200, 5),)


def test_malformed_support_fails_instead_of_promoting(tmp_path):
    gtf, support = _files(tmp_path)
    support.write_text(support.read_text().replace("8\t4", "1\t4"))
    with pytest.raises(ValueError, match="counts"):
        load_assembly_junctions(gtf, support)
