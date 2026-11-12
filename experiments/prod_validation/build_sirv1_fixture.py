#!/usr/bin/env python3
"""Freeze one SIRV locus, WITH its raw signal, as the end-to-end test fixture.

Takes every alignment record on SIRV1 of SG-NEx H9 direct RNA replicate4 run2
(aligned to the SIRV Set 4 reference with minimap2 `-ax splice -uf -k14`), the
raw signal of exactly those reads, the SIRV1 contig and its eight annotated
transcripts. Nothing is edited: coordinates, CIGARs, tags, read sequences and
signal samples are the originals. The locus is small enough to run in seconds
and still has exact alignment-score ties that the signal has to break, which is
the path a sequence-only fixture cannot reach.

Run inside the PyFIN image (needs pysam and pyslow5):

    python3 experiments/prod_validation/build_sirv1_fixture.py
"""
from pathlib import Path

import pysam
import pyslow5

GRID = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex")
SAMPLE = "SGNex_H9_directRNA_replicate4_run2"
SRC = GRID / "results" / "sirv4_full_sweep_ratio_grid_v2" / SAMPLE
REF = GRID / "refs" / "sirv"
DST = Path(__file__).resolve().parents[2] / "tests" / "data" / "sirv1_mini"
CONTIG = "SIRV1"


def main():
    DST.mkdir(parents=True, exist_ok=True)

    with pysam.FastaFile(str(REF / "SIRV_Set4.fa")) as fasta, (DST / "reference.fa").open("w") as out:
        seq = fasta.fetch(CONTIG)
        out.write(f">{CONTIG}\n")
        for i in range(0, len(seq), 80):
            out.write(seq[i:i + 80] + "\n")
    pysam.faidx(str(DST / "reference.fa"))

    n_tx = 0
    with (REF / "SIRV_Set4.gtf").open() as src, (DST / "annotation.gtf").open("w") as out:
        for line in src:
            if line.split("\t", 1)[0] == CONTIG:
                out.write(line)
                n_tx += line.split("\t")[2] == "transcript"

    read_ids = []
    with pysam.AlignmentFile(str(SRC / "align" / f"{SAMPLE}.sorted.bam"), "rb") as reader:
        header = reader.header.to_dict()
        header["SQ"] = [sq for sq in header["SQ"] if sq["SN"] == CONTIG]
        # the original @PG line carries the aligner command with site-local paths
        header.pop("PG", None)
        with pysam.AlignmentFile(str(DST / "reads.bam"), "wb", header=header) as writer:
            for record in reader.fetch(CONTIG):
                clone = pysam.AlignedSegment.from_dict(record.to_dict(), writer.header)
                writer.write(clone)
                if record.query_name not in read_ids:
                    read_ids.append(record.query_name)
    pysam.index(str(DST / "reads.bam"))

    source = pyslow5.Open(str(SRC / "subset" / "mapped.blow5"), "r")
    target_path = DST / "reads.blow5"
    for stale in (target_path, Path(str(target_path) + ".idx")):
        if stale.exists():
            stale.unlink()
    target = pyslow5.Open(str(target_path), "w", rec_press="zlib", sig_press="svb-zd")
    # only the attributes the source actually carries: an attribute left at
    # None cannot be written
    target.write_header(
        {key: str(value) for key, value in source.get_all_headers().items() if value is not None},
        end_reason_labels=source.get_aux_enum_labels("end_reason"),
    )
    aux_names = source.get_aux_names()
    written = 0
    for read_id in sorted(read_ids):
        read = source.get_read(read_id, aux="all")
        if read is None:
            raise SystemExit(f"{read_id} is aligned on {CONTIG} but absent from the signal file")
        record, aux = target.get_empty_record(aux=True)
        for key in record:
            record[key] = read[key]
        for key in aux_names:
            aux[key] = read[key]
        target.write_record(record, aux)
        written += 1
    target.close()
    source.close()
    # opening for read writes the index next to the file; ship it, so a
    # read-only checkout can run the fixture
    check = pyslow5.Open(str(target_path), "r")
    ids, count = check.get_read_ids()
    check.close()
    if count != written or sorted(ids) != sorted(read_ids):
        raise SystemExit("the written signal file does not hold exactly the aligned reads")

    print(f"{DST}: {CONTIG}, {n_tx} annotated transcripts, {len(read_ids)} reads, {written} signals")


if __name__ == "__main__":
    main()
