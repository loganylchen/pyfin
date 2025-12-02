import os
import fin
from fin.io.io_read_manager import create_subset_manager,ReadSubsetManager
from fin.io import generate_isolated_intervals, extract_reads_for_interval

TEST_DIR='./testdata'
bam_path=os.path.join(TEST_DIR,'test.bam')
gtf_path=os.path.join(TEST_DIR,'reference.gtf')
fasta_path= os.path.join(TEST_DIR,'reference.fa')

result = generate_isolated_intervals(
      bam_path,
      gtf_path
  )
intervals = result['intervals']
fusion_ids = result['fusion_read_ids']

print(intervals)
print(fusion_ids)