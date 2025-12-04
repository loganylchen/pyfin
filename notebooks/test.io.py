import os
import fin

from fin.io.io_read_manager import create_subset_manager,ReadSubsetManager
from fin.io import generate_isolated_intervals, extract_reads_for_interval
import logging
from fin.analysis import ThreePrimePositionClustering

# Set debug level (choose one of these options):

# Option 1: Set global debug level
# logging.basicConfig(level=logging.DEBUG)

# Option 2: Use the package logger with debug level
from fin.utils.log_config import setup_logger
logger = setup_logger(__name__, level='DEBUG', log_file='eventalign_debug.log')

TEST_DIR='./testdata'
bam_path=os.path.join(TEST_DIR,'test.bam')
gtf_path=os.path.join(TEST_DIR,'reference.gtf')
fasta_path= os.path.join(TEST_DIR,'reference.fa')
transcript_path=os.path.join(TEST_DIR,'transcript.fa')

result = generate_isolated_intervals(
      bam_path,
      gtf_path,
      tmp_dir='tmp_dir'
  )
intervals = result['intervals']
fusion_ids = result['fusion_read_ids']

# print(intervals)
# print(fusion_ids)

tppc = ThreePrimePositionClustering(bam_path,transcript_path)
# clustering_result =tppc.cluster_three_prime_positions(intervals[9].attrs,intervals[9].three_prime_pos)

for i in intervals:
  if i.read_count > 20 and i.read_count < 40 and i.strand=='-':
    break


for j in tppc.iter_interval([i]):
  read_seqs = j[0]
  ref_seqs = j[1]
  for region in j[2]:
    ids = region['ids']
    three_prime_positions = region['3_positions']
    break
