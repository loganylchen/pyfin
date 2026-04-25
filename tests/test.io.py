import os

from fin.io import generate_isolated_intervals
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
    break
#   read_seqs = j[0]
#   ref_seqs = j[1]
#   for region in j[2]:
#     ids = region['ids']
#     three_prime_positions = region['3_positions']
#     break

# strand='-'
# read_seq_list = [(j, read_seqs[j],three_prime_positions[i]) for i,j in enumerate(ids) if not j.startswith('gtf_')]
# ref_seq_list = [(j.replace('gtf_',''),ref_seqs[j.replace('gtf_','')],three_prime_positions[i]) for i,j in enumerate(ids) if j.startswith('gtf_')]
# contained_read = set()

# for read_i in read_seq_list:
#     read_id, read_seq, read_3_end = read_i
#     for ref_i in ref_seq_list:
#         ref_id, ref_seq, ref_3_end = ref_i
#         if strand == '+':
#             end_dif = read_3_end- ref_3_end
#         else:
#             end_dif =  ref_3_end - read_3_end
#         if end_dif > 0:
#             if read_seq[:-end_dif] in ref_seq:
#                 print(f'----{read_id}:{end_dif}')
#                 print(read_seq[:-end_dif])
#                 print('===============')
#                 print(ref_seq)
#                 print(read_seq[:-end_dif] in ref_seq)
#                 contained_read.add(read_id)
#                 break
#         else:
#             if read_seq in ref_seq[:end_dif]:
#                 print(f'----{read_id}:{end_dif}')
#                 print(read_seq )
#                 print('===============')
#                 print(ref_seq[:end_dif])
#                 print(read_seq[:-end_dif] in ref_seq)
#                 contained_read.add(read_id)
#                 break
            


# potential_novels = sorted([i for i in read_seq_list if i[0] not in contained_read],key=lambda x: len(x[1]))
# print(len(contained_read))

# print(len(potential_novels))
# for i, read_i in enumerate(potential_novels):
#     read_id, read_seq, read_3_end = read_i
#     for other_i in potential_novels[i+1:]:
#         ref_id, ref_seq, ref_3_end = other_i
#         if ref_id == read_id:
#             continue
#         if strand == '+':
#             end_dif = read_3_end- ref_3_end
#         else:
#             end_dif =  ref_3_end - read_3_end
#         if end_dif > 0:
#             if read_seq[:-end_dif] in ref_seq:
#                 print(f'----{read_id}:{end_dif}')
#                 print(read_seq[:-end_dif])
#                 print('===============')
#                 print(ref_seq)
#                 print(read_seq[:-end_dif] in ref_seq)
#                 contained_read.add(read_id)
#                 break
#         else:
#             if read_seq in ref_seq[:end_dif]:
#                 print(f'----{read_id}:{end_dif}')
#                 print(read_seq )
#                 print('===============')
#                 print(ref_seq[:end_dif])
#                 print(read_seq[:-end_dif] in ref_seq)
#                 contained_read.add(read_id)
#                 break


# potential_novels = sorted([i for i in read_seq_list if i[0] not in contained_read],key=lambda x: len(x[1]))
# print(len(contained_read))

# print(len(potential_novels))