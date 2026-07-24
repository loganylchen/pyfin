#!/bin/bash
#SBATCH --job-name=stepA_gather
#SBATCH --output=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs/%x_%j.out
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
L=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/slurm_logs
echo "########## BYTE-IDENTITY (diff job 211279) ##########"
cat $L/m2a_diff_211279.out 2>/dev/null || echo "NO DIFF OUTPUT"
echo; echo "########## UNIT TESTS (job 211281) ##########"
grep -E "passed|failed|error|PYTEST_EXIT" $L/stepA_test_211281.out 2>/dev/null | tail -15 || echo "NO TEST OUTPUT"
echo "########## END ##########"
