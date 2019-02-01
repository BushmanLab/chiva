#!/usr/bin/env bash
set -e

# Input arguments
__CHIVA_ENV=${1-chiva}
__CORES=${2-1}

# Clear test directory
rm -rf etc/tests/HIV_test

# Test script
source ${HOME}/miniconda3/etc/profile.d/conda.sh
conda activate ${__CHIVA_ENV}

# Create test analysis directory
chiva setup configs/HIV_test.config.yml -- -np
chiva setup configs/HIV_test.config.yml -- --nolock

# Generate test DAG graph
chiva run configs/HIV_test.config.yml -- -np

#chiva run configs/HIV_test.config.yml -- --dag --nolock | dot -Tsvg > \
#    etc/tests/HIV_test/output_data/HIV_test.dag.svg

chiva run configs/HIV_test.config.yml -- \
    -w 30 --nolock --keep-going --cores ${__CORES}

chiva report etc/tests/HIV_test/output_data/standardized_uniq_sites.rds \
    -c configs/HIV_test.config.yml \
    -o etc/tests/HIV_test/output_data/report.HIV_test.html \
    -t html
    
# Test for precise outputs
Rscript tools/rscripts/check_file_digests.R etc/tests/HIV_test.digests.yml -v
