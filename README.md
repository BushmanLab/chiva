# cHIVa
Combined HIV integration site and genomic variant analysis pipeline.

# Install
`git clone https://github.com/BushmanLab/chiva.git`

`cd chiva`

`./install.sh`

# Testing
`chiva setup configs/HIV_test.config.hml`

`chiva run configs/HIV_test.config.hml`

# Processing Overview

This section provides an overview of each step in the processing pipeline.

## Demultiplexing

Reads are demultiplexed into R1 and R2 `.fastq.gz` files for each sample using a standard demultiplexing algorithm.

## Trimming

R1 reads for each sample are trimmed on the 5' end to remove the linker sequence used in that sample's library construction.  Reads not matching the linker sequence are filtered out.  The remaining reads are further trimmed on the 3' end to remove the reverse complement of the expected LTR sequence.

R2 reads for each sample are trimmed in three steps: First the 8nt "primer bit" is removed from the 5' end, and reads without a perfect match are filtered out.  Second, an approximate match of the expected LTR sequence is removed from the 5' end, and reads without an approximate match are filtered out.  Lastly, a CA sequence is removed from the 5' end, and reads without a CA at the 5' end are filtered out.  This last bit is done because the CA is known to be the end of the LTR and thus remaining sequence should be human (or internal viral genome).

## Filtering

Read pairs that did not pass BOTH sets of trimming filters (R1 and R2) are removed from downstream analysis.

## Consolidation

Unique R1 and R2 sequences are identified.  A key file is generated to map sequencer IDs onto each unique sequence identified.

## Mapping

R1 and R2 reads are mapped, independently, against the hg38 genome using BLAT.  Default parameters can be found in the `configs/HIV_test.config.yml` file.

## Int site identification

For each sample, BLAT output is analyzed using `tools/rscripts/couple.R` and a list of unique integration sites is generated.  We define a "unique site" as a unique pair of "anchor" (R2) and "adrift" (R1) sequences.

## Condensing int sites

The list of unique sites for each sample is condensed across replicates and regions (U3 and U5) for each sample.  Sites found in both U3 and U5 regions are marked as "Dual Detect".

## Filtering int sites

The final processing step is filtering of "crossover" sites which are sites suspected of spilling over from one sample to another.

