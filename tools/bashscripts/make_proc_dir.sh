# This is a short script to generate processing directories
# First (and only) argument should be the config file
#
set -e

RUN=$(grep 'Run_Name :' ${1} | sed 's/\s//g' | sed 's/"//g' | sed 's/Run_Name://')
DATADIR=$(grep 'Processing_Path' ${1} | sed 's/\s//g' | sed 's/"//g' | sed 's/Processing_Path://')
PROCDIR="${DATADIR}/${RUN}"

mkdir ${PROCDIR}
ln -s $(readlink -f ${1}) ${PROCDIR}/${RUN}.config.yml
mkdir ${PROCDIR}/input_data
mkdir ${PROCDIR}/logs
mkdir ${PROCDIR}/analysis_data
mkdir ${PROCDIR}/analysis_data/uniqSites
mkdir ${PROCDIR}/analysis_data/chimeras
mkdir ${PROCDIR}/analysis_data/multihits
mkdir ${PROCDIR}/output_data

echo "  ${RUN} run has been set up."
echo "  Please place input fastq.gz files if path is not available in config"
echo "  into the following directory:"
echo "    ${PROCDIR}/input_data"
