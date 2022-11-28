#!/bin/bash
# Usage: bias_field_start.sh input_file output_directory

# Safety precautions
# Using a variable that is not set will raise an error
#make sure file ends in T1
module load spm/12
module load afni
# load module for 3dcalc command! Forgot what it is -

set -o nounset
# Exit if a command fails.
set -o errexit

# Assumes that the FreeSurfer recon follows the same name as the initials of the input file.
# And also assumes that SUBJECTS_DIR is set correctly.

usage() {
    echo "Bias field correct an image using SPM."
    echo ""
    echo "Usage: bias_field_start.sh biasedimg outputdir"
    echo ""
    echo "Assumes:"
    echo "SPM loaded on your matlab startup.m file."
}

if [ $# -lt 2 ] ; then
    usage
    exit 0
fi

if [ "${1: -3}" == ".gz" ]
then
    cp "$1" uncorr.nii.gz
    # The SPM batch script expects uncorr.nii
    gzip -d uncorr.nii.gz
else
    cp "$1" uncorr.nii
fi

corrname=$(basename "$1" | sed "s/T1w.nii.gz/preproc-bico_T1w.nii/")
echo "$1" "$2/$corrname"

echo "SPM motion batch starting."

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

if ! test -f ./Bias_field_script_job.m; then
    cp "${DIR}"/Bias_field_script_job.m ./Bias_field_script_job.m
fi

# Run the batch but afterwards change the output from double to short int.
matlab -nodesktop -nosplash -r "Bias_field_script_job"
3dcalc -a muncorr.nii -prefix muncorr.nii -overwrite -expr 'a' -datum short

mv muncorr.nii "$2"/"$corrname"
gzip "$2"/"$corrname"

rm BiasField_uncorr.nii uncorr_seg8.mat uncorr.nii
rm c*uncorr.nii

echo "und tschuess"
