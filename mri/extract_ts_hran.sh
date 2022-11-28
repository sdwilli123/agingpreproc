# run from inside the stcfsl_mc2_hran_edithd directory 

export masknumbers=("11" "14" "20" "23" "26")
export bp=export bp=${PWD%/*}
mkdir ${bp}/ts_hran
export fsn=fs_recon_biascorr
export fp=${bp}/$fsn
for fn in *.nii; do
fname=`$FSLDIR/bin/remove_ext ${fn}`; # store filename
r="$(echo ${fname} | cut -d'_' -f1)"
r="${r:3:1}"
echo ${r}

export masknum=${masknumbers[$r-1]}
echo ${masknum}

fslmeants -i ${bp}/stcfsl_mc2_hran_edithd/${fn} -m ${bp}/masks/run${masknum}_cortex.nii -o ${bp}/ts_hran/run${r}_cortex.txt

fslmeants -i ${bp}/stcfsl_mc2_hran_edithd/${fn} -m ${bp}/masks/r${masknum}_v.nii -o ${bp}/ts_hran/run${r}_v.txt

fslmeants -i ${bp}/stcfsl_mc2_hran_edithd/${fn} -m ${bp}/masks/run${masknum}_v1.nii -o ${bp}/ts_hran/run${r}_v1.txt


done
