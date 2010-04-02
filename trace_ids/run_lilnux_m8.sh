path="/net/direct/dnas01/arm2arm/DATA/LIA/SPHBAR/NFW/MODELS/"
snap=$1
model="MODEL8"
echo ./bin/trace_id.x --snapshotList=${path}/${model}/RUNG2/SNAPS/snap_m8_gal_sfr_bh_${snap} ${path}/${model}_BH/RUNG2/SNAPS/snap_m8_gal_sfr_bh_${snap} --out-file=${path}/${model}/ANALYSIS/OTHER/${model}_${snap}.idx 
