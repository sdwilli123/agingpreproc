%%
runs = [14 20 23 26];
stems= {'Match', 'Breathhold', 'Match', 'Match'};
mkdir stcfsl_mc2_hran_edithd
for jj=1:length(runs)
runnum=runs(jj); 
nme = cell2mat(stems(jj)); 
base = '/projectnb2/fastfmri/sdwilli/aging/221104_ag115a/';
fn = ['run_' num2str(runnum) '_' nme '_SMS_CMRR_2.5mm_S8pe4_TR378_stc_mc2.nii'];
niftiFile=[ base 'stcfsl_mc2/' fn];
V_info = niftiinfo(niftiFile);
r =num2str(jj+1); 
%
newInfo = niftiinfo([base 'stcfsl_mc2_hran/hran_run0' r '.nii']);
V_array = double(niftiread([base 'stcfsl_mc2_hran/hran_run0' r '.nii']));

newInfo.PixelDimensions=V_info.PixelDimensions(1:4);
newInfo.BitsPerPixel=V_info.BitsPerPixel;
newInfo.TimeUnits=V_info.TimeUnits;
newInfo.SpaceUnits=V_info.SpaceUnits;
newInfo.SpatialDimension=V_info.SpatialDimension;
newInfo.TransformName=V_info.TransformName;
newInfo.Transform=V_info.Transform;

%
new_volume_name = [ 'run' r '_stc_mc2_hran.nii'];
niftiwrite(V_array, [ base 'stcfsl_mc2_hran_edithd/' new_volume_name], newInfo);
end 
%end
