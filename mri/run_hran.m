%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%   new function version 
%   version sep 2022 SW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%inputs 
% r - specifies output run number (1 ... N) 
% fn =  %input functional sequence name as a string
% wm = 1; %boolean that specifies whether to use WM  (1) or brain mask (0)

%output
%saves clean nifti file as 'hran_run#r.nii' + 2 figs showing res
%mood17 run 11 only has 200 timepoints because subject was confused 
runs = [ 11 13 17 19]; 
typeR = { 'checker', 'rest', 'rest', 'rest'}; 
newruns =1:5;
base = '/projectnb2/fastfmri/sdwilli/aging/ag106b/';

for i =1:4
close all;clc
runnum =runs(i); 
r = num2str(newruns(i)); 
fn = ['run_' num2str(runnum) '_' cell2mat(typeR(i)) '_SMS_CMRR_2.5mm_S8pe4_TR378_stc_mc2.nii'];
WM = 0; 
hran_(base, fn, r, runnum, WM)
end 









function hran_(base, fn, r, runnum, WM)
%%
addpath(genpath('/usr2/postdoc/sdwilli/Downloads/HRAN-master'))
addpath(genpath('/usr2/postdoc/sdwilli/Downloads/chronux_2_12'))
%
delete(gcp('nocreate'))
NSLOTS = str2num(getenv('NSLOTS'))
maxNumCompThreads(NSLOTS); 
parpool('local',20)
%%

niftiFile=[ base 'stcfsl_mc2/' fn];
V_info = niftiinfo(niftiFile);

V = double(niftiread(V_info));
%
addpath(genpath('/ad/eng/research/eng_research_lewislab/users/bsetzer/scripts/Preprocessing/EEG/fieldtrip-20191025'));
vox =load_nifti(niftiFile);
v1=reshape(vox.vol,[prod(vox.dim(2:4)),vox.dim(5)]);
meanv=squeeze(mean(v1));
%
if WM
   WM_mask = boolean(niftiread([base 'masks/run' num2str(runnum) '_WM_mask.nii']));
    physiologicalData = squeeze(sum(WM_mask.*V,[1 2 3]))./sum(WM_mask,[1 2 3]);
else 
    physiologicalData =double(meanv)'; 
end 

rmpath(genpath('/ad/eng/research/eng_research_lewislab/users/bsetzer/scripts/Preprocessing/EEG/fieldtrip-20191025'));
addpath(genpath(fullfile('/ad/eng/research/eng_research_lewislab/users/ziy027/scripts/chronux_2_12')));
inputParams = struct;

%  Set TR, moving window length, percent overlap
TR = .378; %****have to set it manually bc header is WRONG**** V_info.PixelDimensions(4); %0.347s in the demo data
inputParams.TR = TR; %TR of fMRI data (s)
inputParams.windowLength = 40; % length of moving window (s) (e.g. 24-30s)
inputParams.percentOverlap = .75; %.75 % percent overlap of windows (e.g. .5 or .75)

% Set neural regressors 

time = [0:inputParams.TR:size(V,4)*inputParams.TR-inputParams.TR];
neuralFrequency = 1.2; % Hz
neuralSignal = 60*cos(2*pi*neuralFrequency.*time);
neuralZ = [cos(2*pi*neuralFrequency.*time)';sin(2*pi*neuralFrequency.*time)'];
inputParams.neuralZ = neuralZ;

% Set physiological frequency range
inputParams.cardiacFreqRange = [40:85]; % cardiac fundamental freq range (bpm)
inputParams.respFreqRange = [9:21];  % bpm, up to .4 is generous 

%Set model orders
inputParams.P_freq = 1; % AR order used to estimate physio frequencies (e.g. 1-4)
inputParams.R_freq = 1; % Resp order used to estimate physio frequencies (e.g. 1-4)
inputParams.C_freq = 1; % Cardiac order used to estimate physio frequencies (e.g. 1-3)
inputParams.N_freq = 0; % Number of neural regressors used to estimate physio freqs (e.g. 0 bc assuming no neural signal in this region)
inputParams.X_freq = 0; % Interaction order used to estimate physio frequencies (e.g. 0-1)

% 5) Set model orders used to create design matrix to regress physio noise
% from each voxel in the brain
inputParams.P_Z = 4; % AR order used to regress physio noise (e.g. 1-4)
inputParams.R_Z = 4; % Resp order used to regress physio noise (e.g. 1-4)
inputParams.C_Z = 4; % Cardiac order used to regress physio noise (e.g. 1-3)
inputParams.N_Z =0; % size(neuralZ,2); % Number of neural regressors used to regress physio noise (e.g. number of columns in neuralZ)
inputParams.X_Z = 0; % Interaction order used to regress physio noise (e.g. 0-1)

% Toggle for waitbar on or off
waitbarBoolean = 1; % 1 - on, 0 - off

% Run function to estimate physio!!
HRAN_estFreqs_outputParams = HRAN_estimatePhysFreqsSW(physiologicalData,inputParams,waitbarBoolean);
%this only took ~ 2 min to run 
% Plot the estimates
fig = figure('Position',[1 1 800 600]);
% Spectrogram of physiological data
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(physiologicalData), movingwin, params);
hold on
imagesc(stime, sfreq, 10*log10(spec_original)')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_hr_hat,'w','filled')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_rr_hat,'w','filled')
hold off
set(gca,'YDir','normal')
colormap('jet')
colorbar
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Physiological Noise ROI')
%
if WM 
    saveas(gcf, [ base 'physio/run' r 'estimatesWM.jpg'])
else 
   saveas(gcf, [ base 'physio/run' r 'estimatesGM.jpg'])
end
%% GET READY TO ITERATE THROUGH DATA 
% Store dimensions
[xDim,yDim,zDim,tDim] = size(V);

% Initialize arrays to store de noised data
deNoisedV = zeros(size(V));

disp(runnum)
% also need to register the brainmask to this run 
brainMask = boolean(niftiread([base 'masks/run' num2str(runnum) '_brain_mask.nii']));


% Iterate through data
% I should have edited this to display which z dimension we are on to have an idea
% of how long things will take - has been running for 30 min 
tic
parfor z = 1:zDim % NOTE: this is a parfor loop to take advantage of parallelization, can replace with regular loop
    disp(z)
    for y = 1:yDim
        for x = 1:xDim
            if brainMask(x,y,z) % NOTE: this line not required, but can speed up de-noising (by only de-noising voxels in the brain)
                %['Voxel: ' num2str(x) ',' num2str(y) ',' num2str(z)] % should comment out this line when running!
                deNoisedVoxel = HRAN_regressPhysNoise(squeeze(V(x,y,z,:)),HRAN_estFreqs_outputParams);
                deNoisedV(x,y,z,:) = deNoisedVoxel;
            end
        end
    end
end
toc

% Display amygdala with and without physio noise
% runnum=19;
% amgMask = boolean(niftiread('/projectnb/fastfmri/sdwilli/mood/mood5sd/masks/r19_lamyg.nii.gz'));
% %amgMask = boolean(niftiread( [baseM 'masks/R_amygdala_run' num2str(runnum) '.nii.gz']));
 originalamg = squeeze(sum(brainMask.*V,[1 2 3]))./sum(brainMask,[1 2 3]);
 deNoisedamg = squeeze(sum(brainMask.*deNoisedV,[1 2 3]))./sum(brainMask,[1 2 3]);

% Plot time series and spectra 
fig = figure('Position',[1 1 800 600]);
%sgtitle(['Run ' num2str(runnum-7) ' physio removal amygdala'])
% Time Series
subplot(2,3,[1 2 3])
hold on
plot(time(15:end), detrend(originalamg(15:end)/mean(originalamg(15:end))), 'k','LineWidth',2)
plt = plot(time(15:end), detrend(deNoisedamg(15:end)/mean(deNoisedamg(15:end))),'Color',[152,78,163]./256,'LineWidth',2);
plt.Color(4) = .6;
hold off
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
legend('Original brain','De-noised brain')

% Power Spectra
subplot(2,3,4)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[P_orig,f] = mtspectrumc(detrend(originalamg),params);
[P_deNoised,f] = mtspectrumc(detrend(deNoisedamg),params);
hold on
plot(f,10*log10(P_orig),'k','LineWidth',2)
plt = plot(f,10*log10(P_deNoised),'Color',[152,78,163]./256,'LineWidth',2);
plt.Color(4) = .6;
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
hold off

% Spectrograms
% Original
subplot(2,3,5)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(originalamg), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
%colorbar
%caxis([-30 15])
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Original')

% De-noised
subplot(2,3,6)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(deNoisedamg), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
%colorbar
%caxis([-30 15])
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('De-noised')
%
if WM
    saveas(gcf, [ base 'physio/run' r 'final_resultsWM.pdf'])
else
    saveas(gcf, [ base 'physio/run' r 'final_resultsGM.pdf'])
end 

% Save de-noised file 
if WM
    fileName_denoised = [base 'stcfsl_mc2_hran/hranWM_run0' r]; 
    niftiwrite(deNoisedV,fileName_denoised );
else 
    fileName_denoised = [base 'stcfsl_mc2_hran/hranGM_run0' r]; 
    niftiwrite(deNoisedV,fileName_denoised );

end 









end

