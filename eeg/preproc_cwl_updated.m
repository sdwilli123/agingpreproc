%version sept 2022
%sw


%   preprocesses CWL EEG data collected inside the scanner
%       GAC: AAS
%       BCG: CWL regression
%       REREF
%% Load dependencies ------------
addpath(genpath('/ad/eng/research/eng_research_lewislab/users/sdwilli/scripts/'))
addpath(genpath('/usr2/postdoc/sdwilli/Downloads/chronux_2_12'))
addpath('/usr2/postdoc/sdwilli/Desktop/eeglab14_1_2b')
eeglab
%% make output folders
mkdir QC_figures
mkdir gac
mkdir gac_cwlreg
mkdir gac_cwlreg_aligned

%% call preproc function on brainvision files 
runnumber= 2; 
high =-1; % 1= high, 0 = low -1 = if not HIGH VS.LOW
d ='/ad/eng/research/eng_research_lewislab/data/emot/mood23_TSD_110522/eeg/scan/';
%'/ad/eng/research/eng_research_lewislab/data/CSF_aging/ag108a_daytime_221012/eeg/';
if high==1
    filename = sprintf('run%i_high.vhdr', runnumber); 
elseif high==0
    filename = sprintf('run%i_low.vhdr', runnumber); 
else 
    %manually enter filename:
    runnumber= 6; 
    filename = 'run6-rest.vhdr'; %sprintf('run%i.vhdr', runnumber);
    %for MRS we are just using CWL 
end
preproc_cwl(filename, runnumber, high,d) 


function preproc_cwl(filename, runnumber, high, d)
EEG =pop_loadbv(sprintf(d), filename);
r = runnumber; 

if high <2
trigs = []; 
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type, 'B  1')
  % if strcmp(EEG.event(i).type, 'S  1')
        trigs(end+1) = EEG.event(i).latency; 
    end
end

disp(['NUM TRIGS: ' num2str(length(trigs))])
save(['gac/run' num2str(r) '_trig_fs5000'], 'trigs');

trigs = trigs(1:end);
startime = trigs(1)/5000; %units = sec
stoptime = trigs(end)/5000;
tottime = (stoptime-startime)/60;


% store basic information of raw data
eeg = EEG.data; 
trig =trigs;
numtrig = length(trig);
t0=1;
tf = length(eeg);


% load scan triggers and prepare vars ------------------------
vars.UseGAC = 1;
vars.ntr = numtrig;
vars.currentPosition = 1;
vars.OrigChunk = eeg(:, t0:tf); 
[vars.SamplesInChunk, ChansInChunk] = size(vars.OrigChunk');
vars.currentPosition_orig = t0; 


% GAC using custom func
% how many TRs to use as template
trss = 20; %20 %def: number of TRs to average across 
vars.nTRtemp =trss;
vars.trGap= 1890; %units = samplesn(fs*TR = 5000*.378)
vars.OrigChunk = eeg;
vars.Clean = vars.OrigChunk';
%
%trig_uniform = trigs(1): vars.trGap: (trigs(1)+(1890*numtrig)-1);
%
if vars.UseGAC % Perform GAC
    tic
    for startTR = 1:((vars.ntr)-vars.nTRtemp)
        % get the data from the last n TRs, and reshape so its vars.trGap
        % samples X M channels X vars.ntr epochs
        vars.Positionstrat = trigs(startTR); %units = samples 
        vars.Positionend = trigs(startTR) +(vars.trGap*vars.nTRtemp); %changed this line 
        vars.chan = size(eeg,1);
        EEG.Recording_orig = eeg(:,vars.Positionstrat:vars.Positionend-1)';
        template = reshape(EEG.Recording_orig, vars.trGap, vars.nTRtemp, vars.chan); %
        template = permute(template, [1, 3, 2]);

        % get template average across n TRs
        MeanTemplate = mean(template, 3); % across TRS 

        % basically just subtract the template from our data - this
        % implemnation works in real time, your implementation will differ
        % offline
        if vars.SamplesInChunk <= vars.trGap %if for some reason smaller than a TR, then shrinking template to match the chunk
            vars.OrigChunk(1:end - 1, :) = vars.OrigChunk(1:end - 1, :) - ...
                MeanTemplate(end - vars.SamplesInChunk + 1:end, :)';
        else  %normal cases 
%             numReps = floor(vars.SamplesInChunk / vars.trGap) + 1;
%             numReps = vars.n;
            MeanTemplate = repmat(MeanTemplate, vars.nTRtemp, 1);
%             vars.OrigChunk = EEG.Recording_orig;
            for i = 1:vars.chan
                vars.Clean(vars.Positionstrat:vars.Positionend-1,i) = vars.OrigChunk(i,vars.Positionstrat:vars.Positionend-1)' - MeanTemplate(:,i);end 
        end
    end
end
toc


%
params = struct;
params.Fs = 5000;
params.tapers=[3 5]; 
params.fpass = [ 0 20];
movingwin = [10 1];
for ch= [33]

clear S t f
[S, t, f] = mtspecgramc_SW(vars.Clean(t0:tf,ch)  , movingwin, params);
figure; subplot(122); imagesc(t/60,f,10*log10(S)');axis xy;colormap jet;colorbar;set(colorbar,'visible','on'); caxis([-10 30])
title('GAC')
clear S t f
[S, t, f] = mtspecgramc_SW(eeg(ch,t0:tf), movingwin, params);
subplot(121); imagesc(t/60,f,10*log10(S)');axis xy;colormap jet;colorbar;set(colorbar,'visible','on'); caxis([-10 30])
title('raw')

sgtitle([ EEG.chanlocs(ch).labels])

end
%
%
saveas(gcf, ['QC_figures/run' num2str(r)  EEG.chanlocs(ch).labels 'gac3TRS.jpg']);
% save out GAC .mat file as intermediate step

%

EEG1 = EEG; 
EEG1.data = vars.Clean'; 
EEG=EEG1;

save(['gac/run' num2str(r) '_gac_fs5000'], 'EEG');

%
EEG = pop_loadset(['gac/run' num2str(r) '_gac_fs5000']);
end

EEG_res = pop_resample(EEG, 500);

%

% 
% cwl regression
win =25 ; %14 % 40 %14 % 25 %14.5 ; 
lag = 0.09 ;%0.126 %0.09 %0.126 %021; %0.09% 0.126;
fs = 500; 
eeg_cwl = pop_cwregression(EEG_res,fs,win,lag,1,'hann', [32:36]  , 1:32,'taperedhann',0);


save([ 'gac_cwlreg/run' num2str(r) '_gac_fs500_cwlreg'], 'eeg_cwl');


%

wi = [ 25 2]; 
for cnum =  [9 10 20 22 ] 
params=struct;
params.Fs =500; 
params.tapers=[5 9];
params.fpass=[0 25];
cleandat = eeg_cwl.data; 
ch = cleandat(cnum,:); t0=1;
[spec t f]=mtspecgramc_SW(ch(t0:end)  , wi,params);
figure
subplot(121)
imagesc(t/60,f,log10(abs(spec'))); 
axis xy; title([ 'gradient + cwl corrected win ' num2str(win) 's lag ' num2str(lag*1000) 'ms']); 
colormap('jet')
caxis([-2 2])

%
params.Fs =500; 
ch =EEG_res.data(cnum,:); 
[spec t f]=mtspecgramc_SW(ch(t0:end) , wi,params);
subplot(122)
imagesc(t/60,f,log10(abs(spec'))); 
colormap('jet')
xlabel('Time(min)')
sgtitle(num2str(cnum))
axis xy; title( 'gradient corr');
caxis([-2 2])

%sgtitle(['ch' num2str(cnum) ' mood11  run' num2str(r)])
%colorbar()

saveas(gcf, [ 'QC_figures/run' num2str(r) 'gac20TRS__cwlwecg'  num2str(cnum) '.jpg']);
end
%end
%  Aligned 
trig = trigs;
dat =eeg_cwl.data;
fs = 500; 
samplestart = trig(1)/5000*500;
%samplestart = startime*fs; 
samplestop = trig(end)/5000*500;
%stoptime*fs; 
%

aligned_cwl = dat(:, samplestart:samplestop);
%mkdir gac_cwlreg_aligned
save(['gac_cwlreg_aligned/r' num2str(r) '_aligned_gac_cwl_fs500'], 'aligned_cwl')



%
%aligned_cwl = dat(:,130420:end);

%
averef = mean(aligned_cwl(1:22,:));
m1 = aligned_cwl(25,:); 
m2 = aligned_cwl(26,:);
m = mean([m1;m2]);
%
t = (1:length(averef))/500; 
eyes = []; cnt=0;
figure; hold on; 

    eyes(1,:) = smoothdata(aligned_cwl(29,:) - averef, 'movmean', 500); 
    plot(t/60, eyes(1,:)' ); 
     eyes(2,:) = smoothdata(aligned_cwl(28,:) - averef, 'movmean', 500); 
    plot(t/60, eyes(2,:)' + cnt*200); 


%
for c = [ 9 20 22 ] %[3 5 9 10 20 22 ]
params = struct;
params.Fs = 500;
params.fpass = [ 0 25];
params.tapers=[2 3];  %aligned_cwl(6,:)
movingwin = [5 2]; 
ch = aligned_cwl(c,:) -averef;

[S, t, f] = mtspecgramc_SW(ch , movingwin, params);
figure;
imagesc(t/60,f,log10(S)');axis xy;colormap jet;colorbar;set(colorbar,'visible','on'); 
%caxis([-2 2])
title(num2str(c))
end
%
mkdir(sprintf('gac_cwlreg_aligned_averef/run0%i',r))


%c=20;
for c=1:22
ch = aligned_cwl(c,:) -averef;
save(['gac_cwlreg_aligned_averef/run0' num2str(r) '/ch' num2str(c) ] ,'ch' )
end



end




%%%%%%%%%%%%%%%%%%%%%% end %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Eye channels in EEG
% figure; hold on
% for ch  = [29 30]
% tt = (1:length(aligned_cwl))/500; 
%  plot(tt(100:end)/60, smoothdata(aligned_cwl(ch,100:end), 'movmean', 9000))
% end
%%
%r = 7;
%ave = mean(aligned_cwl([6:8 11:19 21:25],:)); %
%ave = mean(aligned_cwl([1:5],:)); 
% % params=struct;
% % params.Fs =500; 
% % params.tapers=[2 3];
% % params.fpass=[0 16];
% % movingwin = [25 2];
% % % Fz Fc1 Cz C3 p3 O1 o2  pz poz oz
% % for ch =  [3 4 5 9 10 20 22 ]
% %     %[ 1 3 5 4 6 8 17 18 19 9 10 20 ] %17 21 18 5 7 9 10 19 20]
% %         params.Fs = 500; 
% %         [S, t, f] = mtspecgramc_SW(aligned_cwl(ch,:) -m   , movingwin, params);
% %         figure;
% %         imagesc(t/60,f,10*log10(S)');axis xy;
% %         colormap jet;colorbar;set(colorbar,'visible','on');
% %         caxis([-10 20])
% %         title(num2str(ch))
% %         %saveas(gcf, [ 'QC_figures/run' num2str(r) 'gac20TRS__cwl_aligned'  num2str(ch) '.jpg']);
% % 
% % end
% % %% createa matrix to score
% % %  fz cz pz o2 e1 e2 
% % dat = [];
% % %chans =[17 18 19 10 29 30]; mood14
% % chans = [3 5 20]; 
% % for ch = chans
% % dat(end+1,:) = aligned_cwl(ch,:) - averef;
% % end
% % %
% % dat = vertcat(dat, eyes);
% % %
% % mkdir toscore
% % save([ 'toscore/r' num2str(r) '_eeg' ], 'dat')
% % %%
% % t = (100:length(aligned_cwl))/500; 
% % for ch = 27:31
% % figure; plot(t, smoothdata(aligned_cwl(ch,100:end), 'movmean', 1500))
% % title(num2str(ch))
% % end
%%

         