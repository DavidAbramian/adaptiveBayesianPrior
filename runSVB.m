%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Run SVB estimation for BIDS format fMRI data
%
% INPUT:        dS - data settings struct
%               VBMethod - string containing estimation method,
%                                  possible options: SPMsVB2d, SPMsVB3d,
%                                  ImprovedVB2d, ImprovedVB3d.
%               
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      2017-10-26
%
%
% Modified:     2019-10 David Abramian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runSVB(dS,VBMethod)
%% Run

% Read settings
dataPath = dS.dataPath;
outputPath = dS.outputPath;
subjStr = dS.subjStr;
fileStr = dS.fileStr;
TR = dS.TR;
condStr = dS.condStr;
contrastStr = dS.contrastStr;

actualSubjStr = ['sub-',subjStr];

% % SPM Paths
% SPMPath2 = strcat(SPMPath,'/svb');
% SPMPath3 = strcat(SPMPath2,'/',VBMethod);
% addpath(SPMPath);
% addpath(SPMPath2);
% addpath(SPMPath3);

% BOLD path
dirBOLD = fullfile(dataPath, actualSubjStr, 'func');

% Results path
resultsPath = fullfile(outputPath, actualSubjStr, VBMethod);
if ~exist(resultsPath,'dir')
    mkdir(resultsPath);
end

% Delete old SPM file
SPMFile = fullfile(resultsPath,'SPM.mat'); 
if exist(SPMFile,'file')
    delete(SPMFile) 
end

% Initialise SPM defaults
spm('Defaults','fMRI');
spm_jobman('initcfg');

% Number of fMRI volumes
load(fullfile(dirBOLD, [fileStr,'bold.mat']),'mat');
nScans = size(mat,3); 
clear mat;

% BOLD data
clear matlabbatch
matlabbatch{1}.spm.stats.fmri_spec.dir = cellstr(resultsPath);
matlabbatch{1}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
scans = cell(nScans,1);
for t = 1:nScans
  scans{t} = fullfile(dirBOLD, ['r',fileStr,'bold.nii,',num2str(t)]);
end
matlabbatch{1}.spm.stats.fmri_spec.sess.scans = scans;

% Conditions
condCell = readTaskFile(fullfile(dataPath, condStr)); 
eventCell = readTaskFile(fullfile(dirBOLD, [fileStr,'events.tsv'])); 

conditions = string(deblank(condCell(:,3)));
events = string(deblank(eventCell(:,3)));

nCond = size(condCell,1);
for k = 1:nCond
    index = strcmp(conditions(k), events);
    
    onsets = str2double(string(eventCell(index,1)));
    durations = str2double(string(eventCell(index,2)));

    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).name = char(conditions(k));
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).onset = onsets;
    matlabbatch{1}.spm.stats.fmri_spec.sess.cond(k).duration = durations;
end

% Contrasts
fid = fopen(fullfile(dataPath, contrastStr));
commandString = ['%s%s',repmat('%f',1,nCond)]; % Always start by reading two strings
contrastCell = textscan(fid,commandString);  
fclose(fid);

contrasts = cell2mat(contrastCell(:,3:2+nCond));
nContrasts = size(contrasts,1);
contrastNames = contrastCell{2};

contrastVectors = zeros(nContrasts,2*nCond);
contrastVectors(:,2*(1:nCond) - 1) = contrasts;
contrastVectors = contrastVectors ./ sum(abs(contrastVectors),2);

% Head motion parameters
hmpFile = fullfile(dirBOLD, ['rp_',fileStr,'bold.txt']);
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = cellstr(hmpFile);

% High-pass filter
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = 128;

matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 8;

% HRF with one temporal derivative
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [1,0];

% Masking threshold
if isfield(dS, 'mThresh')
    matlabbatch{1}.spm.stats.fmri_spec.mthresh = dS.mThresh;
end

% Timing units (seconds/scans)
if isfield(dS, 'timingUnits')
    matlabbatch{1}.spm.stats.fmri_spec.timing.units = dS.timingUnits;
end

% Explicit mask
if isfield(dS, 'explicitMask')
    matlabbatch{1}.spm.stats.fmri_spec.mask = dS.explicitMask;
end

% Estimate model
matlabbatch{2}.spm.stats.fmri_est.spmmat = cellstr(SPMFile);    
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.space.volume.block_type = 'Slices';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.signal = 'UGL';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.ARP = 3;
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.noise.UGL = 1;
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.LogEv = 'No';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.first = 'No';
matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.anova.second = 'Yes';
for k = 1:nContrasts
    contrastName = contrastNames{k};
    contrastVector = contrastVectors(k,:);
    matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.gcon(k).name = contrastName; 
    matlabbatch{2}.spm.stats.fmri_est.method.Bayesian.gcon(k).convec = contrastVector;
end
spm_jobman('run',matlabbatch);
