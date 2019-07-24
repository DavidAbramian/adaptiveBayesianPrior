
%% Settings

dS.dataPath = 'HCP';
dS.outputPath = 'Results-onlytop';

dS.explicitMask = {'/home/davab27/Work/BayesianProject/Bayesian2/HCP/sub-1/anat/sub-1_T1w.nii,1'};
% dS.mThresh = -Inf;
dS.mThresh = 0.2;

dS.subjStr = '1';
dS.taskStr = 'motor';
dS.runStr = '_run-01_';
dS.fileStr = ['sub-',dS.subjStr,'_task-',dS.taskStr,dS.runStr];
dS.TR = 0.72;
dS.condStr = 'condition_key.txt';
dS.contrastStr = 'task_contrasts.txt';