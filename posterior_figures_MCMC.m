
maskFile = fullfile('Results-onlytop','sub-1','SVB2D_Per','mask.nii');
[~, maskVol] = ml_load_nifti(maskFile);

sliceNbr = 68;
mask = logical(maskVol(:,:,sliceNbr));

I = find(mask);
N = length(I);

[~, T1vol] = ml_load_nifti(fullfile('HCP','sub-1','anat','sub-1_T1w.nii'));

%% Generate precision matrix Dw

for setting = 1:3
    switch setting
        case 1
            dataName = '_Per';
            figString = 'UGL';
        case 2
            dataName = '_simple_model';
            figString = '4DIR';
        case 3
            dataName = '_better_simple_model';
            figString = 'ANYDIR';     
    end

    dataDir = fullfile('Results-onlytop','sub-1', ['MCMC2D' dataName], 'MCMC.mat');
    load(dataDir, 'MCMC')

    %%

    blah = nan(size(mask));
    blah(I) = MCMC.b(9).wPostMean(3,:);
%     blah(I) = MCMC.b(9).aPostMean(1,:);

    figure(1)
    subplot(2,3,setting)

    imagesc(blah,[-3.5,3.5])
    colormap gray
    colorbar
    title(strcat("Posterior mean of W'_{3,\cdot} with ", figString, " model"))
%     title(['Posterior mean of A_1 with ', figString, ' model'])

    
    blah(I) = MCMC.b(9).wPostStd(3,:);
%     blah(I) = MCMC.b(9).aPostStd(1,:);

    subplot(2,3,setting+3)

    imagesc(blah,[0,1])
    colormap gray
    colorbar
    title(strcat("Posterior std of W'_{3,\cdot} with ", figString, " model"))
%     title(['Posterior std of A_1 with ', figString, ' model'])
end

%%
set(gcf, 'Position', [0,0,2000,1200])

