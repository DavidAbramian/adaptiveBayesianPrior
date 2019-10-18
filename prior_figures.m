
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
            D = spm_svb_setupPrecMats_Per({'LI'}, [], size(mask), I, 2);
            figString = 'UGL';
        case 2
            D = spm_svb_setupPrecMats_simple_model({'LI'}, N, size(mask), I, 2, sliceNbr, 0);
            figString = '4DIR';
        case 3
            D = spm_svb_setupPrecMats_better_simple_model({'LI'}, N, size(mask), I, 2, sliceNbr, 0);
            figString = 'ANYDIR';     
    end

    D = D{1};
    
    %%
    [A,B] = eig(full(D));

    nEigs = nnz(diag(B) > 1e-10);
    a = length(B) - (nEigs-1);

    sqrtMat = real(A(:,a:end) * B(a:end,a:end)^(1/2));

    %%
    
    nSamples = 10000;

    q1 = 10;
    q2 = 0.1;

    alpha = gamrnd(q2,q1,1,nSamples);

    y = randn(nEigs,nSamples);
    z = sqrtMat*y;

    z2 = z ./ sqrt(alpha);

    priorMean = mean(z2,2);
    priorStd = var(z2,[],2);

    %%

    blah = nan(size(mask));
    blah(I) = priorMean;

    figure(1)
    subplot(2,3,setting)

    imagesc(blah)
    colormap gray
    colorbar
    title(strcat("Prior mean of W'_{k,\cdot} with ", figString, " model"))
    
    blah(I) = priorStd;

    subplot(2,3,setting+3)

    imagesc(blah)
    colormap gray
    colorbar
    title(strcat("Prior std of W'_{k,\cdot} with ", figString, " model"))

end

%%
set(gcf, 'Position', [0,0,2000,1200])

