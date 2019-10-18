
maskFile = fullfile('Results-onlytop','sub-1','SVB2D_Per','mask.nii');
[~, maskVol] = ml_load_nifti(maskFile);

sliceNbr = 68;
mask = logical(maskVol(:,:,sliceNbr));

I = find(mask);
N = length(I);

%% Generate precision matrix Dw

D = spm_svb_setupPrecMats_Per({'LI'}, [], size(mask), I, 2);
% D = spm_svb_setupPrecMats_simple_model({'LI'}, N, size(mask), I, 2, sliceNbr, 0);
% D = spm_svb_setupPrecMats_better_simple_model({'LI'}, N, size(mask), I, 2, sliceNbr, 0);
D = D{1};

% Cholesky doesn't seem to work when matrices are non-singular

[A,B] = eig(full(D));

nEigs = nnz(diag(B) > 1e-10);
a = length(B) - (nEigs-1);

sqrtMat = real(A(:,a:end) * B(a:end,a:end)^(1/2));

%% Sample from the distribution

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

blah = zeros(size(mask));
blah(I) = priorMean;

figure(1)
imagesc(blah)
colormap gray
colorbar
title('Prior mean of W with UGL model')
set(gcf, 'Position', [0,0,600,500])
% axis equal

blah(I) = priorStd;

figure(2)
imagesc(blah)
colormap gray
colorbar
title('Prior std of W with UGL model')
set(gcf, 'Position', [0,0,600,500])
% axis equal

