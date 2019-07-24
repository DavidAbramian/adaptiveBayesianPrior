% Better simple model
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PURPOSE:      Setup precision matrices
%               Q = G'*G
%               G not properly updated for all dimensions/types
%               
% AUTHOR:       Per Siden
%               Division of Statistics and Machine Learning
%               Department of Computer and Information Science
%               Linkoping University      
%
% FIRST VER.:   2016-06-09
% REVISED:      
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [QList,GList] = spm_svb_setupPrecMats(QTypes,N,sz,bmask,ndim,sliceNbr,plotResult)

QTypes = strcat(QTypes,'adapt');

K = length(QTypes);
QList = cell(K,1);
GList = cell(K,1);

for k = 1:K
                
    if strcmp(QTypes{k},'LIadapt') % adaptive spatial model
             
flipControlTensor = 0;

tensorSmoothing = 5;
smoothTensor = 0;
filterSize = 3;

% Read T1 data
T1dataPath = '/home/davab27/Work/BayesianProject/Bayesian/HCP/sub-1/anat/sub-1_T1w.nii'; % HCP data
% T1dataPath = '/home/davab27/Bayesian/br-tum-bids/sub-01/anat/sub-01_T1w_3.9949.nii'; % Brain tumor data

[~, nifti] =  ml_load_nifti(T1dataPath);
nifti = double(nifti);
nifti = nifti / max(nifti(:));
noise_im = nifti(:,:,sliceNbr);

[sy, sx] = size(noise_im);
if length(bmask) == numel(noise_im)
    noIndex = 1;
else
    noIndex = 0;
end

%% Optimization of quadrature filters

% Center frequency for quadrature filters
u0 = pi/3;

% Bandwidth (in octaves) for quadrature filters
B = 2;

% Size of filter in spatial domain
spatial_size = [13 13];

% Size of filter in frequency domain,
% larger than spatial size to get a better estimate of the error
% filter is continuous in frequency domain
frequency_size = 2*spatial_size+1;

% Frequency weight function, punish errors for low frequencies
frequency_rexp = -1;
cosexp = 1;
SNRest = 10;

Fw = F_weightgensnr(frequency_size,frequency_rexp,cosexp,SNRest);
% figure; wamesh(log(Fw)); title('Log of frequency weight function')


% Spatial weight function, punish filter coefficients far from origin in spatial domain
spatial_rexp = 2;
fw0 = goodsw(spatial_size,spatial_rexp);
% figure; wamesh(fw0); title('Spatial weight function')

% Spatial weight scalar in optimization. Should we care more about error in spatial domain or frequency domain?
fw_amp = 300;
fw = fw_amp.*fw0;

% Spatial ideal filter, Dirac impulse
fi = wa(spatial_size,0);
fi = putorigo(fi,1);

% Spatial mask
fm = wa(spatial_size,0);
fm = putdata(fm,ones(spatial_size));

% Filter directions ( [y x] due to Matlab notation )
for otherK = 1:4
    dir{otherK} = [sind((otherK-1)*45)  cosd((otherK-1)*45) ]';
end

% Get the ideal quadrature filters in the frequency domain
for otherK = 1:4
    Fi{otherK} = quadrature(frequency_size,u0,B,dir{otherK});
%     figure; wamesh(Fi{k}); title(['Ideal quadrature filter in frequency domain for direction ' num2str(k)])
end

% Optimize the quadrature filters
for otherK = 1:4
    [f{otherK},F{otherK}] = krnopt(Fi{otherK},Fw,fm,fi,fw);
%     figure; wamesh(F{k}); title(['Optimized quadrature filter in frequency domain for direction ' num2str(k)])
end

% for otherK = 1:4
%     figure; wamesh(real(f{otherK})); title(['Real part of optimized quadrature filter in spatial domain for direction ' num2str(otherK)])
% end
% 
% for otherK = 1:4
%     figure; wamesh(imag(f{otherK})); title(['Imaginary part of optimized quadrature filter in spatial domain for direction ' num2str(otherK)])
% end

% filterr(F{1},Fi{1},Fw,fi,fw,1);

%% Apply quadrature filters to image and estimate tensor

% Do convolution
for otherK = 1:4
    q{otherK} = abs(imfilter(noise_im,getdata(f{otherK}),'conv','replicate'));
end

% figure; imagesc(q{1}); title('Magnitude of first quadrature filter response'); colormap gray

% Calculate projection tensors
for otherK = 1:4
    M{otherK} = 4/4 * dir{otherK} * dir{otherK}' - 1/4*eye(2);
end

% Calculate tensor components
t1 = zeros(size(noise_im));
t2 = zeros(size(noise_im));
t3 = zeros(size(noise_im));

for otherK = 1:4
    Mk = M{otherK};
    qk = q{otherK};
    t1 = t1 + qk * Mk(1,1);
    t2 = t2 + qk * Mk(1,2);
    t3 = t3 + qk * Mk(2,2);
end

% % Look at tensor components
% figure; mygopimage(t1,-1)
% figure; mygopimage(t2,-1)
% figure; mygopimage(t3,-1)

if smoothTensor
    % Lowpass filter the tensor components
    t1 = aver2(t1,tensorSmoothing);
    t2 = aver2(t2,tensorSmoothing);
    t3 = aver2(t3,tensorSmoothing);
end

%--------------------------------------
% Map and flip tensor
%--------------------------------------

% Calculate eigen values

lambda_1 = (t1 + t3)/2 + sqrt(t2.^2 + 1/4*t1.^2 + 1/4*t3.^2 - t1.*t3/2);
lambda_2 = (t1 + t3)/2 - sqrt(t2.^2 + 1/4*t1.^2 + 1/4*t3.^2 - t1.*t3/2);

% Look at eigen values

% figure; mygopimage(lambda_1,-1)
% figure; mygopimage(lambda_2,-1)

% Make sure second eigen value is positive
lambda_2 = abs(lambda_2);


% Map norm of tensor
x = sqrt(lambda_1.^2 + lambda_2.^2);
x = x / max(x(:));

% Parameters for m-func

% Depend on noise level!

m_sigma = 0.2; % noise threshold
% m_alpha = 0.2; % overshoot, to amplify weak structures
m_alpha = 0; % overshoot, to amplify weak structures
m_beta = 4; % steepness of curve
m_j = 1;

gamma_1 = m_func(x, m_sigma, m_alpha, m_beta, m_j, 0);

% figure; imagesc(gamma_1); title('Mapped tensor magnitude'); colormap gray

% Map anisotropy of tensor
x = lambda_2 ./(lambda_1 + eps);

% Parameters for mu-func
mu_alpha = 0.4;
mu_beta = 3;
mu_j = 1;

mu_2 = mu_func(x, mu_alpha, mu_beta, mu_j, 0);
gamma_2 = mu_2 .* gamma_1;

% Calculate outer product of eigen vectors
% not necessary to know eigen vectors
ee1_1 = t1 - lambda_2;
ee1_2 = t2;
ee1_3 = t3 - lambda_2;

% Calculate norm, instead of dividing
norm = sqrt(ee1_1.^2 + 2*ee1_2.^2 + ee1_3.^2);
norm = 1./(norm + eps);

% Normalize outer product
ee1_1 = norm .* ee1_1;
ee1_2 = norm .* ee1_2;
ee1_3 = norm .* ee1_3;

ee2_1 = 1 - ee1_1;
ee2_2 = 0 - ee1_2;
ee2_3 = 1 - ee1_3;

% Compute the control tensor
if flipControlTensor == 0
    c1 = gamma_1 .* ee1_1 + gamma_2 .* ee2_1;
    c2 = gamma_1 .* ee1_2 + gamma_2 .* ee2_2;
    c3 = gamma_1 .* ee1_3 + gamma_2 .* ee2_3;
else
    c1 = gamma_2 .* ee1_1 + gamma_1 .* ee2_1;
    c2 = gamma_2 .* ee1_2 + gamma_1 .* ee2_2;
    c3 = gamma_2 .* ee1_3 + gamma_1 .* ee2_3;
end

%% Find eigenvectors of all tensors

eigVals = zeros(2, sy, sx);
eigVecs = zeros(2, 2, sy, sx);

for y = 1:sy
    for x = 1:sx
%         tensor = [t1(y,x) t2(y,x); t2(y,x) t3(y,x)];
        tensor = [c1(y,x) c2(y,x); c2(y,x) c3(y,x)];
        [a,b] = eigs(tensor);
        eigVals(:,y,x) = diag(b);
        eigVecs(:,:,y,x) = a;
    end
end

% Correct 2nd eigenvalue
eigVals(2,:) = abs(eigVals(2,:));

%% Test filters on brain image

method = 'anyDir';

n = (filterSize-1)/2;
x = (-n:n) .* ones(2*n+1,1);
y = x';
phi = atan2(y,x);
phi(n+1,n+1) = NaN;
r = sqrt(x.^2 + y.^2);
        
% p1 = 8;
% p2 = 2;

% p1 = 20;
% p2 = 4;

p1 = 12;
p2 = 5;

fprintf('p1 = %i, p2 = %i \n', p1, p2)

filters = zeros(sy,sx,filterSize,filterSize);
% figure
for x = 1:sx
    for y = 1:sy
        % Calculate angle from first eigenvector
        ang = atan(eigVecs(1,1,y,x) / eigVecs(2,1,y,x));
        
        % Filter is cosine around this angle
%         filt = -cos(phi - ang).^20 ./ r.^2;
        filt = -cos(phi - ang).^p1 ./ r.^p2;
        const = -sum(filt(~isnan(filt)));
        filt(n+1,n+1) = const;
        
        % Normalize filter to have 2 at center
%         filters(y,x,:,:) = rot90(filt) * (2/const);
        filters(y,x,:,:) = rot90(filt);
        
%         tempIm = noise_im;
%         tempIm(y,x) = 1;
%         subplot(121)
%         imagesc(tempIm,[0,1]), axis image, colormap gray
%         
%         subplot(122)
%         imagesc(filt), axis image, colormap gray
%         
%         pause
        
    end
end


%% Calculate filter matrix

filterMatrix = zeros(sy*sx);

offset = (filterSize - 1)/2;

for i = 1:sx
    for j = 1:sy
        imPad = zeros(sy + filterSize -1, sx + filterSize -1);
        imPad( j:j+2*offset, i:i+2*offset ) = squeeze(filters(j,i,:,:)); 
        
        temp = imPad(offset+1:offset+sy, offset+1:offset+sx);

        filterMatrix(sy*(i-1)+j,:) = temp(:);
    end
end

if noIndex
    mat2 = sparse(filterMatrix);
else
    mat2 = sparse(filterMatrix(bmask,bmask));
end

%% Try to fix filter matrix (make each row sum to 0)


%% Filter image

imFiltered = zeros(size(noise_im));
imFiltered(bmask) = mat2 * noise_im(bmask);

a = issymmetric(mat2);
b = max(abs(sum(mat2,2)));

% figure
% % subplot(131)
% maxAbs = max(abs(imFiltered(:)));
% imagesc(imFiltered,[-maxAbs,maxAbs]), colorbar, colormap gray, axis image


% Make the matrix symmetric
mat3 = (mat2 + mat2')/2;

imFiltered2 = zeros(size(noise_im));
imFiltered2(bmask) = mat3 * noise_im(bmask);

% figure
% % subplot(132)
% maxAbs = max(abs(imFiltered2(:)));
% imagesc(imFiltered2,[-maxAbs,maxAbs]), colorbar, colormap gray, axis image


% Fix symmetric matrix
mat4 = mat3;
d = find(speye(size(mat4)));
fix = sum(mat4,2);
mat4(d) = mat4(d)-fix;

imFiltered3 = zeros(size(noise_im));
imFiltered3(bmask) = mat4*noise_im(bmask);

if plotResult
    figure
    % subplot(133)
    maxAbs = max(abs(imFiltered3(:)));
    imagesc(imFiltered3,[-maxAbs,maxAbs]), colorbar, colormap gray, axis image

    MSE = sqrt(mean(imFiltered3.^2,'all'));
    title(sprintf('%s, MSE = %f, p1 = %.1f, p2 = %.1f \n', method, MSE, p1, p2))
end

%% Create G matrix from Q

Q = mat4;

% Q = QList{1};

blah = triu(Q,1);

nEdges = nnz(blah);
I = find(blah);

[y,x] = ind2sub(size(Q), I);

xx = [(1:nEdges)', x];
yy = [(1:nEdges)', y];

G = zeros(nEdges, N);

ind1 = sub2ind(size(G), xx(:,1), xx(:,2));
ind2 = sub2ind(size(G), yy(:,1), yy(:,2));

G(ind1) = sqrt(-blah(I));
G(ind2) = -sqrt(-blah(I));

QList = repmat({sparse(Q)}, 1, K);
GList = repmat({sparse(G)}, 1, K);

break

    else
        fprintf('ERROR: Invalid Precision Matrix Type! \n');
    end
    
end


