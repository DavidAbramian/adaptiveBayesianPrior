% Per
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
%
% Modified:     2019-10 David Abramian
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function [QList,GList] = spm_svb_setupPrecMats(QTypes,N,sz,bmask,ndim,sliceNbr,plotResult)

K = length(QTypes);
QList = cell(K,1);
GList = cell(K,1);

for k = 1:K
    if strcmp(QTypes{k},'LI') % intrinsic Laplacian
        if ndim == 2
            Dr = spdiags([-1*ones(sz(1)-1,1), 1*ones(sz(1)-1,1)], [0,1], sz(1)-1,sz(1));
            Dc =  spdiags([-1*ones(sz(2)-1,1), 1*ones(sz(2)-1,1)], [0,1], sz(2)-1,sz(2));
            Drs = kron(speye(sz(2)), Dr);
            Dcs = kron(Dc, speye(sz(1)));
            
            bmask_mat = zeros(sz); 
            bmask_mat(bmask) = 1;
            
            hasColNeigh = bmask_mat(1:(end-1),:) & ~diff(bmask_mat);
            hasRowNeigh = bmask_mat(:,1:(end-1)) & ~diff(bmask_mat')';
            
            GList{k} = [Drs(hasColNeigh(:), bmask); Dcs(hasRowNeigh(:), bmask)];
            QList{k} = GList{k}'*GList{k};
        elseif ndim == 3
            Dx = spdiags([-1*ones(sz(1)-1,1), 1*ones(sz(1)-1,1)], [0,1], sz(1)-1,sz(1));
            Dy = spdiags([-1*ones(sz(2)-1,1), 1*ones(sz(2)-1,1)], [0,1], sz(2)-1,sz(2));
            Dz = spdiags([-1*ones(sz(3)-1,1), 1*ones(sz(3)-1,1)], [0,1], sz(3)-1,sz(3));
            Dxs = kron(speye(sz(3)), kron(speye(sz(2)),Dx));
            Dys = kron(speye(sz(3)), kron(Dy,speye(sz(1))));
            Dzs = kron(Dz, kron(speye(sz(2)),speye(sz(1))));
            
            bmask_mat = zeros(sz); 
            bmask_mat(bmask) = 1;
            
            hasXNeigh = bmask_mat(1:(end-1),:,:) & ~diff(bmask_mat);
            hasYNeigh = bmask_mat(:,1:(end-1),:) & ~diff(bmask_mat,1,2);
            hasZNeigh = bmask_mat(:,:,1:(end-1)) & ~diff(bmask_mat,1,3);
            
            GList{k} = [Dxs(hasXNeigh(:),bmask);Dys(hasYNeigh(:),bmask);Dzs(hasZNeigh(:),bmask)];
            QList{k} = GList{k}'*GList{k};
        end            
    else
        fprintf('ERROR: Invalid Precision Matrix Type! \n');
    end
    
end


