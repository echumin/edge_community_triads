function [tri_counter,trifrac_by_comm,cid_by_comm,tri_by_hemi,ntriad,cedge_by_comm] = edge_triads(input_mat,hemiID,ref_nodes,touch_nodes)
%% --- EDGE TRIADS ---%%
% Using an edge community affiliation node-by-node matrix (mat) count the
% proportion of each type of edge community triads.

% INPUTS
%   input_mat    - a node-by-node-by-1 edge community affiliation matrix
%                - a 3D matrix where the third dimention is either null
%                  permutations or subjects will also be accepted.
%   hemiID       - an index of hemispheres (1 = left ; 2 = right)
%   ref_nodes    - [optional] a vector of node indices for which to 
%                  count triads.
%   touch_nodes  - [optional; requires ref_nodes] a vector of node indices
%                  which the ref nodes must touch in the triads.
%
% IN THE CURRENT IMPLEMENTATION PROVIDING REFERENCE NODES AS INPUT WILL
% EXCLUDE TRIADS WITH >1 REFERENCE NODE. ADDITIONALLY, PROVIDING TOUCH
% NODES WILL EXCLUDE TRIADS WITH <2 OR >2 TOUCH NODES.
% IDEA IS TO LOOK AT TRIADS WITH SINGLE SUBCORTICAL (REF) AND TWO CORTICAL
% (TOUCH) NODES.

% OUTPUTS
%   tri_counter     - a node-by-4 triad type matrix with proportions (out of
%                   total nubmer of possible triads around a node) for each 
%                   node and each triad type.
%   trifrac_by_comm - a node-by-# of communities-by-3 triad type matrix
%                   (excluding diverse). A sorting of tri_counder by edge
%                   communities.
%   cid_by_comm     - 3 (triad-type) element cell with ref_nodes on the
%                   rows and counts of each touch node appearing in triads 
%                   with that ref node, separated by edge community on the 
%                   3rd dimention.  
%   tri_by_hemi
%   ntriad
%   cedge_by_comm
%
% TRIAD TYPES
%   *denotes edges that touch the 'reference' (R) node around which the 
%   triads are being counted.
%
%   1 - A*:A*:A   - a single edge community (CLOSED LOOP)
%   2 - A*:A*:B   - two edge communities, with single community touching 
%                   the R node (FORKED)
%   3 - B*:A*:A   - two edge communities, with two communties touching
%                   the R node (L-SHAPE)
%   4 - A*:B*:C   - three edge communties (DIVERSE)


% 2021 - Evgeny Chumin, IUB: Original Version 
%%
if size(input_mat,1)~=size(input_mat,2)
    fprintf(2,'Input matrix is not square.\n')
    return
end
if size(input_mat,3)>1
    fprintf('Input mat is 3D. Assuming null model matrices are being input.\n')
    Np = size(input_mat,3);
else
    Np=1;
end
if length(hemiID)~=size(input_mat,1)
    fprintf(2,'Hemisphere indiced do not match number of nodes.\n')
    return
elseif length(unique(hemiID)) > 2
    fprintf(2,'More than 2 hemisphere indices present.')
    return
end
% build matrix of within/between hemisphere edges
% 0=between hemispheres 1=within hemisphere
for i=1:length(hemiID)
    for j=1:length(hemiID)
        if hemiID(i)==hemiID(j)
            hemiMAT(i,j)=1;
        else
            hemiMAT(i,j)=0;
        end
    end
end

%% build indices of all triangles
% all unique node triad labels
%{ 
    pq is an index of all possible node triads. This index takes a while
    to build (~90min for ~250 nodes), so its recommended to prebuild this
    index and save it. 
    This will save a great deal of time when running multiple matrices.
    Example for building the index:
        tri = nchoosek(1:1:N,3);
        pq=double.empty;
        for ie=1:length(tri)
            pq((ie*3)-2:ie+(ie*2),:)= [tri(ie,1) tri(ie,2);tri(ie,1) tri(ie,3);tri(ie,2) tri(ie,3)];
        end
        save(['pq' num2str(N) '.mat'],'pq')
%}
load(['pq' num2str(size(input_mat,1)) '.mat'])

if exist('ref_nodes','var')
    fprintf('Counting triads for nodes in index\n')
    N = length(ref_nodes);
    ref_idx=sum(pq(:,1)==ref_nodes,2)+sum(pq(:,2)==ref_nodes,2);
    % remove connection among ref nodes
    ref_idx(ref_idx==2)=nan;
    % index of triads with a single ref node
    ref_idx=squeeze(sum(permute(reshape(ref_idx',[1,3,length(ref_idx)/3]),[2 1 3])))==2;

else
    fprintf('Counting triads for all nodes in the network.\n')
    N = size(input_mat,1);
    ref_nodes = 1:1:N;
end
if exist('touch_nodes','var')
    N2 = length(touch_nodes);
    touch_idx=sum(pq(:,1)==touch_nodes,2)+sum(pq(:,2)==touch_nodes,2);
    % remove points where only one or no touch nodes are present
    touch_idx(touch_idx==0)=nan;
    % index of triads with 2 touch nodes
    touch_idx=squeeze(sum(permute(reshape(touch_idx',[1,3,length(touch_idx)/3]),[2 1 3])))==4;
else
    N2 = size(input_mat,1);
    touch_nodes = 1:1:N2;
end 
nc = max(input_mat(:));

% preallocate
tri_counter = zeros(N,4,Np);       % 4 triad types
trifrac_by_comm = zeros(N,nc,3,Np);    % excluding diverse triads 
cid_by_comm = cell(3,Np);          % excluding diverse triads
tri_by_hemi = cell(3,Np);            % excluding diverse triads
cedge_by_comm = cell(3,Np);
    for np=1:Np
        for i =1:3
            cid_by_comm{i,np} = zeros(N,N2,nc);
            tri_by_hemi{i,np} = zeros(N,3,nc); % ipsilateral, bilateral, or contralateral triad relative to reference node
            cedge_by_comm{i,np} = zeros(N2,N2,N,nc);
        end
    end
    
%%
pq2=permute(reshape(pq',[2,3,length(pq)/3]),[2 1 3]);
if exist('ref_idx','var')
    if exist('touch_idx','var')
        pq2=pq2(:,:,(ref_idx & touch_idx));
    else
        pq2=pq2(:,:,ref_idx);
    end
end
pq = reshape(permute(pq2,[2 1 3]),[2,length(pq2)*3])';
        
if Np>1
    fprintf('Percent of mats completed...0')
end
for p=1:Np
    mat = input_mat(:,:,p);
%---------------------------%
% indices of edges in matrix
idx = sub2ind(size(mat),pq(:,1),pq(:,2));
% pull out community labels
vec = mat(idx);
vec2 = hemiMAT(idx);
% reshape into each column a triad
ecom=reshape(vec,[3,length(vec)/3]);
hemiEDGE=reshape(vec2,[3,length(vec2)/3]);
%----------------------------%
% counting number of communities in triads
counts = arrayfun(@(C) numel(unique(ecom(:,C))), 1:size(ecom,2));
idx_one = counts'==1;
idx_two = counts'==2;
idx_three = counts'==3;

for n=1:N
    % find indices of trids that contain the reference node
    idx2=squeeze(sum(sum(pq2==ref_nodes(n))))>0;
    nt = nnz(idx2); % number of triads
    ntriad(n,1)=nt;
    %% closed loops
    % Triads with reference node that are also closed loops
    idx2c = idx_one & idx2;
    cm = ecom(1,idx2c); % edge communities of closed loops
    % fraction of total triads around a node that are closed loops 
    tri_counter(n,1,p) = length(cm)./nt;
    % fraction of total triads around a node that are closed loops
    % by edge community
    cl_com = cell2mat(arrayfun(@(x)length(find(cm == x)), 1:1:nc, 'Uniform', false));
    trifrac_by_comm(n,:,1,p) = cl_com./nt;
   
    hm = hemiEDGE(:,idx2c); %hemisphere info of closed loop triads
   
    % Identify touch_nodes
    corners=pq2(:,:,idx2c); % node indices of closed loop triads with ref_node
    idxe = find(cl_com>0);  % find edge communities with nonzero triad counts
    % for each of those edge communties
    for c=1:length(idxe) 
        tmp = corners(:,:,cm==idxe(c)); % node indices of triads in edge community idxe(c)
        tmp2 = hm(:,cm==idxe(c));       % hemisphere indices of edges in edge community idxe(c)
        idxr = squeeze(sum(tmp==ref_nodes(n),2)); % reference node positions in triad
        % hemisphere indices of edges from the reference node
        re = reshape(tmp2(logical(idxr)),2,size(idxr,2)); 
        % count the number of times touch nodes apprear in triads
        cid_by_comm{1,p}(n,:,idxe(c))=cell2mat(arrayfun(@(x)sum(tmp(:) == x)/2, 1:1:N2, 'Uniform', false));
        % count triad laterality relative to reference node
        tri_by_hemi{1,p}(n,1,idxe(c))= sum(sum(re)==2); % ipsilateral
        tri_by_hemi{1,p}(n,2,idxe(c))= sum(sum(re)==1); % bilateral
        tri_by_hemi{1,p}(n,3,idxe(c))= sum(sum(re)==0); % contralateral
        if size(tmp,3)>1
            tmp_v = reshape(permute(tmp,[2 1 3]),[2,size(tmp,3)*3])';
        else
            tmp_v = tmp;
        end
        idxr_v = idxr(:);
        crt_edges = tmp_v(idxr_v==0,:);
        crt_ind = sub2ind([N2 N2],crt_edges(:,1),crt_edges(:,2));
        crtmt = zeros(N2,N2);
        crtmt(crt_ind)=1;
        cedge_by_comm{1,p}(:,:,n,idxe(c))=crtmt;
        clear tmp_v idxr_v crtmt tmp tmp2 idxr re
    end
    clear cm corners c idxe hm cl_com 
    %% forked and l-shape
    % Triads with reference node that are either forked or L
    idxfl=idx_two & idx2;
    cm=ecom(:,idxfl);       % edge community labels
    corners=pq2(:,:,idxfl); % node indices of triads
    idx5=logical.empty;
    for f=1:size(cm,2)
        [r,~]=find(corners(:,:,f)==ref_nodes(n));
        idx5(f,1)=cm(r(1),f)==cm(r(2),f); % forked = 1; l-shape = 0    
    end
    
    % fractions of triads
    tri_counter(n,2,p)=sum(idx5)./nt;   % forked triads
    tri_counter(n,3,p)=sum(idx5==0)./nt;% L triads
    
    hm = hemiEDGE(:,idxfl); %hemisphere info of triads 
    %
    corners_frk = corners(:,:,idx5);
    hm_frk = hm(:,idx5);
    cmf = cm(:,idx5);  % edge communities of forked triads
    cmf = mode(cmf,1); % dominant (fork) community
    cf_com = cell2mat(arrayfun(@(x)sum(cmf == x), 1:1:nc, 'Uniform', false));
    trifrac_by_comm(n,:,2,p) = cf_com./nt;
    %
    idxe1 = find(cf_com>0); % find edge communities with nonzero triad counts
    for c=1:length(idxe1)
        tmp = corners_frk(:,:,cmf==idxe1(c));  % node indices of triads in edge community idxe(c)
        tmp2 = hm_frk(:,cmf==idxe1(c));        % hemisphere indices of edges in edge community idxe(c)
        idxr = squeeze(sum(tmp==ref_nodes(n),2)); % reference node positions in triad
        cid_by_comm{2,p}(n,:,idxe1(c))=cell2mat(arrayfun(@(x)sum(tmp(:) == x)/2, 1:1:N2, 'Uniform', false));
        % hemisphere indices of edges from the reference node
        re = reshape(tmp2(logical(idxr)),2,size(idxr,2)); 
        % count triad laterality relative to reference node
        tri_by_hemi{2,p}(n,1,idxe1(c))= sum(sum(re)==2); % ipsilateral
        tri_by_hemi{2,p}(n,2,idxe1(c))= sum(sum(re)==1); % bilateral
        tri_by_hemi{2,p}(n,3,idxe1(c))= sum(sum(re)==0); % contralateral
        if size(tmp,3)>1
            tmp_v1 = reshape(permute(tmp,[2 1 3]),[2,size(tmp,3)*3])';
        else
            tmp_v1 = tmp;
        end
        idxr_v1 = idxr(:);
        crt_edges1 = tmp_v1(idxr_v1==0,:);
        crt_ind1 = sub2ind([N2 N2],crt_edges1(:,1),crt_edges1(:,2));
        crtmt1 = zeros(N2,N2);
        crtmt1(crt_ind1)=1;
        cedge_by_comm{2,p}(:,:,n,idxe1(c))=crtmt1;
        clear tmp_v1 idxr_v1 crtmt1 tmp tmp2 idxr re
    end
    
    corners_lsp = corners(:,:,idx5==0);
    hm_lsp = hm(:,idx5==0);
    cml = cm(:,idx5==0); % edge communities of L triads
    cml = mode(cml,1);      % dominant (L) community
    cls_com = cell2mat(arrayfun(@(x)sum(cml == x), 1:1:nc, 'Uniform', false));
    trifrac_by_comm(n,:,3,p) = cls_com./nt;
    %
    idxe2 = find(cls_com>0); % find edge communities with nonzero triad counts
    for c=1:length(idxe2)
        tmp = corners_lsp(:,:,cml==idxe2(c)); % node indices of triads in edge community idxe(c)
        tmp2 = hm_lsp(:,cml==idxe2(c));        % hemisphere indices of edges in edge community idxe(c)
        idxr2 = squeeze(sum(tmp==ref_nodes(n),2)); % reference node positions in triad
        cid_by_comm{3,p}(n,:,idxe2(c))=cell2mat(arrayfun(@(x)sum(tmp(:) == x)/2, 1:1:N2, 'Uniform', false));
        % hemisphere indices of edges from the reference node
        re2 = reshape(tmp2(logical(idxr2)),2,size(idxr2,2)); 
        % count triad laterality relative to reference node
        tri_by_hemi{3,p}(n,1,idxe2(c))= sum(sum(re2)==2); % ipsilateral
        tri_by_hemi{3,p}(n,2,idxe2(c))= sum(sum(re2)==1); % bilateral
        tri_by_hemi{3,p}(n,3,idxe2(c))= sum(sum(re2)==0); % contralateral
        if size(tmp,3)>1
            tmp_v2 = reshape(permute(tmp,[2 1 3]),[2,size(tmp,3)*3])';
        else
            tmp_v2 = tmp;
        end
        idxr_v2 = idxr2(:);
        crt_edges2 = tmp_v2(idxr_v2==0,:);
        crt_ind2 = sub2ind([N2 N2],crt_edges2(:,1),crt_edges2(:,2));
        crtmt2 = zeros(N2,N2);
        crtmt2(crt_ind2)=1;
        cedge_by_comm{3,p}(:,:,n,idxe2(c))=crtmt2;
        clear tmp_v2 idxr_v2 crtmt2 tmp tmp2 idxr2 re
    end
    
    clear cm corners c idxe1 idxe2 tmp tmp2 idxr idxr2 re re2 cls_com cf_com
    %% diverse
    idxd = idx_three & idx2;
    cm = ecom(1,idxd);
    tri_counter(n,4,p) = length(cm)./nt;
end
clear mat
    prg = (p/Np)*100;
    if rem(prg,10)==0
        if prg == 100
            fprintf('..%d\n Done!\n',prg)
        else
            fprintf('..%d',prg)
        end
    end
end



