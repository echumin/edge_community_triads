% Evgeny J Chumin, Indiana University, 2022


% load edge community matrix and hemisphere indices.
% This is a 132 node matrix with 100 cortical Schaefer nodes
% and 32 subcortical Tian Scale 2 nodes.
load('edge_comm_mat.mat')

% subcortical nodes 101 through 132 are used as reference nodes, computing
% triads around these nodes
reference_nodes = 101:1:132;

% cortical nodes 1 through 100 are used as touch nodes, such that only
% triads that connect a single reference node to two touch nodes are
% counted.
touch_nodes = 1:1:100;

[tri_counter,trifrac_by_comm,cid_by_comm,tri_by_hemi,ntriad,cedge_by_comm] = ...
    edge_triads(edge_comm_mat,hemiID,reference_nodes,touch_nodes);

% OUTPUTS:

%% tri_counter
% for each reference node (row) and the four triad types (colum) contains
% the fraction of triads of that type. Rows will sum to 1 that is all
% triads around a node fall into one of four categories.
figure; imagesc(tri_counter)
xticks(1:1:4)
xticklabels({'closed loop','forked','L-shape','diverse'})
xlabel('Triad Type')
ylabel('Reference Nodes')
title('Triad Distributions by Type')

%% trifrac_by_comm
% Excluding the diverse triad type, which has no dominant community, the
% other three triad type distributions are broken up by which edge
% community has the most edges.
% This is a 3D matrix of nodes(rows) by community label (colums), by triad
% type (3rd dimention), ordered as closed-loop, forked, and l-shape.

%Example community 5:
figure; imagesc(squeeze(trifrac_by_comm(:,5,:)))
xticks(1:1:3)
xticklabels({'closed loop','forked','L-shape'})
xlabel('Triad Type')
ylabel('Reference Nodes')
title({'Fraction of Total Triads within edge community 5','broken down by type'})

%% cid_by_comm
% A 3 by 1 cell (closed loop, forked, L-shape, in that order), 
% where each triad type (cell element) contains a 3D matrix of refrence
% nodes (rows), by touch nodes (colums), by edge communities (3rd dimention).
% The values are the number of times that a particular touch node was in a 
% triad with a particular reference node. 

%Example cortical enpoints in community 5 forked triads:
figure; imagesc(cid_by_comm{2}(:,:,5))
xlabel('Schaefer 100 nodes (touch nodes)')
ylabel('Tian II nodes (reference nodes)')
title({'Counts of reference-touch pairs','within edge community 5 forked triads'})

%% tri_by_hemi
% A 3 by 1 cell (closed loop, forked, L-shape, in that order), 
% where each triad type (cell element) contains a 3D matrix of refrence
% nodes (rows), by triad lateralization type (columns), by edge communities
% (3rd dimention).
% Lateralization Types: 1 - ipsilateral (touch nodes in the same hemisphere
%                           as reference node).
%                       2 - bilateral (the two touch nodes are in
%                           differenct hemispheres).
%                       3 - contralateral (touch nodes are in the opposite
%                           hemisphere relative to the reference node).
% The values are counts of triads that fell within a particular type (cell
% element), a particular lateralization (columns in matrix), for a
% particular edge community (3rd dimention in matrix).

%Example of triad laterality in community 5 forked triads:
figure; plot(tri_by_hemi{2}(:,:,5)')
xlim([.5 3.5])
xticks(1:1:3)
xticklabels({'ipsilateral','bilateral','contralateral'})
xlabel('Triad Hemisphere Localizaiton')
ylabel('Triad Counts')
title({'Counts of Forked Triads in Edge Community 5',...
    'by Lateralization',...
    '(each line is a Tian II reference node)'})

%% ntriad
% A vector containing the total number of possible triads for each reference
% node,for a given reference and touch node set of inputs.
% Useful for potential scaling/normalization of other outputs.

%% cedge_by_comm
% A set of binary matrices showing the cortical edges for a particular triad
% type, of a particular edge community. (one for each reference node).

% A 3 by 1 cell (closed loop, forked, L-shape, in that order), 
% where each triad type (cell element) contains a 4D matrix, of touch-by-touch 
% node binary matrices (row, column), for each reference node (3rd dimention),
% for each edge community (4th dimention).

% Example of Cortical edges of Forked Triads in Edge Community 5 that
% connect to the 10th reference node.

binMat = cedge_by_comm{2}(:,:,10,5);
figure; spy(binMat)
xlabel('Touch Nodes (Schaefer 100)')
ylabel('Touch Nodes (Schaefer 100)')
title({'Touch Node Edges (Schaefer 100) of',...
        'Forked Triads in Edge Community 5',...
            'Connected to the 10th Reference Node (Tian II node 10)'})






























