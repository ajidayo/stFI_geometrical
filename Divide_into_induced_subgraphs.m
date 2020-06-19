function [subG_bin,subG_sizes,allIdx_stFI,UpdateNum] ...
    = Divide_into_induced_subgraphs(sC,UpdateNum,MeshNum,att)

global DISPDEBUGGINGMESSAGE

if DISPDEBUGGINGMESSAGE
    disp('Divide_into_induced_subgraphs:CALLED')
end
% this unit divides the spatial mesh into subgraphs 
% after omitting the elements nearby the timestep-varing boundaries
% within the subgraphs the timemarching can be written in the spatial FI formulation



%% omit the elements nearby the timestep-varing boundaries
% the edges surrounding the boundary faces is omitted

%case1(defalt)
%sComit=logical(sC); 
%case2(for matlab ver R2018a?)
sComit=abs(sC); 


%disp('checkpoint alpha')

allIdx_stFI.f=find(att.f.inc_bounde==true);
for ff=1:size(allIdx_stFI.f,1)
    f=allIdx_stFI.f(ff);
    %case1(defalt)
    %sComit(f,:)=false;
    %case2(for matlab ver R2018a?)
    sComit(f,:)=0;
end

    
%% generate graph object (not digraph)

%AdjF_omit = logical(sComit * (sComit.'));
%case1(defalt)
%G_f_omit = graph(sComit * sComit.','omitselfloops');
%case2(for matlab ver R2018a?)
G_f_omit = graph(sComit * sComit.','omitselfloops');


%% divide the spatial mesh into subgraphs
% conncomp reterns the connecting components of the graph
% subG_bin.f is the bin of f: the number shows the index of the subgraph
% that f is included.
% subG_sizes.f shows the number of faces included in each subgraph.
[subG_bin.f,subG_sizes.f] = conncomp(G_f_omit);
clearvars G_f_omit

[subG_bin.f,subG_sizes.f] = count_FItask(subG_bin.f,subG_sizes.f,att.f.inc_bounde);

subG_bin.e = zeros(MeshNum.E,1);
subG_sizes.e = zeros(size(subG_sizes.f,2),1);

%disp('checkpoint echo')



for e_tar=1:MeshNum.E
    if att.e.adj_bounde(e_tar)==true || att.e.boundaryedge(e_tar)==true
        continue;
    end
    
    % find a face incident to e_target and name it as "f_target"
    sComit_e_tar=sComit(:,e_tar);
    f_tar = find(sComit_e_tar,1);
    % find the index of the subgraph which includes f_target 
    subGIdx = subG_bin.f(f_tar);
    % the bin of e_target is identical to the bin of f_target
    subG_bin.e(e_tar) = subGIdx;
    % count the number of the edges included in the corresponding subdomain
    subG_sizes.e(subGIdx) = subG_sizes.e(subGIdx)+1;
end

%disp('checkpoint foxtrot')

%% limit the face-edge incidence matrices for the subgraphs
% N.B. in 2-D space case a graph and a face-edge incidence matrix is convertible 

% sC_SG=sparse(MeshNum.F,MeshNum.E,size(subG_sizes.f,1));
% sG_SG=sparse(MeshNum.E,MeshNum.N,size(subG_sizes.f,1));
% tildesG_SG=sparse(MeshNum.E,MeshNum.F,size(subG_sizes.f,1));

UpdateNum.SG=zeros(size(subG_sizes.f,2),1);
for subGIdx=1:size(subG_sizes.f,2)
     UpdateNum.SG(subGIdx)=UpdateNum.f(find(subG_bin.f==subGIdx,1));
end

allIdx_bound_e=find(att.e.boundaryedge==true);
%[size(allIdx_bound_e,1) size(allIdx_bound_e,2)]
allIdx_adjbounde_e=find(att.e.adj_bounde==true);
%[size(allIdx_adjbounde_e,1) size(allIdx_adjbounde_e,2)]
allIdx_stFI.e=[allIdx_bound_e;allIdx_adjbounde_e];
%[size(allIdx_stFI_e,1) size(allIdx_stFI_e,2)]
clearvars allIdx_bound_e allIdx_adjbounde_e

%disp('checkpoint golf')
if DISPDEBUGGINGMESSAGE
    disp('Divide_into_induced_subgraphs:ENDED')
end
end
