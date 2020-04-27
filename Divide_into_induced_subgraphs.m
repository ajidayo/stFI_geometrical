function [subG_bin,subG_sizes,allIdx_stFI,denominator] ...
    = Divide_into_induced_subgraphs(sC,denominator,MeshNum,att)
disp('Divide_into_induced_subgraphs:CALLED')

%global EPSILON

denominator_face=denominator.f;

%%% input sC, sG, denominator (ADD att.bound_e ??? ; calculated previously in Constitutive)
% output denominator, subG_f_bin, subG_f_sizes, subG_e_bin, subG_e_sizes,
% allIdx_stFI_f, allIdx_stFI_e


% this unit divides the spatial mesh into subgraphs 
% by omitting the elements nearby the timestep-varing boundaries
% within the subgraphs the timemarching can be written in the spatial FI formulation

% inputs: timestep-varing boundary information, incedence matrix sC
% outputs: (small) incedence matrices sml_sC for the subgraphs


%% omit the elements nearby the timestep-varing boundaries
% the edges surrounding the boundary faces is omitted

%sComit=logical(sC);
sComit=logical(sC); 

% make adjacent matrix
%%%%%%% may be used in Constitutive; DO check later on %%%%%%%%
%AdjE=logical(sG*sG.');

%disp('checkpoint alpha')

% for e_target=1:MeshNum.E    
%     %%%%%% CAPSULIZE FROM %%%%%%
%     sC_e_target=sC(:,e_target);
%     row_sC_e_target=find(sC_e_target);
%     % check boundary-attribute
%     check_bound=0.0;
%     for ff=1:size(row_sC_e_target,1)
%         f=row_sC_e_target(ff);
%         check_bound=check_bound+sC(f,e_target)*denominator_face(f);
%     end
%     if abs(check_bound)<EPSILON % not a boundary
%          continue;
%     end
%     %%%%%% CAPSULIZE TO %%%%%%    
%     
%     %disp(e_target)
%     
%     % if e_target is a dt-varying boundary
%     att.bound_e(e_target)=true;
%     % for all f incident to e_target
%     for ff=1:size(row_sC_e_target,1)
%         f=row_sC_e_target(ff);
%         att.inc_bound_f(f)=true;
%         % for all e incident to f, detach from other faces
%         sComit(f,:)=false;
%     end
% end % for e_target
% att_bound_e=att.bound_e;

allIdx_stFI.f=find(att.inc_bound_f==true);
for ff=1:size(allIdx_stFI.f,1)
    f=allIdx_stFI.f(ff);
    sComit(f,:)=false;
end

% allIdx_bound_e=find(att_bound_e==true);
% for ee=1:size(allIdx_bound_e,1)
%     e_bound=allIdx_bound_e(ee);
%     AdjE_e_bound=AdjE(e_bound,:);
%     col_AdjE_e_bound=find(AdjE_e_bound);
%     for eee=1:size(col_AdjE_e_bound,2)
%         e_target=col_AdjE_e_bound(eee);
%         if att.bound_e(e_target)==false
%             att.adj_bounde_e(e_target)=true;
%         end
%     end
% end % for e_target
% 
% clearvars AdjE

%disp('checkpoint bravo')

% % att_inc_bound_f=att.inc_bound_f;
% allIdx_stFI.f=find(att.inc_bound_f==true);
% allIdx_stFI_f=allIdx_stFI.f;
% 
% for ff=1:size(allIdx_stFI_f,1)
%     %disp(ff)
%     f_target=allIdx_stFI_f(ff);
%     sC_f_target=sC(f_target,:);
%     col_sC_f_target=find(sC_f_target);
%     for ee=1:size(col_sC_f_target,2)
%         e=col_sC_f_target(ee);
%         if att.bound_e(e)==false && att.adj_bounde_e(e)==false
%             %disp(e)
%             att_bound_SG_e(e)=true;
%         end
%     end
% end
% 
% att.bound_SG_e=att_bound_SG_e;

%disp('checkpoint charlie')

    
%% generate graph object (not digraph)

%AdjF_omit = logical(sComit * (sComit.'));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G_f_omit = graph(sComit * sComit.','omitselfloops'); % OK ??? •Ó‚ªƒ_ƒu‚ç‚È‚¢H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%disp('checkpoint delta')

%% divide the spatial mesh into subgraphs
% conncomp reterns the connecting components of the graph
% subG_f_bin is the bin of f: the number shows the index of the subgraph
% that f is included.
% subG_f_sizes shows the number of faces included in each subgraph.
[subG_f_bin,subG_f_sizes] = conncomp(G_f_omit);
clearvars G_f_omit

[subG_f_bin,subG_f_sizes] = count_FItask(subG_f_bin,subG_f_sizes,att.inc_bound_f);

subG_e_bin = zeros(MeshNum.E,1);
subG_e_sizes = zeros(size(subG_f_sizes,2),1);

%disp('checkpoint echo')



for e_target=1:MeshNum.E
    if att.adj_bounde_e(e_target)==true || att.bound_e(e_target)==true
        continue;
    end
    
    % find a face incident to e_target and name it as "f_target"
    sComit_e_target=sComit(:,e_target);
    f_target = find(sComit_e_target,1);
    % find the index of the subgraph which includes f_target 
    subGIdx = subG_f_bin(f_target);
    % the bin of e_target is identical to the bin of f_target
    subG_e_bin(e_target) = subGIdx;
    % count the number of the edges included in the corresponding subdomain
    subG_e_sizes(subGIdx) = subG_e_sizes(subGIdx)+1;
end

%disp('checkpoint foxtrot')

%% limit the face-edge incidence matrices for the subgraphs
% N.B. in 2-D space case a graph and a face-edge incidence matrix is convertible 

% sC_SG=sparse(MeshNum.F,MeshNum.E,size(subG_f_sizes,1));
% sG_SG=sparse(MeshNum.E,MeshNum.N,size(subG_f_sizes,1));
% tildesG_SG=sparse(MeshNum.E,MeshNum.F,size(subG_f_sizes,1));

denominator_SG=zeros(size(subG_f_sizes,2),1);
for subGIdx=1:size(subG_f_sizes,2)
    % make the list of f, e whose bin is equal to subGIdx
%    f_lst_subG = find(subG_f_bin==subGIdx);
%    e_lst_subG = find(subG_e_bin==subGIdx);
%    denominator_SG(subGIdx)=denominator_edge(e_lst_subG(1));
     denominator_SG(subGIdx)=denominator_face(find(subG_f_bin==subGIdx,1));
    
%     % generate sC for subdomain subGIdx
%     LimtoSG_sC = sparse(MeshNum.F);
%     for ff=1:subG_f_sizes(subGIdx)
%         f=f_lst_subG(ff);
%         LimtoSG_sC(f,f)=1;
%     end   
%     sC_SG(:,:,subGIdx) = LimtoSG_sC * sC;
%     
    %   limtoSG_f= subG_f_bin==subGIdx;
    %   sC_SG(:,:,subGIdx) = diag(limtoSG_f)*sC
    
    % N.B. sC_subG includes the outer-boundary-edges of the subdomain.
    % Also, sC_subG does not include the row-entry for the faces which is 
    % incident to the dt-varying boundary edges 
    
    % generate sG for subgraph 
    % NOT neccesary if the st-grid corresponds to the FI method
    %     LimtoSG_sG = sparse(MeshNum.E);
    %     for ee=1:subG_e_sizes(subGIdx)
    %         e=e_lst_subG(ee);
    %         LimtoSG_sG(e,e)=1;
    %     end
    %     sG_SG(:,:,subGIdx) = LimtoSG_sG * sG;
    % N.B. sG_SG includes the outer-boundary-edges of the subdomain.
    
    % generate tildesG for subgraph
%     LimtoSG_sG = sparse(MeshNum.E);
%     for ee=1:subG_e_sizes(subGIdx)
%         e=e_lst_subG(ee);
%         LimtoSG_sG(e,e)=1;
%     end
%     tildesG_SG(:,:,subGIdx) = LimtoSG_sG * sC.';
    % N.B. tildesG_SG includes the outer-boundary-edges of the subdomain
end

allIdx_bound_e=find(att.bound_e==true);
%[size(allIdx_bound_e,1) size(allIdx_bound_e,2)]
allIdx_adjbounde_e=find(att.adj_bounde_e==true);
%[size(allIdx_adjbounde_e,1) size(allIdx_adjbounde_e,2)]
allIdx_stFI.e=[allIdx_bound_e;allIdx_adjbounde_e];
%[size(allIdx_stFI_e,1) size(allIdx_stFI_e,2)]
clearvars allIdx_bound_e allIdx_adjbounde_e

% return allIdx_stFI_f allIdx_stFI_e

subG_bin.f=subG_f_bin;
subG_bin.e=subG_e_bin;

subG_sizes.f=subG_f_sizes;
subG_sizes.e=subG_e_sizes;
denominator.SG=denominator_SG;

%disp('checkpoint golf')
disp('Divide_into_induced_subgraphs:ENDED')
end
