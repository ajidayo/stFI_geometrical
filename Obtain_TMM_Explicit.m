function [TMM_Explicit] ...
    = Obtain_TMM_Explicit(kappaoverZ,sC,UpdateNum,allIdx_stFI,subG_bin,subG_sizes,att,first_pIdx,MeshNum)

disp('Obtain_TMM_Explicit: CALLED')

% this unit generates the explicit time-marching matrix combining the
% space-time and spatial FI formulation.
% the time-marching within the subgraph is expressed in the FI formulation 
% and treated as one whole task.
% Therefore the task-dependence graph is altered from the full-space-time
% FI version.

task=struct('typ',[],'p_tgt',[],'omega_tgt',[],'SG_tgt',[],'Tsec_tgt',[]);

sum=0;
offset_FItask=zeros(size(subG_sizes.f,1),1);
for SGIdx=1:size(subG_sizes.f,2)
    offset_FItask(SGIdx)=sum;
    sum=sum+UpdateNum.SG(SGIdx);
end
FItaskNum=sum;

SG_for_FItask=zeros(FItaskNum,1);
Timesection_for_FItask=zeros(FItaskNum,1);
FItask=0;
edgecounter=0;
for SGIdx=1:size(subG_sizes.f,2)
    for i=1:UpdateNum.SG(SGIdx)
        FItask=FItask+1;
        task(FItask).typ="FI";
        SG_for_FItask(FItask)=SGIdx;
        Timesection_for_FItask(FItask)=i;
        if i>=2
            edgecounter=edgecounter+1;
            s(edgecounter)=FItask-1;
            t(edgecounter)=FItask;
        end
    end
end

taskIdx_for_p = sparse(MeshNum.P,1);% gives the taskIdx corresponding to the calculation of p
p_for_taskIdx = sparse(FItaskNum,1);% gives the p corresponding to taskIdx
% (iff taskIdx=1,...,FItaskNum, equals zero)
omega_for_taskIdx = sparse(FItaskNum,1);% gives the omega corresponding to each taskIdx
% (iff taskIdx=1,...,FItaskNum, equals zero)
taskIdx_newest=FItaskNum;

allIdx_bound_SG_e=find(att.e.bound_of_SG==true);
for ee=1:size(allIdx_bound_SG_e,1) %SGboundary edges
    e=allIdx_bound_SG_e(ee);
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1) % faces incident to edge
        f=row_sC_e(ff);
        if att.f.inc_bounde(f) == true %not SG face 
            SG_tgt=subG_bin.e(e);
            %disp(SG_tgt)
            for i=1:UpdateNum.SG(SG_tgt)
                edgecounter=edgecounter+1;
                if i==1
                    p_is_init=true;
                else
                    p_is_init=false;
                end
                FItask_tgt=offset_FItask(SG_tgt)+i;
                p_s=first_pIdx.f(f)+i-1;
                [s(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
                    = semi_add_dep_edge(p_s,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,p_is_init);               
                t(edgecounter)=FItask_tgt;
            end
        end
    end
end

%disp('checkpoint alpha')

%% generate D for the boundary edges/faces and generate task-dependence graph

disp('Constructing Taskdependence graph, D and Ctrans')
% N.B. taskIdx =1,...,FItaskNum corresponds to FItask
% and taskIdx =FItaskNum+1,... corresponds to p-calculation task
%%%%%% DO use structures to represent tasks %%%%%
% like; task.type = {stFI, FI, PML, PEC, ...}
% task.Idx = hoge

D=sparse(1,MeshNum.P);
Ctrans=sparse(1,MeshNum.P);

% omega is the row-index of D_tildeD_Zinverse; 
% each omega corresponds to an calculation of stFI variable
omega=0;
% Each edgecounter corresponds to a dependence relation between tasks

for ff=1:size(allIdx_stFI.f,1)
   f=allIdx_stFI.f(ff);
   sC_f=sC(f,:);
   col_sC_f=find(sC_f);
   for i=1:UpdateNum.f(f)
       omega=omega+1;
       edgecounter=edgecounter+1;
       if i==1
           p_s_is_init=true;
       else
           p_s_is_init=false;
       end
       p_s=first_pIdx.f(f)+i-1;
       [s(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
           = semi_add_dep_edge(p_s,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,p_s_is_init);
       p_t=first_pIdx.f(f)+i;
       [t(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
           = semi_add_dep_edge(p_t,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,false      );
       D(omega,p_s)=-1;
       D(omega,p_t)= 1; %the variable to be calculated at line Omega
       omega_for_taskIdx(taskIdx_for_p(p_t))=omega;
       for ee=1:size(col_sC_f,2)
           e=col_sC_f(ee);
           r_e_over_f=UpdateNum.e(e)/UpdateNum.f(f);
           for j=r_e_over_f*(i-1)+1:r_e_over_f*i
               edgecounter=edgecounter+1;
               p_s=first_pIdx.e(e)+j;
               if att.e.bound_of_SG(e)==true % If e is an outer boundary of a FI region,
                   s(edgecounter)=offset_FItask(subG_bin.e(e))+j;
                   %[f,e,p_s] %disp
               elseif att.e.adj_bounde(e)==true || att.e.boundaryedge(e)==true
                   [s(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
                       = semi_add_dep_edge(p_s,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,false);
                   %[f,e,p_s] %disp
               else
                   disp('error')
                   %[f,e,p_s] %disp
                   pause
               end
               p_t=first_pIdx.f(f)+i;
               [t(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
                   = semi_add_dep_edge(p_t,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,false);
               D(omega, p_s)=sC(f,e);
           end % j
       end % ee
   end % i
end %ff
OmegaNum=omega;

%disp('checkpoint bravo')

omega=0;
for ee=1:size(allIdx_stFI.e,1)
   e=allIdx_stFI.e(ee);
   sC_e=sC(:,e);
   row_sC_e=find(sC_e);
   for i=0:UpdateNum.e(e)-1
       if i==0
           p_s_is_init=true;
       else
           p_s_is_init=false;
       end
       omega=omega+1;
       edgecounter=edgecounter+1;
       p_s=first_pIdx.e(e)+i;
       [s(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
           = semi_add_dep_edge(p_s,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,p_s_is_init);
       p_t=first_pIdx.e(e)+i+1;
       [t(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
           = semi_add_dep_edge(p_t,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,false);
       Ctrans(omega,p_s)=-1;
       Ctrans(omega,p_t)=1; %the variable to be calculated at line Omega
       omega_for_taskIdx(taskIdx_for_p(p_t))=omega+OmegaNum;
       for ff=1:size(row_sC_e,1)
           f=row_sC_e(ff);
           % for debugging
%           if att.f.inc_bounde(f)~=true
%               disp('error')
%               pause
%           end
           for j=0:UpdateNum.f(f)-1
               if UpdateNum.e(e)*j==UpdateNum.f(f)*i
                   edgecounter=edgecounter+1;
                   p_s=first_pIdx.f(f)+j;
                   [s(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
                       = semi_add_dep_edge(p_s,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,p_s_is_init);
                   p_t=first_pIdx.e(e)+i+1;
                   [t(edgecounter),taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task]...
                       = semi_add_dep_edge(p_t,taskIdx_for_p,p_for_taskIdx,taskIdx_newest,task,false);
                   Ctrans(omega,p_s)=sC(f,e);
                   break;
               end %if
           end %j
       end %ff
   end %i
end %ee

OmegaDualNum = omega;

taskNum = taskIdx_newest;

% omega is the row-index of D_tildeD_Zinverse; 
% each omega corresponds to an calculation of stFI variable
D_tildeD_Zinv=[D;Ctrans * spdiags(kappaoverZ,0,MeshNum.P,MeshNum.P)];

%clearvars taskIdx_for_p
clearvars D Ctrans

%disp('checkpoint charlie')


%% obtain the topological order of the tasks
Taskdependence = digraph;
Taskdependence = addedge(Taskdependence,s,t);

%% generate TMM_Explicit; Notes
% the task-dependence graph, D and TMM_FI(i) 
% (i=1,2,...,Number of subgraphs)
% is transformed into TMM_Explicit
% the transformation is alike to that in "planB"
% cf. Obtain_TimemarchingMatrix_Explicit_planB

disp('Calculating orders')

Taskorder=toposort(Taskdependence);
%Taskorder=toposort(Taskdependence,'Order','stable');

clearvars Taskdependence s t

disp('Constructing Time Marching Matrix')

TMM_Intermidiate=spalloc(MeshNum.P,MeshNum.P,5*MeshNum.P);
TMM_Intermidiate=TMM_Intermidiate+speye(MeshNum.P);

%disp("checkpoint charlie-2")

for taskIdx=1:taskNum
    %disp(taskIdx)
    task_tgt=Taskorder(taskIdx);
    if task(task_tgt).typ == "InitVal"
        continue
    elseif task(task_tgt).typ == "FI"
%        disp('FI')
        SG_tgt = SG_for_FItask(task_tgt);
        Timesection_tgt = Timesection_for_FItask(task_tgt);
        % call FI for SG_tgt, Timesection_tgt
        TMM_Intermidiate= ...
            TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappaoverZ,subG_bin.f,subG_bin.e,first_pIdx)...
            *TMM_Intermidiate;
       % T=TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappatimesz,subG_bin.f,subG_bin.e,first_p);
    elseif task(task_tgt).typ=="stFI"
        %disp('stFI')
        p_tgt = p_for_taskIdx(task_tgt);
        omega_tgt = omega_for_taskIdx(task_tgt);
        % call stFI for p_tgt, Omega_to_calculate
        TMM_Intermidiate= ...
            TMM_stFI_single_p(p_tgt, omega_tgt, D_tildeD_Zinv) ...
            * TMM_Intermidiate;
    else
        disp('error: task-type undefined')
        pause
    end
end

clearvars Taskorder


%disp('checkpoint delta')


%% generate TMM_Explicit
% Note that the (intermidiate?) variables are already bypassed in
% TMM_Intermidiate.
% For p which is an (intermidiate?) variable, delete the p-th row.

Store=logical(sparse([],[],[],MeshNum.F+MeshNum.E,MeshNum.P,MeshNum.F+MeshNum.E));
for f=1:MeshNum.F
    p_tgt=first_pIdx.f(f)+UpdateNum.f(f);
    Store(f,p_tgt)=true;
end 
for e=1:MeshNum.E
    p_tgt=first_pIdx.e(e)+UpdateNum.e(e);
    Store(MeshNum.F+e,p_tgt)=true;
end
Store_Init=logical(sparse([],[],[],MeshNum.F+MeshNum.E,MeshNum.P,MeshNum.F+MeshNum.E));
for f=1:MeshNum.F
    p_tgt=first_pIdx.f(f);
    Store_Init(f,p_tgt)=true;
end 
for e=1:MeshNum.E
    p_tgt=first_pIdx.e(e);
    Store_Init(MeshNum.F+e,p_tgt)=true;
end

TMM_Explicit = Store*TMM_Intermidiate*Store_Init.';

clearvars TMM_Intermidiate

% return TMM_Explicit Store (clearvars -except TMM_Explicit Store)

disp('TMM_Explicit calculated')

%disp('checkpoint echo')

disp('Obtain_TMM_Explicit: ENDED')
