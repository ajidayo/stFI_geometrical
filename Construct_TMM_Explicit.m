function [TMM_Explicit] ...
    = Construct_TMM_Explicit(Taskorder,task,D_tildeD_Zinv,kappaoverZ,sC,UpdateNum,subG_bin,first_pIdx,MeshNum)
TMM_Intermidiate=spalloc(MeshNum.P,MeshNum.P,5*MeshNum.P);
TMM_Intermidiate=TMM_Intermidiate+speye(MeshNum.P);

%disp("checkpoint charlie-2")
taskNum=size(task,2);
for taskIdx=1:taskNum
    %disp(taskIdx)
    task_tgt=Taskorder(taskIdx);
    if task(task_tgt).typ == "InitVal"
        continue
    elseif task(task_tgt).typ == "FI"
%        disp('FI')
        SG_tgt = task(task_tgt).SG_tgt;
        Timesection_tgt = task(task_tgt).Tsec_tgt;
        % call FI for SG_tgt, Timesection_tgt
        TMM_Intermidiate= ...
            TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappaoverZ,subG_bin,first_pIdx,MeshNum)...
            *TMM_Intermidiate;
       % T=TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappatimesz,subG_bin.f,subG_bin.e,first_p);
    elseif task(task_tgt).typ=="stFI"
        %disp('stFI')
        p_tgt = task(task_tgt).p_tgt;
        omega_tgt = task(task_tgt).omega_tgt;
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

%disp('TMM_Explicit calculated')

%disp('checkpoint echo')

%disp('Obtain_TMM_Explicit: ENDED')
end