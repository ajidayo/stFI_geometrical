function TMM = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem)
TMM_Redundant = speye(Num_of_Elem.STP);
for nth_task = 1:size(TaskOrder,2)
    taskIdx = TaskOrder(nthtask); 
    switch Task(taskIdx).type
        case "ST_FI_Prim"
            TMM_Redundant = SingleTask_STFI_Prim(Task(taskIdx),D2,TMM_Redundant);
        case "ST_FI_Dual"
            TMM_Redundant = SingleTask_STFI_Dual(Task(taskIdx),D1,Kappa_over_Z,TMM_Redundant);
        case "Sp_FI"
            TMM_Redundant = SingleTask_SpFI(Task(taskIdx),sC,Kappa_over_Z,Source,SpElemProperties,Num_of_Elem,TMM_Redundant);
        otherwise
            disp("No such Tasktype defined.")
    end
end

TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties);

end

%%
function TMM_Redundant = SingleTask_STFI_Prim(STFI_PrimTask,D2,TMM_Redundant)

TMM_Onestep=speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(STFI_PrimTask.p_tgt,:) = ...
    [-(1.0/D2(STFI_PrimTask.omega_tgt,STFI_PrimTask.p_tgt))*D2(STFI_PrimTask.omega_tgt,:),zeros()];
TMM_Onestep(STFI_PrimTask.p_tgt,STFI_PrimTask.p_tgt)=0;
TMM_Redundant=TMM_Onestep*TMM_Redundant;
end

%%
function TMM_Redundant= SingleTask_STFI_Dual(STFI_DualTask,D1,Kappa_over_Z,TMM_Redundant)
TMM_Onestep=speye(size(TMM_Redundant,1),size(TMM_Redundant,2));
TMM_Onestep(STFI_DualTask.p_tgt,:) = ...
    [-(1.0/(D1(STFI_DualTask.omega_tgt,STFI_DualTask.p_tgt)*Kappa_over_Z(STFI_DualTask.p_tgt)))...
    *(D1(STFI_DualTask.omega_tgt,:).*Kappa_over_Z.'),...
    zeros()];
TMM_Onestep(STFI_DualTask.p_tgt,STFI_DualTask.p_tgt)=0;
TMM_Redundant=TMM_Onestep*TMM_Redundant;
end

%%
function TMM_Redundant = ...
    SingleTask_SpFI(SpFI_Task,sC,Kappa_over_Z,Source,SpElemProperties,Num_of_Elem,TMM_Redundant)
%TMM_Onestep      = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));

Sp_FI_Region_tgt = SpFI_Task.Sp_FI_Region_tgt;
TimeSection_tgt  = SpFI_Task.TimeSection_tgt;

STP_to_Conserve  = true(Num_of_Elem.STP,1);
Map_STPnext_SpP = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpP));
Map_SpP_STPpres = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));
Map_STPnext_SpS = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpP));
Map_SpS_STPpres = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));

Map_SpPinctoSpS_STPpres=logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));

% f_lst_SG is the index of "equal-to-SG_tgt" elements 
SpP_lst_SpFI = SpFI_TaskInfo(SpFI_TaskIdx).ElemNum.SpP;
SpS_lst_SpFI = SpFI_TaskInfo(SpFI_TaskIdx).ElemNum.SpS;
for nthSpP  = 1:size(SpP_lst_SpFI,1)
    SpP     = SpP_lst_SpFI(nthSpP);
    STPnext = SpElemProperties.SpP.FirstSTPIdx(SpP)+TimeSection_tgt;
    STPpres = SpElemProperties.SpP.FirstSTPIdx(SpP)+TimeSection_tgt-1;
    Map_STPnext_SpP(STPnext,SpP)    = true;
    Map_SpP_STPpres(SpP,STPpres)    = true;
    STP_to_Conserve(STPnext)        = false;
end
for nthSpS  = 1:size(SpS_lst_SpFI,1)
    SpS     = SpS_lst_SpFI(nthSpS);
    STPnext = SpElemProperties.SpS.FirstSTPIdx(SpS)+TimeSection_tgt;
    STPpres = SpElemProperties.SpS.FirstSTPIdx(SpS)+TimeSection_tgt-1;
    Map_STPnext_SpS(STPnext,SpS)    = true;
    Map_SpS_STPpres(SpS,STPpres)    = true;
    STP_to_Conserve(STPnext)        = false;
    row_sC_e    = find(sC(:,SpS));
    for nthSpP  = 1:size(row_sC_e,1)%faces incident to e
        SpP     = row_sC_e(nthSpP);
        STPpres = SpElemProperties.SpP.FirstSTPIdx(SpP)+TimeSection_tgt-1;
        Map_SpPinctoSpS_STPpres(SpP,STPpres)    = true;
    end
end

z       = spdiags(Kappa_over_Z.^(-1),0,Num_of_Elem.STP,Num_of_Elem.STP);
zinv    = spdiags(Kappa_over_Z      ,0,Num_of_Elem.STP,Num_of_Elem.STP);

ST_SourceIdx=0;
for nth_source = 1:size(Source,2) 
    ST_SourceIdx=ST_SourceIdx+Source(nth_source).UpdNum;
end
SourcePattern = sparse(Num_of_Elem.SpS,ST_SourceIdx);
ST_SourceIdx=0;
for nth_source = 1:size(Source,2)
    for CurrentTimeSec = 1:Source(nth_source).UpdNum 
        ST_SourceIdx=ST_SourceIdx+1;
        if Source(nth_source).SpFI_RegionIdx==Sp_FI_Region_tgt && CurrentTimeSec == TimeSection_tgt
            SourcePattern(Source(nth_source).DualFace_tgt,ST_SourceIdx)=1;
        end
    end
end

% the minus for the dual grid is because [z] includes the sign-flipping
%% Are these signes definitely right? check.
TMM_Onestep     = ...
    [z*Map_STPnext_SpS*(Map_SpS_STPpres -sC.'*Map_SpPinctoSpS_STPpres)*zinv, ...
    z*Map_STPnext_SpS*SourcePattern];
TMM_Onestep     = TMM_Onestep + ...
    [spdiags(STP_to_Conserve,0,Num_of_Elem.STP,Num_of_Elem.STP), ...
    zeros(Num_of_Elem.STP,size(SourcePattern,2))];
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

TMM_Onestep     = [Map_STPnext_SpP*(Map_SpP_STPpres - sC*(Map_STPnext_SpS.') ), ...
    zeros(Num_of_Elem.STP,size(SourcePattern,2))];
TMM_Onestep     = TMM_Onestep +[spdiags(STP_to_Conserve,0,Num_of_Elem.STP,Num_of_Elem.STP), ...
    zeros(Num_of_Elem.STP,size(SourcePattern,2))];
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

end

%%
function TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties)
Store=logical(sparse([],[],[],Num_of_Elem.SpP+Num_of_Elem.SpS,Num_of_Elem.STS,Num_of_Elem.SpP+Num_of_Elem.SpS));
for nth_SpP=1:Num_of_Elem.SpP
    p_tgt=SpElemProperties.SpP.FirstSTPIdx(SpP)+SpElemProperties.SpP.UpdNum(nth_SpP);
    Store(nth_SpP,p_tgt)=true;
end 
for nth_SpS=1:Num_of_Elem.SpS
    p_tgt=SpElemProperties.SpS.FirstSTPIdx(SpS)+SpElemProperties.SpS.UpdNum(nth_SpS);
    Store(Num_of_Elem.SpP+nth_SpS,p_tgt)=true;
end
Store_Init=logical(sparse([],[],[],Num_of_Elem.SpP+Num_of_Elem.SpS,Num_of_Elem.STS,Num_of_Elem.SpP+Num_of_Elem.SpS));
for nth_SpP=1:Num_of_Elem.SpP
    p_tgt=SpElemProperties.SpP.FirstSTPIdx(SpP);
    Store_Init(nth_SpP,p_tgt)=true;
end 
for nth_SpS=1:Num_of_Elem.SpS
    p_tgt=SpElemProperties.SpS.FirstSTPIdx(SpS);
    Store_Init(Num_of_Elem.SpP+nth_SpS,p_tgt)=true;
end

TMM = Store*TMM_Redundant*Store_Init.';
end