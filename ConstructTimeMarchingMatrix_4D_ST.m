function TMM = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem)
Sum=0;
for nth_source = 1:size(Source,2) 
    Sum=Sum+Source(nth_source).UpdNum;
end
TMM_Redundant = speye(Num_of_Elem.STP+Sum);
for nth_task = 1:size(TaskOrder,2)
    taskIdx = TaskOrder(nth_task); 
    switch Task(taskIdx).Type
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

TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties,Num_of_Elem,Sum);

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
    SingleTask_SpFI(SpFI_TaskInfo,sC,Kappa_over_Z,Source,SpElemProperties,Num_of_Elem,TMM_Redundant)
disp('SingleTask_SpFI CALLED')
%TMM_Onestep      = speye(size(TMM_Redundant,1),size(TMM_Redundant,2));

Sp_FI_Region_tgt = SpFI_TaskInfo.Sp_FI_Region_tgt;
TimeSection_tgt  = SpFI_TaskInfo.TimeSection_tgt;

STP_to_Conserve_Prim    = true(Num_of_Elem.STP,1);
STP_to_Conserve_Dual    = true(Num_of_Elem.STP,1);
Map_STPFutr_SpP         = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpP));
Map_SpP_STPPast         = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));
Map_STPFutr_SpS         = logical(sparse(Num_of_Elem.STP,Num_of_Elem.SpS));
Map_SpS_STPPast         = logical(sparse(Num_of_Elem.SpS,Num_of_Elem.STP));
Map_SpPinctoSpS_STPPast = logical(sparse(Num_of_Elem.SpP,Num_of_Elem.STP));
Map_SpSinctoSpP_STPFutr = logical(sparse(Num_of_Elem.SpS,Num_of_Elem.STP));

for SpPIdx  = SpFI_TaskInfo.ElemIdx.SpP
    STPFutr = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt;
    STPPast = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt-1;
    Map_STPFutr_SpP(STPFutr,SpPIdx)    = true;
    Map_SpP_STPPast(SpPIdx,STPPast)    = true;
    STP_to_Conserve_Prim(STPFutr)      = false;
    for SpSIdx  = find(sC(SpPIdx,:))
        STPFutr = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt;
        Map_SpSinctoSpP_STPFutr(SpSIdx,STPFutr)    = true;
    end
end
for SpSIdx  = SpFI_TaskInfo.ElemIdx.SpS
    STPFutr = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt;
    STPPast = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+TimeSection_tgt-1;
    Map_STPFutr_SpS(STPFutr,SpSIdx)    = true;
    Map_SpS_STPPast(SpSIdx,STPPast)    = true;
    STP_to_Conserve_Dual(STPFutr)      = false;
    for SpPIdx  = find(sC(:,SpSIdx)).'
        STPPast = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+TimeSection_tgt-1;
        Map_SpPinctoSpS_STPPast(SpPIdx,STPPast)    = true;
    end
end

z       = spdiags(Kappa_over_Z.^(-1),0,Num_of_Elem.STP,Num_of_Elem.STP);
zinv    = spdiags(Kappa_over_Z      ,0,Num_of_Elem.STP,Num_of_Elem.STP);

Sum=0;
for SourceIdx = 1:size(Source,2) 
    Sum=Sum+Source(SourceIdx).UpdNum;
end
SourcePattern = sparse(Num_of_Elem.SpS,Sum);
ST_SourceIdx = 0;
for SourceIdx = 1:size(Source,2)
    for CurrentTimeSec = 1:Source(SourceIdx).UpdNum 
        ST_SourceIdx = ST_SourceIdx+1;
        if SpElemProperties.SpS.SpFI_TaskIdx(Source(SourceIdx).DualFace_tgt)==Sp_FI_Region_tgt ...
                && CurrentTimeSec == TimeSection_tgt
            SourcePattern(Source(SourceIdx).DualFace_tgt,ST_SourceIdx)=1;
        end
    end
end
% the minus for the dual grid is caused by the sign [z] has.
TMM_Onestep     = ...
    [z*Map_STPFutr_SpS*(Map_SpS_STPPast - sC.'*Map_SpPinctoSpS_STPPast)*zinv, ...
    z*Map_STPFutr_SpS*SourcePattern;...
    sparse(Num_of_Elem.STP,size(SourcePattern,2)).',sparse(size(SourcePattern,2),size(SourcePattern,2))...
    ];
TMM_Onestep     = TMM_Onestep ...
    + ...
    [spdiags(STP_to_Conserve_Dual,0,Num_of_Elem.STP,Num_of_Elem.STP), ...
    sparse(Num_of_Elem.STP,size(SourcePattern,2));...
    sparse(Num_of_Elem.STP,size(SourcePattern,2)).',speye(size(SourcePattern,2))...
    ];
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

TMM_Onestep     = ...
    [  Map_STPFutr_SpP*(Map_SpP_STPPast - sC  *Map_SpSinctoSpP_STPFutr), ...
    sparse(Num_of_Elem.STP,size(SourcePattern,2));...
    sparse(Num_of_Elem.STP,size(SourcePattern,2)).',sparse(size(SourcePattern,2),size(SourcePattern,2))...
    ];
TMM_Onestep     = TMM_Onestep +...
    [spdiags(STP_to_Conserve_Prim,0,Num_of_Elem.STP,Num_of_Elem.STP), ...
    sparse(Num_of_Elem.STP,size(SourcePattern,2));...
    sparse(Num_of_Elem.STP,size(SourcePattern,2)).',speye(size(SourcePattern,2))];
TMM_Redundant   = TMM_Onestep * TMM_Redundant;

end

%%
function TMM = RemoveRedundantDoF(TMM_Redundant,SpElemProperties,Num_of_Elem,Sum)
Store=logical(sparse([],[],[],Num_of_Elem.SpP+Num_of_Elem.SpS,Num_of_Elem.STP+Sum,Num_of_Elem.SpP+Num_of_Elem.SpS));
for SpPIdx=1:Num_of_Elem.SpP
    STP_tgt=SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx);
    Store(SpPIdx,STP_tgt)=true;
end 
for SpSIdx=1:Num_of_Elem.SpS
    STP_tgt=SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx);
    Store(Num_of_Elem.SpP+SpSIdx,STP_tgt)=true;
end
Store_Init=logical(sparse([],[],[],Num_of_Elem.SpP+Num_of_Elem.SpS+Sum,Num_of_Elem.STP+Sum,Num_of_Elem.SpP+Num_of_Elem.SpS));
for SpPIdx=1:Num_of_Elem.SpP
    STP_tgt=SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
    Store_Init(SpPIdx,STP_tgt)=true;
end 
for SpSIdx=1:Num_of_Elem.SpS
    STP_tgt=SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
    Store_Init(Num_of_Elem.SpP+SpSIdx,STP_tgt)=true;
end
for ST_SourceIdx =1:Sum
    Store_Init(Num_of_Elem.SpP+Num_of_Elem.SpS+ST_SourceIdx,Num_of_Elem.STP+ST_SourceIdx)=true;
end
TMM = Store*TMM_Redundant*Store_Init.';
end