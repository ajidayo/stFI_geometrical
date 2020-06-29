function [Task,TaskDependanceGraph,Map_SpElem_to_FirstGlobTask] ...
    = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemPropreties,Task,TaskDependanceGraph,Map_SpElem_to_FirstGlobTask)
sDPattern_PartiallyOmitted=logical(sD);

for SpP = find(SpElemPropreties.SpP.Belong_to_ST_FI)
    sDPattern_PartiallyOmitted(:,SpP)=false;
end
Adj_SpP_PartiallyOmitted = graph(sDPattern_PartiallyOmitted.' * sDPattern_PartiallyOmitted,'omitselfloops');


[SG_bin_SpP,SG_ElemNum_SpP] = conncomp(Adj_SpP_PartiallyOmitted);
clearvars Adj_SpP_PartiallyOmitted

[SpFI_TaskInfo] = EliminateST_FI_SpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties);
clearvars SG_bin_SpP SG_ElemNum_SpP
SpFI_TaskInfo.ElemIdx.SpS = zeros(Num_of_Elem.SpS,1);
SpFI_TaskInfo.ElemNum.SpS = zeros(size(SpFI_TaskInfo.ElemNum.SpP,2),1);
for SpS_tgt=1:Num_of_Elem.SpS
    switch SpElemProperties.SpS.Belong_to_ST_FI(SpS_tgt) || SpElemProperties.SpS.PEC(SpS_tgt)
        case true
            continue;
        case false
            SpP_fetch = find(sC(:,SpS_tgt),1);
            SpFI_TaskIdx = SpElemProperties.SpP.SpFI_TaskIdx(SpP_fetch);
            Num_of_IncludedSpS = SpFI_TaskInfo(SpFI_TaskIdx).ElemNum.SpS+1;
            SpFI_TaskInfo(SpFI_TaskIdx).ElemNum.SpS = Num_of_IncludedSpS;
            SpFI_TaskInfo(SpFI_TaskIdx).ElemIdx.SpS(Num_of_IncludedSpS) = SpS_tgt;
            SpElemProperties.SpS.SpFI_TaskIdx(SpS_tgt)=SpFI_TaskIdx;
    end
end
switch size(Task,2)
    case 1
        GlobalTaskIdx = 0;
    otherwise
        GlobalTaskIdx = size(Task,2);
end
for SpFI_TaskIdx=1:size(SpFI_Task,2)
     SpFI_TaskInfo(SpFI_TaskIdx).UpdNum = ...
         SpElemProperties.SpP.UpdNum(SpFI_TaskInfo(SpFI_TaskIdx).ElemIdx.SpP(1));
     for CurrentTimeSec = 1:SpFI_TaskInfo(SpFI_TaskIdx).UpdNum
         GlobalTaskIdx = GlobalTaskIdx+1;
         if CurrentTimeSec == 1
             Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx) = GlobalTaskIdx;
         end
         Task(GlobalTaskIdx) = SpFI_TaskInfo(SpFI_TaskIdx);
         Task(GlobalTaskIdx).Type = "Sp_FI";
         Task(GlobalTaskIdx).TimeSection_tgt = CurrentTimeSec;
     end
end

EdgeIdx =0;
for SpFI_TaskIdx = size(SpFI_TaskInfo,2)
    for CurrentTimeSec = 2:SpFI_TaskIdx.UpdNum(SpFI_TaskIdx)
        EdgeIdx = EdgeIdx+1;
        StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx)-1+CurrentTimeSec-1;
        TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_TaskIdx)-1+CurrentTimeSec;
    end
end
TaskDependanceGraph = addedge(TaskDependanceGraph,StaTask,TgtTask);
clearvars StaTask TgtTask

end

%% 
function [SpFI_TaskInfo,SpElemProperties] = EliminateST_FI_SpPs(SG_bin_SpP,SG_ElemNum_SpP,SpElemProperties)
SpElemProperties.SpP.SpFI_TaskIdx=zeros(Num_of_Elem.SpP,1);
SpFI_TaskIdx=0;
for SGIdx=1:size(SG_ElemNum_SpP,2)
    LogiIdx=find(SG_bin_SpP==SGIdx);
    SpP_test=LogiIdx(1);
    switch SpElemProperties.SpP.Belong_to_ST_FI(SpP_test)
        case true
            SpElemProperties.SpP.SpFI_TaskIdx(SpP_test)=0;
        case false
            SpFI_TaskIdx=SpFI_TaskIdx+1;
            SpFI_TaskInfo(SpFI_TaskIdx).Sp_FI_Region_tgt = SpFI_TaskIdx;
            SpFI_TaskInfo(SpFI_TaskIdx).ElemIdx.SpP=LogiIdx;
            SpFI_TaskInfo(SpFI_TaskIdx).ElemNum.SpP=SG_ElemNum_SpP(SGIdx);
            SpElemProperties.SpP.SpFI_TaskIdx(LogiIdx)=SpFI_TaskIdx;
    end
end
end