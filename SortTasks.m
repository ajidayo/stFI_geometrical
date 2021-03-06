function TaskOrder = SortTasks(TaskDepGraph,Map_SpElem_to_FirstGlobTask,sC,SpElemProperties)

%% add edges to task dependency graph: Dependencies between SpFI tasks and STFI tasks
EdgeIdx =0;
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    SpFI_RegionIdx = max(SpElemProperties.SpP.SpFI_TaskIdx(logical(sC(:,SpSIdx)).'));
    if SpFI_RegionIdx > 0
        for CurrentTimeSec = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
            EdgeIdx = EdgeIdx+1;
            StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)+CurrentTimeSec-1;
            TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_RegionIdx)-1+CurrentTimeSec;
        end
        for CurrentTimeSec = 2:SpElemProperties.SpS.UpdNum(SpSIdx)
            EdgeIdx = EdgeIdx+1;
            StaTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpFI_RegionIdx(SpFI_RegionIdx)-1+CurrentTimeSec-1;
            TgtTask(EdgeIdx) = Map_SpElem_to_FirstGlobTask.SpSIdx(SpSIdx)+CurrentTimeSec;
        end
    end
end
if exist('StaTask','var') == 1
    TaskDepGraph = addedge(TaskDepGraph,StaTask,TgtTask);
end
clearvars StaTask TgtTask
if numedges(TaskDepGraph) == 0
    TaskOrder = 1;
else
    TaskOrder = toposort(TaskDepGraph);
end
end