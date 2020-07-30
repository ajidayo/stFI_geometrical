function [SpElemProperties,STElemProperties,Num_of_Elem,PrimFacePos] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos)
%% UpdNum
UpdNum_SpV = SpElemProperties.SpV.UpdNum;
UpdNum_SpP = zeros(Num_of_Elem.SpP,1);
UpdNum_SpS = zeros(Num_of_Elem.SpS,1);
UpdNum_SpN = zeros(Num_of_Elem.SpN,1);
for SpPIdx=1:Num_of_Elem.SpP
    UpdNum_SpP(SpPIdx) = max(UpdNum_SpV(logical(sD(:,SpPIdx))));
end
for SpSIdx=1:Num_of_Elem.SpS
    UpdNum_SpS(SpSIdx) = max(UpdNum_SpP(logical(sC(:,SpSIdx))));
end
for SpNIdx=1:Num_of_Elem.SpN
    UpdNum_SpN(SpNIdx) = max(UpdNum_SpS(logical(sG(:,SpNIdx))));
end
SpElemProperties.SpP.UpdNum=UpdNum_SpP;
SpElemProperties.SpS.UpdNum=UpdNum_SpS;
SpElemProperties.SpN.UpdNum=UpdNum_SpN;

%% FirstSTNIdx, RefSpN
CurrentSTNIdx=1;
for SpNIdx=1:Num_of_Elem.SpN
    SpElemProperties.SpN.FirstSTNIdx(SpNIdx) = CurrentSTNIdx;
    for STNIdx = SpElemProperties.SpN.FirstSTNIdx(SpNIdx):...
            SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+UpdNum_SpN(SpNIdx)+1
        STElemProperties.STN.RefSpN(STNIdx) = SpNIdx;
    end
    CurrentSTNIdx = CurrentSTNIdx +UpdNum_SpN(SpNIdx) +2;
end
Num_of_Elem.STN = CurrentSTNIdx -1;
%% FirstSpSIdx
CurrentSTSIdx = 1;
for SpSIdx = 1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTSIdx(SpSIdx) = CurrentSTSIdx;
    CurrentSTSIdx = CurrentSTSIdx + UpdNum_SpS(SpSIdx) +1;
end
for SpNIdx=1:Num_of_Elem.SpN
    SpElemProperties.SpN.FirstSTSIdx(SpNIdx) = CurrentSTSIdx;
    CurrentSTSIdx = CurrentSTSIdx + UpdNum_SpN(SpNIdx) +1;    
end
Num_of_Elem.STS = CurrentSTSIdx -1;

%% FirstSTPIdx
CurrentSTPIdx=1;
for SpPIdx=1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTPIdx(SpPIdx) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpP(SpPIdx) +1;
end
for SpSIdx=1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTPIdx(SpSIdx) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpS(SpSIdx) +1;
end

Num_of_Elem.STP = CurrentSTPIdx -1;

%% FirstSTOmegaIdx
CurrentSTOmegaIdx =1;
for SpVIdx =1:Num_of_Elem.SpV
    SpElemProperties.SpV.FirstSTOmegaIdx(SpVIdx) = CurrentSTOmegaIdx;
    CurrentSTOmegaIdx = CurrentSTOmegaIdx + SpElemProperties.SpV.UpdNum(SpVIdx) +1;
end
for SpSIdx =1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTOmegaIdx(SpSIdx) = CurrentSTOmegaIdx;
    CurrentSTOmegaIdx = CurrentSTOmegaIdx + SpElemProperties.SpP.UpdNum(SpSIdx) +1;
end
Num_of_Elem.STOmega = CurrentSTOmegaIdx-1;

%% STV
CurrentSTVIdx=0;
for SpVIdx=1:Num_of_Elem.SpV
    UpdNum = UpdNum_SpV(SpVIdx);
    for TimeIdx = -0.5/UpdNum
        CurrentSTVIdx = CurrentSTVIdx+1;
        SpElemProperties.SpV.FirstSTVIdx(SpVIdx) = CurrentSTVIdx;
        STElemProperties.STV.TimeIdx(CurrentSTVIdx) = TimeIdx;
        STElemProperties.STV.RefSpV(CurrentSTVIdx)  = SpVIdx;
    end
    for TimeIdx = 0.5/UpdNum:1/UpdNum:(UpdNum-0.5)/UpdNum
        CurrentSTVIdx = CurrentSTVIdx+1;
        STElemProperties.STV.TimeIdx(CurrentSTVIdx) = TimeIdx;
        STElemProperties.STV.RefSpV(CurrentSTVIdx)  = SpVIdx;
    end
end
Num_of_Elem.STV = CurrentSTVIdx;

%% tasktype
% Usage Guide: SpElemPropreties.SpP.Belong_to_ST_FI(SpPIdx)
% First check if SpPs are on dt-Interfaces. Interface SpPs belongs to ST_FI-region.
% SpPs adjacent to interface-SpPs belongs to ST_FI-region.
SpElemProperties.SpP.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpP));
SpElemProperties.SpS.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpS));
SpElemProperties.SpS.UpdNumBoundary  = logical(sparse(1,Num_of_Elem.SpS));
AdjM_SpP_via_SpV = sD.'*sD;
for SpPIdx = 1:Num_of_Elem.SpP
    Sum = 0;
    IncSpV_List = find(sD(:,SpPIdx)).';
    for IncSpV = IncSpV_List
        Sum = Sum + SpElemProperties.SpV.UpdNum(IncSpV);
    end
    if Sum == size(IncSpV_List,2)*SpElemProperties.SpV.UpdNum(IncSpV_List(1))
        SpElemProperties.SpP.Belong_to_ST_FI(SpPIdx) = false;
    else
        SpElemProperties.SpP.Belong_to_ST_FI(SpPIdx) = true;
        SpElemProperties.SpP.Belong_to_ST_FI(logical(AdjM_SpP_via_SpV(SpPIdx,:))) = true;
        for IncSpS = find(sC(SpPIdx,:))
            SpElemProperties.SpS.UpdNumBoundary(IncSpS)  = true;
            SpElemProperties.SpS.Belong_to_ST_FI(IncSpS) = true;
            for AdjSpS = find(logical(sG(IncSpS,:))*logical(sG).')
                SpElemProperties.SpS.Belong_to_ST_FI(AdjSpS) = true;
            end
        end
    end
end

% SpSs incident to ST_FI-SpPs belongs to ST_FI-region.
% for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
%     for IncSpS = find(sC(SpPIdx,:))
%         SpElemProperties.SpS.Belong_to_ST_FI(IncSpS) = true;
%     end
% end
SpElemProperties.SpN.Belong_to_ST_FI = logical(sparse(1,Num_of_Elem.SpN));
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    for IncSpN = find(sG(SpSIdx,:))
        SpElemProperties.SpN.Belong_to_ST_FI(IncSpN) = true;
    end
end
% Other SpPs, and SpSs belongs to Sp_FI-region.
% ST_FI-SpSs which are incident to Sp_FI SpPs have 'SpFI_BoundarySpS' attribute (in this neccesary??) 

%dummy 
SpElemProperties.SpS.UpdNumCorner = logical(sparse(1,Num_of_Elem.SpS));

%% position of Primal Faces
PrimFacePos(Num_of_Elem.SpP).Vec = [0;0;0];
for SpPIdx = 1:Num_of_Elem.SpP
    Nodes = find(( logical(sG).'*logical(sC(SpPIdx,:)).' ).');
    PosVec_SpP = [0;0;0];
    for SpNIdx = Nodes
        PosVec_SpP = PosVec_SpP+NodePos.Prim(SpNIdx).Vec;
    end
    PrimFacePos(SpPIdx).Vec = PosVec_SpP/size(Nodes,2);
end
end