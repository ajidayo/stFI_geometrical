function [SpElemProperties,Num_of_STP] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos)
%% UpdNum
UpdNum_SpV = SpElemProperties.SpV.UpdNum;
UpdNum_SpP = zeros(Num_of_Elem.SpP,1);
UpdNum_SpS = zeros(Num_of_Elem.SpS,1);
UpdNum_SpN = zeros(Num_of_Elem.SpN,1);
for SpP=1:Num_of_Elem.SpP
    UpdNum_SpP(SpP) = max(UpdNum_SpV(find(sD(:,SpP))));
end
for SpS=1:Num_of_Elem.SpS
    UpdNum_SpS(SpS) = max(UpdNum_SpP(find(sC(:,SpS))));
end
for SpN=1:Num_of_Elem.SpN
    UpdNum_SpN(SpN) = max(UpdNum_SpS(find(sG(:,SpN))));
end
SpElemProperties.SpP.UpdNum=UpdNum_SpP;
SpElemProperties.SpS.UpdNum=UpdNum_SpS;
SpElemProperties.SpN.UpdNum=UpdNum_SpN;
%% FirstSTPIdx
CurrentSTPIdx=1;
for nth_SpP=1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTPIdx(nth_SpP) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpP(nth_SpP) +1 ;
end
for nth_SpS=1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTPIdx(nth_SpS) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpS(nth_SpS) +1 ;
end
Num_of_STP = CurrentSTPIdx -1;
%% tasktype
% Usage Guide: SpElemPropreties.SpP.Belong_to_ST_FI(SpPIdx)
% First check if SpPs are on dt-Interfaces. Interface SpPs belongs to ST_FI-region.
% SpPs adjacent to interface-SpPs belongs to ST_FI-region.
SpElemProperties.SpP.Belong_to_ST_FI = logical(sparse(Num_of_Elem.SpP,1));
AdjM_SpP_via_SpV = sD.'*sD;
for nth_SpP = 1:Num_of_Elem.SpP
    Sum = 0;
    IncSpV_List = find(sD(:,nth_SpP));
    for IncSpV = IncSpV_List.'
        Sum = Sum + SpElemProperties.SpV.UpdNum(IncSpV);
    end
    if Sum == size(IncSpV_List,1)*SpElemProperties.SpV.UpdNum(IncSpV_List(1))
        SpElemProperties.SpP.Belong_to_ST_FI(nth_SpP) = false;
    else
        SpElemProperties.SpP.Belong_to_ST_FI(nth_SpP) = true;
        SpElemProperties.SpP.Belong_to_ST_FI(find(AdjM_SpP_via_SpV(SpP,:))) = true;
    end
end
SpElemProperties.SpS.Belong_to_ST_FI = logical(sparse(Num_of_Elem.SpS,1));
% SpSs incident to ST_FI-SpPs belongs to ST_FI-region.
for nth_SpP = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpS = find(sC(nth_SpP,:))
        SpElemProperties.SpS.Belong_to_ST_FI(IncSpS) = true;
    end
end
% Other SpPs, and SpSs belongs to Sp_FI-region.
% ST_FI-SpSs which are incident to Sp_FI SpPs have 'SpFI_BoundarySpS' attribute (in this neccesary??) 

%% position of Primal Faces
for SpPIdx = 1:Num_of_Elem.SpP
    Nodes = find(logical(sC(SpPIdx,:))*logical(sG));
    PosVec_SpP = [0;0;0];
    for SpNIdx = Nodes
        PosVec_SpP = PosVec_SpP+NodePos.Prim(SpNIdx).Vec;
    end
    SpElemProperties.SpP.Position.x(SpPIdx) = PosVec_SpP(1)/size(Nodes,2); 
    SpElemProperties.SpP.Position.y(SpPIdx) = PosVec_SpP(2)/size(Nodes,2); 
    SpElemProperties.SpP.Position.z(SpPIdx) = PosVec_SpP(3)/size(Nodes,2); 
end
end