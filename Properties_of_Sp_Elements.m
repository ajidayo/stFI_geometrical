function [SpElemProperties,Num_of_STP] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem)
%% UpdNum
UpdNum_SpV = SpElemProperties.SpV.UpdNum;
UpdNum_SpP = zeros(Num_of_Elem.SpP,1);
UpdNum_SpS = zeros(Num_of_Elem.SpS,1);
UpdNum_SpN = zeros(Num_of_Elem.SpN,1);
for SpP=1:Num_of_Elem.SpP
    UpdNum_SpP(SpP) = max(UpdNum_SpV(sD(:,SpP)));
end
for SpS=1:Num_of_Elem.SpS
    UpdNum_SpS(SpS) = max(UpdNum_SpP(sC(:,SpS)));
end
for SpN=1:Num_of_Elem.SpN
    UpdNum_SpN(SpN) = max(UpdNum_SpS(sG(:,SpN)));
end
SpElemProperties.SpP.UpdNum=UpdNum_SpP;
SpElemProperties.SpS.UpdNum=UpdNum_SpS;
SpElemProperties.SpN.UpdNum=UpdNum_SpN;
%% FirstSTPIdx
CurrentSTPIdx=1;
for nthSpP=1:Num_of_Elem.SpP
    SpElemProperties.SpP.FirstSTPIdx(nthSpP) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpP +1 ;
end
for nth_SpS=1:Num_of_Elem.SpS
    SpElemProperties.SpS.FirstSTPIdx(nth_SpS) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpS +1 ;
end
Num_of_STP = CurrentSTPIdx -1;
%% tasktype
% Usage Guide: SpElemPropreties.SpP.Belong_to_ST_FI(SpPIdx)
% First check if SpPs are on dt-Interfaces. Interface SpPs belongs to ST_FI-region.
% SpPs adjacent to interface-SpPs belongs to ST_FI-region.
AdjM_SpP_via_SpV = sD.'*sD;
for nthSpP = 1:Num_of_Elem.SpP
    Sum = 0;
    IncSpV_List = find(sD(:,nthSpP));
    for IncSpV = IncSpV_List
        Sum = Sum + SpElemProperties.SpV.UpdNum(IncSpV);  
    end
    switch Sum == size(IncSpV_List,1)*SpElemProperties.SpV.UpdNum(IncSpV_List(1))
        case true
            SpElemProperties.SpP.Belong_to_ST_FI(nthSpP) = false;
        case false 
            SpElemProperties.SpP.Belong_to_ST_FI(nthSpP) = true;
            SpElemProperties.SpP.Belong_to_ST_FI(find(AdjM_SpP_via_SpV(SpP,:))) = true;
    end
end
% SpSs incident to ST_FI-SpPs belongs to ST_FI-region.
for nth_SpP = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpS = find(sC(nthSpP,:))
        SpElemProperties.SpS.Belong_to_ST_FI(IncSpS) = true;
    end
end
% Other SpPs, and SpSs belongs to Sp_FI-region.
% ST_FI-SpSs which are incident to Sp_FI SpPs have 'SpFI_BoundarySpS' attribute (in this neccesary??) 

%% other attributes
% SpSs on the outer boundary of the calculation region has 'PEC' attribute.
end