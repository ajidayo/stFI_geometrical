function [SpElemPropreties,Num_of_STP] = Properties_of_Sp_Elements(sG,sC,sD,SpElemPropreties,Num_of_Elem)
%% UpdNum
UpdNum_SpV = SpElemPropreties.SpV.UpdNum;
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
SpElemPropreties.SpP.UpdNum=UpdNum_SpP;
SpElemPropreties.SpS.UpdNum=UpdNum_SpS;
SpElemPropreties.SpN.UpdNum=UpdNum_SpN;
%% FirstSTPIdx
CurrentSTPIdx=1;
for nth_SpP=1:Num_of_Elem.SpP
    SpElemPropreties.SpP.FirstSTPIdx(nth_SpP) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpP +1 ;
end
for nth_SpS=1:Num_of_Elem.SpS
    SpElemPropreties.SpS.FirstSTPIdx(nth_SpS) = CurrentSTPIdx;
    CurrentSTPIdx = CurrentSTPIdx +UpdNum_SpS +1 ;
end
Num_of_STP = CurrentSTPIdx -1;
%% tasktype
% Usage Guide: SpElemPropreties.SpP.Belong_to_ST_FI(SpPIdx)
% first check if SpPs are on dt-Interfaces. Interface SpPs belongs to ST_FI-region.
% SpPs adjacent to interface-SpPs belongs to ST_FI-region.
% SpSs incident to ST_FI-SpPs belongs to ST_FI-region.
% Other SpPs, and SpSs belongs to Sp_FI-region.
% ST_FI-SpSs which are incident to Sp_FI SpPs have 'SpFI_BoundarySpS' attribute. 
%% other attributes
% SpSs on the outer boundary of the calculation region has 'PEC' attribute.
end