function [D0,D1,D2,D3] = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem)
Num_of_Elem.STV     = sum(SpElemProperties.SpV.UpdNum);
Num_of_Elem.STOmega = sum(SpElemProperties.SpV.UpdNum) + 2*Num_of_Elem.SpV +sum(SpElemProperties.SpP.UpdNum);
Num_of_Elem.STS     = sum(SpElemProperties.SpS.UpdNum) + 2*Num_of_Elem.SpS + sum(SpElemProperties.SpN.UpdNum) + Num_of_Elem.SpN;
Num_of_Elem.STN     = sum(SpElemProperties.SpN.UpdNum) + 2*Num_of_Elem.SpN;
D0 = sparse(Num_of_Elem.STS     ,Num_of_Elem.STN    );
D1 = sparse(Num_of_Elem.STP     ,Num_of_Elem.STS    );
D2 = sparse(Num_of_Elem.STOmega ,Num_of_Elem.STP    );
D3 = sparse(Num_of_Elem.STV     ,Num_of_Elem.STOmega);

if size(find(SpElemProperties.SpS.Belong_to_ST_FI),1)==0
    return
end

%% D0
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    for IncSpNIdx = find(sG(SpSIdx,:))
        UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(SpSIdx);
        for  CurrentTimeSection = -1
            STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx);
            STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx);
            D0(STS_tar,STN_tar) = sG(SpSIdx,IncSpNIdx);
        end
        for CurrentTimeSection = 0:SpElemProperties.SpS.UpdNum(SpSIdx)
            STS_tar = SpElemProperties.SpS.FirstSTSIdx(SpSIdx)+1+CurrentTimeSection;
            STN_tar = SpElemProperties.SpN.FirstSTNIdx(IncSpNIdx)+1+UpdRatio*CurrentTimeSection;
            D0(STS_tar,STN_tar) = sG(SpSIdx,IncSpNIdx);
        end
    end
end
for SpNIdx = find(SpElemProperties.SpN.Belong_to_ST_FI)
    for CurrentTimeSection = 0:SpElemProperties.SpN.UpdNum(SpNIdx)
        STS_tar = SpElemProperties.SpN.FirstSTSIdx(SpNIdx)+CurrentTimeSection;
        STN_tar = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+CurrentTimeSection;
        D0(STS_tar,STN_tar) = -1;
        STN_tar = SpElemProperties.SpN.FirstSTNIdx(SpNIdx)+CurrentTimeSection+1;
        D0(STS_tar,STN_tar) =  1;
    end
end

%% D1
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpSIdx = find(sC(SpPIdx,:))
        UpdRatio = SpElemProperties.SpS.UpdNum(IncSpSIdx)/SpElemProperties.SpP.UpdNum(SpPIdx);
        for CurrentTimeSection = 0:SpElemProperties.SpP.UpdNum(SpPIdx)
            STP_tar = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTimeSection;
            STS_tar = SpElemProperties.SpS.FirstSTSIdx(IncSpNIdx)+1+UpdRatio*CurrentTimeSection;
            D1(STP_tar,STS_tar) = sC(SpPIdx,IncSpSIdx);
        end
    end
end
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    for IncSpNIdx = find(sG(SpSIdx,:))
        UpdRatio = SpElemProperties.SpN.UpdNum(IncSpNIdx)/SpElemProperties.SpS.UpdNum(SpSIdx);
        for CurrentTimeSection = 0
            STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
            STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx);
            D1(STP_tar,STS_tar) = sG(SpSIdx,IncSpNIdx);
        end
        for CurrentTimeSection = 1:SpElemProperties.SpS.UpdNum(SpSIdx)
            for LocalTimeSec_for_SpN = 1:UpdRatio
                STP_tar = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+CurrentTimeSection;
                STS_tar = SpElemProperties.SpN.FirstSTSIdx(IncSpNIdx)+LocalTimeSec_for_SpN+UpdRatio*(CurrentTimeSection-1);
                D1(STP_tar,STS_tar) = sG(SpSIdx,IncSpNIdx);
            end
        end
    end
end

%% D2
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpVIdx = find(sD(:,SpPIdx))
        UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(IncSpVIdx);
        for CurrentTime = 0:SpElemProperties.SpV.UpdNum(IncSpVIdx)
            STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx)+CurrentTime;
            STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+UpdRatio*CurrentTime;
            D2(STOmega_tar,STP_tar) = sD(IncSpVIdx,SpPIdx);
        end
    end
    for CurrentTimeSection = 0
        STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpSIdx);
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx);
        D2(STOmega_tar,STP_tar) =  1;
    end
    for CurrentTimeSection = 1:SpElemProperties.SpP.UpdNum(SpSIdx)
        STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpSIdx)+CurrentTimeSection;
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTimeSection-1;
        D2(STOmega_tar,STP_tar) = -1;
        STP_tar     = SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+CurrentTimeSection;
        D2(STOmega_tar,STP_tar) =  1;
    end
end
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    for IncSpPIdx = find(sC(:,SpSIdx))
        UpdRatio = SpElemProperties.SpS.UpdNum(SpSIdx)/SpElemProperties.SpP.UpdNum(IncSpPIdx);
        for CurrentTimeSection = 0
            STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(IncSpPIdx);
            STP_tar     = SpElemProperties.SpS.FirstSTPIdx(SpSIdx);
            D2(STOmega_tar,STP_tar) =  sC(IncSpPIdx,SpSIdx);
        end
        for CurrentTimeSection = 1:SpElemProperties.SpP.UpdNum(SpSIdx)
            for LocalTimeSec_for_SpS = 1:UpdRatio
                STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(IncSpPIdx)+CurrentTimeSection;
                STP_tar     = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+LocalTimeSec_for_SpS+UpdRatio*(CurrentTimeSection-1);
                D2(STOmega_tar,STP_tar) = sC(IncSpPIdx,SpSIdx);
            end
        end
    end
end

%% D3
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    for IncSpVIdx = find(sD(:,SpPIdx))
        UpdRatio = SpElemProperties.SpP.UpdNum(SpPIdx)/SpElemProperties.SpV.UpdNum(IncSpVIdx);
        for CurrentTimeSection = 0
             STV_tar     = SpElemProperties.SpV.FirstSTVIdx(IncSpVIdx);
             STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx);
             D3(STV_tar,STOmega_tar) =  1;
             for LocalTimeSec_for_SpP = 1:UpdRatio
                 STOmega_tar = SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx);
                 D3(STV_tar,STOmega_tar) = sD(IncSpVIdx,SpPIdx);
             end
        end
        for CurrentTimeSection = 1:SpElemProperties.SpV.UpdNum(IncSpVIdx)
            STV_tar     = SpElemProperties.SpV.FirstSTVIdx(IncSpVIdx)+CurrentTimeSection;
            STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx)+CurrentTimeSection-1;
            D3(STV_tar,STOmega_tar) = -1;
            STOmega_tar = SpElemProperties.SpV.FirstSTOmegaIdx(IncSpVIdx)+CurrentTimeSection;
            D3(STV_tar,STOmega_tar) =  1;
            for LocalTimeSec_for_SpP = 1:UpdRatio
                STOmega_tar = ...
                    SpElemProperties.SpP.FirstSTOmegaIdx(SpPIdx)+LocalTimeSec_for_SpP+UpdRatio*CurrentTimeSection;
                D3(STV_tar,STOmega_tar) = sD(IncSpVIdx,SpPIdx);
            end
        end
    end
end


end