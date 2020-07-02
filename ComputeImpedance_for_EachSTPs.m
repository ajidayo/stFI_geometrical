function Z = ComputeImpedance_for_EachSTPs(RefImpedance_SpV,sC,sD,SpElemProperties,Num_of_Elem)
Z = zeros(Num_of_Elem.STP,1);
for SpPIdx = 1:Num_of_Elem.SpP
    if SpElemProperties.SpP.ElecWall(SpPIdx) == true
        for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
            Z(STPIdx) = 1;
        end
    else
        ratio = 0.5;
        IncSpV = find(sD(:,SpPIdx));
        Z_SpP = (ratio*RefImpedance_SpV(IncSpV(1)).^(-1)+(1-ratio)*RefImpedance_SpV(IncSpV(2)).^(-1)).^(-1);
        for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
            Z(STPIdx) = Z_SpP;
        end
    end
end
for SpSIdx = 1:Num_of_Elem.SpS
    if SpElemProperties.SpS.PEC(SpSIdx) == true
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            Z(STPIdx) = 1;
        end
    else
        
        SurroundingSpV = find(sC(:,SpSIdx).'*sD.');
        ratio = 1.0/size(SurroundingSpV,1);
        Sum = 0;
        for SurrSpVIdx = 1:size(SurroundingSpV,1)
            Sum = ratio*RefImpedance_SpV(SurroundingSpV(SurrSpVIdx)).^(-1);
        end
        Z_SpP = Sum.^(-1);
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            Z(STPIdx) = Z_SpP;
        end
    end
end
end