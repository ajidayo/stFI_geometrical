function [kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem)
kappa = zeros(Num_of_Elem.STP,1);

FaceArea.Prim = zeros(Num_of_Elem.SpP,1);
FaceArea.Dual = zeros(Num_of_Elem.SpS,1);
disp('call ComputeKappa_for_SpFI_SpSs')
[kappa,FaceArea.Dual]   = ComputeKappa_for_SpFI_SpSs(kappa,FaceArea.Prim,cdt,sG,sC,sD,NodePos,SpElemProperties);
disp('call ComputeKappa_for_SpFI_SpPs')
[kappa,FaceArea.Prim]   = ComputeKappa_for_SpFI_SpPs(kappa,FaceArea.Dual,cdt,sG,sC,sD,NodePos,SpElemProperties);
end
%%
function [kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==true)
    switch SpElemProperties.SpS.PEC(SpSIdx)
        case true
            for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                    SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
                kappa(STPIdx) = 0;
            end
            FaceArea.Dual(SpSIdx) = 1;
        case false
            %
            %
    end
end
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI==true)
    switch SpElemProperties.SpP.ElecWall(SpPIdx) 
        case true
            for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                    SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
                kappa(STPIdx) = 0;
            end
            FaceArea.Prim(SpPIdx) = 1;
        case false
            %
            %
    end
end
end

%%
function [kappa,FaceAreaPrim] = ComputeKappa_for_SpFI_SpPs(kappa,FaceAreaPrim,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI==false)
    switch SpElemProperties.SpP.ElecWall(SpPIdx)
        case true
            for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                    SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
                kappa(STPIdx) = 0;
            end
            FaceAreaPrim(SpPIdx) = 1;
        case false
            if SpElemProperties.SpP.PrimAreaIsGiven(SpPIdx) == true && SpElemProperties.SpP.DualLengIsGiven(SpPIdx) == true
                Area_PrimSTP = SpElemProperties.SpP.PrimArea(SpPIdx);
                Area_DualSTP = SpElemProperties.SpP.DualLeng(SpPIdx)*cdt/SpElemProperties.SpP.UpdNum(SpPIdx);
            else
                Area_PrimSTP = SpFaceArea(SpPIdx,sG,sC,NodePos.Prim);
                Area_DualSTP = SpEdgeLength(SpPIdx,sD.',NodePos.Dual)*cdt/SpElemProperties.SpP.UpdNum(SpPIdx);
            end
            for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                    SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
                kappa(STPIdx) = Area_DualSTP/Area_PrimSTP;
            end
            FaceAreaPrim(SpPIdx) = Area_PrimSTP;
    end
end
end
function [kappa,FaceAreaDual] = ComputeKappa_for_SpFI_SpSs(kappa,FaceAreaDual,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==false)
    switch SpElemProperties.SpS.PEC(SpSIdx)
        case true
            for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                    SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
                kappa(STPIdx) = 0;
            end
            FaceAreaDual(SpSIdx) = 1;
        case false
            if SpElemProperties.SpS.PrimLengIsGiven(SpSIdx) == true && SpElemProperties.SpS.DualAreaIsGiven(SpSIdx) == true  
                Area_PrimSTP = SpElemProperties.SpS.PrimLeng(SpSIdx)*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
                Area_DualSTP = SpElemProperties.SpS.DualArea(SpSIdx);               
            else
                Area_PrimSTP = SpEdgeLength(SpSIdx,sG,NodePos.Prim)*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
                Area_DualSTP = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
            end
            for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                    SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum
                kappa(STPIdx) = -1*Area_DualSTP/Area_PrimSTP;
            end
            FaceAreaDual(SpSIdx) = Area_DualSTP;
    end
end
end

%%
function Area_Face = SpFaceArea(SpPIdx,sG,sC,NodePos)
Nodes = find(logical(sC(SpPIdx,:))*logical(sG));
Area_Face = ThreeD_SpP_Area(Nodes,sG,NodePos);
end
function Length_Edge = SpEdgeLength(SpSIdx,sG,NodePos)
Endpoints = find(sG(SpSIdx,:));
Length_Edge = norm(NodePos(Endpoints(2)).Vec - NodePos(Endpoints(1)).Vec);
end
function Area_3DFace = ThreeD_SpP_Area(Nodes,sG,NodePos)
StartNode = Nodes(1);
Nodes(1) = [];
LastNode = StartNode;
IsAdj_to_StartNode = logical(sG).'*logical(sG(:,LastNode));
for LeftNodeScan = 1:size(Nodes,2)
    NodeIdx = Nodes(LeftNodeScan);
    if IsAdj_to_StartNode(NodeIdx)==true
        NextNode = NodeIdx;
        Nodes(LeftNodeScan) = [];
        break;
    else
        NextNode = 0;
    end
end
clearvars IsAdj_to_StartNode
Vec2 = NodePos(NextNode).Vec - NodePos(LastNode).Vec;
Vec1 = zeros('like',Vec2);
Num_of_TriangleDivision = size(Nodes,2);
FaceArea_Signed = 0;
for TriDivIdx = 1:Num_of_TriangleDivision
    IsAdj_to_NextNode = logical(sG).'*logical(sG(:,NextNode));
    for LeftNodeScan = 1:size(Nodes,2)
        NodeIdx = Nodes(LeftNodeScan);
        if IsAdj_to_NextNode(NodeIdx)==true && NodeIdx ~= LastNode
            LastNode = NextNode;
            NextNode = NodeIdx;
            Vec1 = Vec1+Vec2;
            Vec2 = NodePos(NextNode).Vec - NodePos(LastNode).Vec;
            Nodes(LeftNodeScan) = [];
            break;
        end
    end
    TriangleArea_Signed = 0.5*norm(cross(Vec1,Vec2));
    FaceArea_Signed = FaceArea_Signed + TriangleArea_Signed;
end
Area_3DFace = abs(FaceArea_Signed);
end
