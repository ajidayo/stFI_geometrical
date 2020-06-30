function [kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem)
kappa = zeros(Num_of_Elem.STP);
[kappa,FaceArea.Dual]   = ComputeKappa_for_SpFI_SpSs(kappa,cdt,sG,sC,sD,NodePos,SpElemProperties);
[kappa,FaceArea.Prim]   = ComputeKappa_for_SpFI_SpPs(kappa,cdt,sG,sC,sD,NodePos,SpElemProperties);
end
%%
function [kappa,FaceAreaPrim] = ComputeKappa_for_SpFI_SpPs(kappa,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI)
    FaceAreaPrim = zeros(Num_of_Elem.SpP);
    Area_PrimSTP = SpFaceArea(SpPIdx,sG,sC,NodePos.Prim);
    Area_DualSTP = SpEdgeLength(SpPIdx,sD.',NodePos.Dual)*cdt/SpElemProperties.SpP.UpdNum;
    for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpPIdx):...
            SpElemProperties.SpS.FirstSTPIdx(SpPIdx)+SpElemProperties.SpS.UpdNum
        kappa(STPIdx) = Area_DualSTP/Area_PrimSTP;
    end
    FaceAreaPrim(SpPIdx) = Area_PrimSTP;
end
end
function [kappa,FaceAreaDual] = ComputeKappa_for_SpFI_SpSs(kappa,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI)
    FaceAreaDual = zeros(Num_of_Elem.SpS);
    Area_PrimSTP = SpEdgeLength(SpSIdx,sG,NodePos.Prim)*cdt/SpElemProperties.SpS.UpdNum;
    Area_DualSTP = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
    for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
            SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum
        kappa(STPIdx) = Area_DualSTP/Area_PrimSTP;
    end
    FaceAreaDual(SpSIdx) = Area_DualSTP;
end
end

%%
function Area_Face = SpFaceArea(SpPIdx,sG,sC,NodePos)
Nodes = find(logical(sC(SpPIdx,:))*logical(sG));
Area_Face = ThreeD_SpP_Area(Nodes,sG,NodePos);
end
function Length_Edge = SpEdgeLength(SpSIdx,sG,NodePos)
Endpoints = find(sG(SpSIdx));
Length_Edge = norm(NodePos(Endpoints(1)).Vec - NodePos(Endpoints(1)).Vec);
end
function Area_3DFace = ThreeD_SpP_Area(Nodes,sG,NodePos)
StartNode = Nodes(1);
Nodes(1) = [];
LastNode = StartNode;
IsAdj_to_StartNode = logical(sG).'*logical(sG(:,LastNode));
for NodeIdx = Nodes
    if IsAdj_to_StartNode(NodeIdx)==true
        NextNode = NodeIdx;
        Nodes(2) = [];
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
    TriangleArea_Signed = norm(cross(Vec1,Vec2));
    FaceArea_Signed = FaceArea_Signed + TriangleArea_Signed;
end
Area_3DFace = abs(FaceArea_Signed);
end
