function [kappa,FaceArea] = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem)
global SpDIM EPSILON
kappa = zeros(Num_of_Elem.STP,1);

FaceArea.Prim = zeros(Num_of_Elem.SpP,1);
FaceArea.Dual = zeros(Num_of_Elem.SpS,1);
DeltaPrimNodePos  = sparse(SpDIM, Num_of_Elem.STN);
disp('call ComputeKappa_for_SpFI_SpSs')
[kappa,FaceArea.Dual]   = ComputeKappa_for_SpFI_SpSs(kappa,FaceArea.Prim,cdt,sG,sC,sD,NodePos,SpElemProperties);
disp('call ComputeKappa_for_SpFI_SpPs')
[kappa,FaceArea.Prim]   = ComputeKappa_for_SpFI_SpPs(kappa,FaceArea.Dual,cdt,sG,sC,sD,NodePos,SpElemProperties);
end
%%
function [kappa,FaceArea] = ComputeKappa_for_STFI_SpPs_and_SpSs(kappa,FaceArea,cdt,sG,sC,sD,NodePos,SpElemProperties)
for SpSIdx = find(SpElemProperties.SpS.UpdNumBoundary==true)
    if SpElemProperties.SpS.PEC(SpSIdx) == true
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            kappa(STPIdx) = 0;
        end
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    elseif SpElemProperties.SpS.UpdNumCorner == true
        continue;
    else
        PrimEdgeTgtLeng = SpEdgeLength(SpSIdx,sG,NodePos.Prim);
        DualFaceTgtArea = SpFaceArea(SpSIdx,sD.',sC.',NodePos.Dual);
        Area_PrimSTP = PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx);
        Area_DualSTP = DualFaceTgtArea;
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum
            kappa(STPIdx) = -1*Area_DualSTP/Area_PrimSTP;
        end
        FaceArea.Dual(SpSIdx) = Area_DualSTP;
        DeltaSTNPos = ...
            Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt);
    end 
end

for SpSIdx = find(SpElemProperties.SpS.Belong_to_ST_FI==true)
    if SpElemProperties.SpS.UpdNumBoundary(SpSIdx)==true
        continue;
    elseif SpElemProperties.SpS.PEC(SpSIdx) == true
        for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx):...
                SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
            kappa(STPIdx) = 0;
        end
        FaceArea.Dual(SpSIdx) = 1;
        continue;
    end
    [kappa, FaceArea.Dual] = ...
        ComputeKappa_for_SingleNonBoundarySpP(kappa,SpSIdx,sG,sC,sD,NodePos,SpElemProperties);
end

for SpPIdx = find(SpElemProperties.SpP.Belong_to_ST_FI == true)
    if SpElemProperties.SpP.ElecWall(SpPIdx) == true
        for STPIdx = SpElemProperties.SpP.FirstSTPIdx(SpPIdx):...
                SpElemProperties.SpP.FirstSTPIdx(SpPIdx)+SpElemProperties.SpP.UpdNum(SpPIdx)
            kappa(STPIdx) = 0;
        end
        FaceArea.Prim(SpPIdx) = 1;
        continue;
    end
    
end
end
%%
function DeltaSTNPos = Compute_DeltaSTNPos_for_SingleBoundaryEdge(SpSIdx,DeltaSTNPos,NodePos,PrimEdgeTgtLeng,DualFaceTgtArea,sG,sC,D2,D3,SpElemProperties,STElemProperties,cdt)
kappa       = -1*DualFaceTgtArea/(PrimEdgeTgtLeng*cdt/SpElemProperties.SpS.UpdNum(SpSIdx));
LocalBasis1 = findDeltaDirection(SpSIdx,sG,sC,NodePos.Prim,SpElemProperties);
SpS_Direc   = sG(SpSIdx,:)*NodePos.Prim;
LocalBasis2 = SpS_Direc - dot(SpS_Direc,LocalBasis1)*LocalBasis1;
LocalBasis2 = LocalBasis2/norm(LocalBasis2);
LocalBasis3 = cross(LocalBasis1,LocalBasis2);

SpEndPoints = find(sG(SpSIdx,:));
STNIdx_Past = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints);
STNIdx_Futr = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)+1;
for STPIdx = SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+1:...
        SpElemProperties.SpS.FirstSTPIdx(SpSIdx)+SpElemProperties.SpS.UpdNum(SpSIdx)
    STNIdx_Past = STNIdx_Past+1;
    STNIdx_Futr = STNIdx_Futr+1;
    
    STNTildeIdx_DefiningDualFace = find(logical(D2(:,STPIdx)).'*logical(D3).');
    STNTildeIdx_DefiningDualFace = SortSTNs_Along_with_FaceOrientation(STNTildeIdx_DefiningDualFace,STPIdx,D3   ,D2   );
    %                              SortSTNs_Along_with_FaceOrientation(STNIdx_DefiningFace         ,STPTar,GradM,CurlM)
    PosVec_Projected_to_03Plane = zeros(2,size(STNTildeIdx_DefiningDualFace,2));
    ColIdx = 0;
    for STNTildeIdx = STNTildeIdx_DefiningDualFace 
        ColIdx = ColIdx +1;
        PosVec_Projected_to_03Plane(0,ColIdx) = cdt *STElemProperties.STV.TimeIdx(STNTildeIdx);
        SpNTildeIdxTemp                       =      STElemProperties.STV.RefSpV(STNTildeIdx);
        PosVec_Projected_to_03Plane(1,ColIdx) = dot(NodePos.Dual(SpNTildeIdxTemp), LocalBasis3)*LocalBasis3;
    end
    Area03_STPtilde       =   SignedArea2D(PosVec_Projected_to_03Plane);
    Delta_of_DeltaNodePos = -(Area03_STPtilde/(PrimEdgeTgtLeng*dot(SpS_Direc, LocalBasis2)))/kappa;
    DeltaSTNPos(:,STNIdx_Futr) = DeltaSTNPos(:,STNIdx_Past) ...
        + Delta_of_DeltaNodePos*LocalBasis1;
end
STNIdx_Last  = STNIdx_Futr;
STNIdx_First = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints);
DeltaSTNPos(STNIdx_First) = DeltaSTNPos(STNIdx_Last);
end
%%
function DeltaDirection = findDeltaDirection(SpSTgt,sG,sC,NodePosPrim,SpElemProperties)
global EPSILON
IncFaceIdx = find(sC(:,SpSTgt)).';
for FaceIdx = IncFaceIdx
    GuessedDirection = [0;0;0];
    IncEdgeIdx = find(sC(FaceIdx,:));
    for EdgeIdx = IncEdgeIdx
        if SpElemProperties.SpS.is_on_UpdNumBoundary == true
            continue;
        elseif EdgeIdx ~=SpSTgt && logical(sG(SpSTgt,:))*logical(sG(EdgeIdx,:)).' ~= 0
            if norm(GuessedDirection) == 0 
                GuessedDirection = sG(EdgeIdx,:)*NodePosPrim;
                GuessedDirection = GuessedDirection/norm(GuessedDirection);
            elseif norm(cross(GuessedDirection,sG(EdgeIdx,:)*NodePosPrim))<EPSILON
                DeltaDirection = GuessedDirection;
                return;
            else
                break;
            end
        end
    end
end
end

%%
function SortedSTNIdx = SortSTNs_Along_with_FaceOrientation(STNIdx_DefiningFace,STPTar,GradM,CurlM)
InputSTNNum = size(STNIdx_DefiningFace,2);
SortedSTNIdx = zeros(1,InputSTNNum);
LeftSTN = STNIdx_DefiningFace;
SortedSTNIdx(1) = STNIdx_DefiningFace(1);
LeftSTN(1) = [];
ProgressIdx=1;
for Counter = 2:InputSTNNum
    for STSIdx = find(CurlM(STPTar,:)).'
        if GradM(STSIdx,SortedSTNIdx(ProgressIdx)) ~= 0 && GradM(STSIdx,SortedSTNIdx(ProgressIdx))*CurlM(STPTar,STSIdx)==-1
            GuideSTSIdx = STSIdx;
            break;
        end
    end
    for LeftSTNScan = 1:size(LeftSTN,2)
        GuessedSTNIdx = LeftSTN(LeftSTNScan);
        if GradM(GuideSTSIdx,GuessedSTNIdx)~=0
            ProgressIdx = ProgressIdx+1;
            SortedSTNIdx(ProgressIdx) = GuessedSTNIdx;
            LeftSTN(LeftSTNScan) = [];
            break;
        end
    end
end
end
%%
function AreaReturn = SignedArea2D(NodePos2D)
Pos0 = NodePos2D(:,1);
Vec1 = zeros('like',Pos0);
Vec2 = NodePos2D(:,2) - NodePos2D(:,1);
NodePos2D(:,1:2) = [];
Num_of_TriangleDivision = size(NodePos2D,2);
AreaReturn = 0;
for TriDivIdx = 1:Num_of_TriangleDivision
    Vec1 = Vec1 + Vec2;
    Vec2 = NodePos2D(:,1) - Pos0;
    NodePos2D(:,1) = [];
    TriangleArea_Signed = 0.5 * Vec1.'*[0 1;-1 0]*Vec2;
    AreaReturn = AreaReturn + TriangleArea_Signed;
end
end

%%
function [kappa, FaceArea_Dual] = ComputeKappa_for_SingleNonBoundarySpS(kappa,SpSTgt,sG,sC,sD,NodePos,SpElemProperties)
PrimEdgeTgtLeng         = SpEdgeLength(SpSTgt,sG,NodePos.Prim);
DualFaceTgtArea         = SpFaceArea(SpSTgt,sD.',sC.',NodePos.Dual);
FaceArea_Dual(SpSTgt)   = DualFaceTgtArea;

SpEndPoints         = find(sG(SpSTgt,:));
UpdNum_SpSTgt       = SpElemProperties.SpS.UpdNum(SpSTgt);
UpdNumEP            = SpElemProperties.SpN.UpdNum(SpEndPoints);

STNIdx_forEPs_Past  = SpElemProperties.SpN.FirstSTNIdx(SpEndPoints)-UpdNumEP/UpdNum_SpSTgt+1;
STPTgt              = SpElemProperties.SpS.FirstSTPIdx(SpSTgt);
for TimeSecTgt=1:UpdNum_SpSTgt
    STPTgt = STPTgt + 1;
    STNIdx_forEPs_Past  = STNIdx_forEPs_Past + UpdNumEP/UpdNum_SpSTgt;
    
    BulgeArea_EP1Side = ...
        CalcBulgeArea(STNIdx_forEPs_Past(1),SpEndPoints(1),SpSTgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
    BulgeArea_EP2Side = ...
        CalcBulgeArea(STNIdx_forEPs_Past(2),SpEndPoints(2),SpSTgt,DeltaSTNPos,NodePos.Prim,sG,SpElemProperties,cdt);
   
    Area_PrimSTP  = ...
        PrimEdgeTgtLeng *cdt/UpdNum_SpSTgt + (BulgeArea_EP1Side+BulgeArea_EP2Side);
    Area_DualSTP  = DualFaceTgtArea;
    
    kappa(STPTgt) = - Area_DualSTP/Area_PrimSTP;    
end
kappa(STPTgt-UpdNum_SpSTgt) = kappa(STPTgt);
end

%%
function BulgeArea = ...
    CalcBulgeArea(STNTgt_Oldest,SpNTgt,SpSTgt,DeltaSTNPos,NodePosPrim,sG,SpElemProperties,cdt)
Basis1 = sG(SpSTgt,:) * NodePosPrim;
Basis1 = Basis1/norm(Basis1);

BulgeArea = 0;
for LocalTimeSec = 1:SpElemProperties.SpN(SpNTgt)/SpElemProperties.SpS.UpdNum(SpSTgt)
    EdgeLengPast = ...
        sG(SpSTgt,SpNTgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec-1), Basis1);
    EdgeLengFutr = ...
        sG(SpSTgt,SpNTgt) * dot(DeltaSTNPos(:,STNTgt_Oldest+LocalTimeSec  ), Basis1);
    BulgeArea = BulgeArea ...
        + (cdt/SpElemProperties.SpN(SpNTgt)) * (EdgeLengPast + EdgeLengFutr)/2;
end
end

%%
function [kappa, FaceAreaPrim] = ...
    ComputeKappa_for_SingleNonBoundarySpP(kappa,SpPTgt,sG,sC,D2,D3,NodePos,SpElemProperties)

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
