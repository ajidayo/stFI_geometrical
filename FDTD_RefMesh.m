function [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties_SpV_UpdNum] ...
    = FDTD_RefMesh(MeshMeasurements,LocalUpdateNum)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

Num_of_Elem.SpV = XSize*YSize*ZSize;
Num_of_Elem.SpP = ...
    (XSize+1)*YSize*ZSize ...
    + XSize*(YSize+1)*ZSize ...
    + XSize*YSize*(ZSize+1);
Num_of_Elem.SpS ...
    =(XSize  )*(YSize+1)*(ZSize+1) ...
    +(XSize+1)*(YSize  )*(ZSize+1) ...
    +(XSize+1)*(YSize+1)*(ZSize  );
Num_of_Elem.SpN = (XSize+1)*(YSize+1)*(ZSize+1);

sG = FDTD_sG(MeshMeasurements,Num_of_Elem);
sC = FDTD_sC(MeshMeasurements,Num_of_Elem);
sD = FDTD_sD(MeshMeasurements,Num_of_Elem);
NodePos = FDTD_NodePos(MeshMeasurements);
SpElemProperties_SpV_UpdNum = LocalUpdateNum*ones(Num_of_Elem.SpV,1);
end

function sG = FDTD_sG(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = XSize+1;
YEdgePerXYPlane     = (XSize+1)*YSize;
YEdgeNum            = (XSize+1)*YSize*(ZSize+1);
ZEdgePerXRow        = XSize+1;
ZEdgePerXYPlane     = (XSize+1)*(YSize+1);
%ZEdgeNum            = (XSize+1)*(YSize+1)*ZSize;

NodePerXRow        = XSize+1;
NodePerXYPlane     = (XSize+1)*(YSize+1);
%NodeNum            = (XSize+1)*(YSize+1)*(ZSize+1);

sG = sparse(Num_of_Elem.SpS,Num_of_Elem.SpN);
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize
            SIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx+1 + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize+1
            SIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx+1 + (YIdx  )*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize+1
            SIdx = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            IncNIdx = XIdx   + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            sG(SIdx,IncNIdx) = -1;
            IncNIdx = XIdx+1 + (YIdx-1)*NodePerXRow + (ZIdx  )*NodePerXYPlane;
            sG(SIdx,IncNIdx) =  1;
        end
    end
end

end

function sC = FDTD_sC(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;
YZFacePerXRow       = XSize+1;
YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum           = (XSize+1)*YSize*ZSize;
ZXFacePerXRow       = XSize;
ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum           = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum          = XSize*YSize*(ZSize+1);

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = XSize+1;
YEdgePerXYPlane     = (XSize+1)*YSize;
YEdgeNum            = (XSize+1)*YSize*(ZSize+1);
ZEdgePerXRow        = XSize+1;
ZEdgePerXYPlane     = (XSize+1)*(YSize+1);
%ZEdgeNum            = (XSize+1)*(YSize+1)*ZSize;

sC = sparse(Num_of_Elem.SpP,Num_of_Elem.SpS);

for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize+1
            PIdx = XIdx   + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx  )*YEdgePerXYPlane+XEdgeNum;  
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx  )*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize
            PIdx = XIdx   + (YIdx-1)*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx  )*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx+1 + (YIdx-1)*ZEdgePerXRow + (ZIdx-1)*ZEdgePerXYPlane+XEdgeNum+YEdgeNum;
            sC(PIdx,IncSIdx) =  1;
        end
    end
end
for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            PIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx-1)*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            IncSIdx = XIdx   + (YIdx-1)*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx   + (YIdx  )*XEdgePerXRow + (ZIdx-1)*XEdgePerXYPlane;
            sC(PIdx,IncSIdx) =  1;
            IncSIdx = XIdx   + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;  
            sC(PIdx,IncSIdx) = -1;
            IncSIdx = XIdx+1 + (YIdx-1)*YEdgePerXRow + (ZIdx-1)*YEdgePerXYPlane+XEdgeNum;
            sC(PIdx,IncSIdx) =  1;
        end
    end
end
end
function sD = FDTD_sD(MeshMeasurements,Num_of_Elem)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
%VolNum              = XSize*YSize*ZSize;

YZFacePerXRow       = XSize+1;
YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum        = (XSize+1)*YSize*ZSize;
ZXFacePerXRow       = XSize;
ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum        = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum        = XSize*YSize*(ZSize+1);

sD = sparse(Num_of_Elem.SpV,Num_of_Elem.SpP);
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            VIdx = XIdx + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
            IncPIdx = XIdx   + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx+1 + (YIdx-1)*YZFacePerXRow + (ZIdx-1)*YZFacePerXYPlane;
            sD(VIdx,IncPIdx) =  1;
            IncPIdx = XIdx   + (YIdx-1)*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx   + (YIdx  )*ZXFacePerXRow + (ZIdx-1)*ZXFacePerXYPlane ...
                +YZFaceNum;            
            sD(VIdx,IncPIdx) =  1;
            IncPIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx-1)*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            sD(VIdx,IncPIdx) = -1;
            IncPIdx = XIdx   + (YIdx-1)*XYFacePerXRow + (ZIdx  )*XYFacePerXYPlane ...
                +YZFaceNum+ZXFaceNum;
            sD(VIdx,IncPIdx) =  1;
        end
    end
end
end

function NodePos = FDTD_NodePos(MeshMeasurements)
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

NodePerXRow        = XSize+1;
NodePerXYPlane     = (XSize+1)*(YSize+1);
%NodeNum            = (XSize+1)*(YSize+1)*(ZSize+1);
VolPerXRow          = XSize;
VolPerXYPlane       = XSize*YSize;
%VolNum              = XSize*YSize*ZSize;



for ZIdx = 1:ZSize+1
    for YIdx = 1:YSize+1
        for XIdx = 1:XSize+1
            SpNIdx = XIdx   + (YIdx-1)*NodePerXRow + (ZIdx-1)*NodePerXYPlane;
            XPos = MeshMeasurements.dx*(XIdx-1);
            YPos = MeshMeasurements.dy*(YIdx-1);
            ZPos = MeshMeasurements.dz*(ZIdx-1);
            NodePos.Prim(SpNIdx).Vec = [XPos;YPos;ZPos]; 
        end
    end
end
for ZIdx = 1:ZSize
    for YIdx = 1:YSize
        for XIdx = 1:XSize
            SpNIdx = XIdx   + (YIdx-1)*VolPerXRow + (ZIdx-1)*VolPerXYPlane;
            XPos = MeshMeasurements.dx*(XIdx-0.5);
            YPos = MeshMeasurements.dy*(YIdx-0.5);
            ZPos = MeshMeasurements.dz*(ZIdx-0.5);
            NodePos.Dual(SpNIdx).Vec = [XPos;YPos;ZPos];
        end
    end
end
end
