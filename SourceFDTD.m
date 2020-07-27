function Source = SourceFDTD(MeshMeasurements,SpElemProperties,sC)
global SourcePeriod
SourcePeriod = 20;

disp("Waveform: Sinusoidal")
Lightspeed = 1;
wavelength = Lightspeed/SourcePeriod;
disp(["Wavelength/meshsize = ", num2str(wavelength/MeshMeasurements.dx) ])

XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

XEdgePerXRow        = XSize;
XEdgePerXYPlane     = XSize*(YSize+1);
XEdgeNum            = XSize*(YSize+1)*(ZSize+1);
YEdgePerXRow        = (XSize+1);
YEdgePerXYPlane     = (XSize+1)*YSize;
%YEdgeNum           = (XSize+1)*YSize*(ZSize+1);
%ZEdgePerXRow       = (XSize+1);
%ZEdgePerXYPlae     = (XSize+1)*(YSize+1);
%ZEdgeNum           = (XSize+1)*(YSize+1)*ZSize;

% YZFacePerXRow       = (XSize+1);
% YZFacePerXYPlane    = (XSize+1)*YSize;
YZFaceNum           = (XSize+1)*YSize*ZSize;
% ZXFacePerXRow       = XSize;
% ZXFacePerXYPlane    = XSize*(YSize+1);
ZXFaceNum           = XSize*(YSize+1)*ZSize;
XYFacePerXRow       = XSize;
XYFacePerXYPlane    = XSize*YSize;
%XYFaceNum          = XSize*YSize*(ZSize+1);

XIdx = round(XSize/2);
YIdx = round(YSize/2);
ZIdx = round(ZSize/2);
SourceSIdx(1) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(2) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(3) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(4) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SpPIdx =  XIdx   + (YIdx-1)*XYFacePerXRow  + (ZIdx-1)*XYFacePerXYPlane + YZFaceNum +ZXFaceNum;
Source = struct;
FirstST_SourceIdx = 1;
for SourceIdx = 1:4
    SpSIdx = SourceSIdx(SourceIdx);
    Source(SourceIdx).UpdNum                    = SpElemProperties.SpS.UpdNum(SpSIdx);
    Source(SourceIdx).DualFace_tgt              = SpSIdx;
    Source(SourceIdx).WaveformFunctionHandle    = @sinewave;
    %Source(SourceIdx).WaveformFunctionHandle    = @zerowave;
    Source(SourceIdx).WaveformSign              = sC(SpPIdx,SpSIdx);
    Source(SourceIdx).Area_TargetDualFace       = 1;
    Source(SourceIdx).FirstST_SourceIdx         = FirstST_SourceIdx;
    FirstST_SourceIdx                           = FirstST_SourceIdx + Source(SourceIdx).UpdNum;
end
end

function val = sinewave(x)
global SourcePeriod
    val = sin(x/SourcePeriod);
end