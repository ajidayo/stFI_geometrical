function Source = SourceFDTD(MeshMeasurements,SpElemProperties)

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

XIdx = round(XSize/2);
YIdx = round(YSize/2);
ZIdx = round(ZSize/2);
SourceSIdx(1) = XIdx   + (YIdx-1)*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(2) = XIdx   + (YIdx  )*XEdgePerXRow  + (ZIdx-1)*XEdgePerXYPlane;
SourceSIdx(3) = XIdx   + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
SourceSIdx(4) = XIdx+1 + (YIdx-1)*YEdgePerXRow  + (ZIdx-1)*YEdgePerXYPlane +XEdgeNum;
Source = struct;
FirstST_SourceIdx = 1;
for SourceIdx = 1:4
    SpSIdx = SourceSIdx(SourceIdx);
    Source(SourceIdx).UpdNum                    = SpElemProperties.SpS.UpdNum(SpSIdx);
    Source(SourceIdx).DualFace_tgt              = SpSIdx;
    Source(SourceIdx).WaveformFunctionHandle    = @sinewave;
    Source(SourceIdx).Area_TargetDualFace       = 1;
    Source(SourceIdx).FirstST_SourceIdx         = FirstST_SourceIdx;
    FirstST_SourceIdx                           = FirstST_SourceIdx + Source(SourceIdx).UpdNum;
end
end