function [MeshParam]=MeshParameters_square_belt(MeshMeasurements)

% MeshMeasurements.XCoord
% MeshMeasurements.YCoord
% MeshMeasurements.FineStartAtYCoord=66;
% MeshMeasurements.TriangleStartAtYCoord=67;
% MeshMeasurements.TriangleEndAtYCoord=85;
% MeshMeasurements.FineEndAtYCoord=86;


MeshParam.Size_X=MeshMeasurements.XCoord;

MeshParam.Fine_Y_from=floor(MeshMeasurements.FineStartAtYCoord)+1;
MeshParam.Fine_Y_to=...
    2*( ceil(MeshMeasurements.FineEndAtYCoord)...
    - floor(MeshMeasurements.FineStartAtYCoord) )...
    +MeshParam.Triangle_Y_to...
    -1;
MeshParam.Fine_Y_to=round(MeshParam.Fine_Y_to);
MeshParam.Size_Y=...
    round(MeshMeasurements.YCoord)...
    -ceil(MeshMeasurements.FineEndAtYCoord)...
    +MeshParam.Fine_Y_to;
MeshParam.Size_Y=round(MeshParam.Size_Y);

% MeshParam.Fine_Y_from=66;
% MeshParam.Fine_Y_to=105;
% MeshParam.Size_Y=120;
