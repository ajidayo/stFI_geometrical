function [MeshParam]=Parameters_Mesh(MeshMeasurements)

% MeshMeasurements.XCoord
% MeshMeasurements.YCoord
% MeshMeasurements.FineStartAtYCoord=66;
% MeshMeasurements.TriangleStartAtYCoord=67;
% MeshMeasurements.TriangleEndAtYCoord=85;
% MeshMeasurements.FineEndAtYCoord=86;


MeshParam.Size_X=MeshMeasurements.XCoord;

MeshParam.Fine_Y_from=floor(MeshMeasurements.FineStartAtYCoord)+1;

MeshParam.Triangle_Y_from=...
    2*( (1/2)*floor(2*MeshMeasurements.TriangleStartAtYCoord)...
    -floor(MeshMeasurements.FineStartAtYCoord) )...
    +MeshParam.Fine_Y_from;
MeshParam.Triangle_Y_from=round(MeshParam.Triangle_Y_from);
MeshParam.Triangle_Y_to=...
    2*( (1/2)*ceil(2*MeshMeasurements.TriangleEndAtYCoord)...
    - (1/2)*floor(2*MeshMeasurements.TriangleStartAtYCoord) )...
    +MeshParam.Triangle_Y_from...
    -1;
MeshParam.Triangle_Y_to=round(MeshParam.Triangle_Y_to);
MeshParam.Fine_Y_to=...
    2*( (1/2)*ceil(2*MeshMeasurements.FineEndAtYCoord)...
    - (1/2)*ceil(2*MeshMeasurements.TriangleEndAtYCoord) )...
    +MeshParam.Triangle_Y_to;
MeshParam.Fine_Y_to=round(MeshParam.Fine_Y_to);
MeshParam.Size_Y=...
    round(MeshMeasurements.YCoord)...
    -(1/2)*ceil(2*MeshMeasurements.FineEndAtYCoord)...
    +MeshParam.Fine_Y_to;
MeshParam.Size_Y=round(MeshParam.Size_Y);

% MeshParam.Fine_Y_from=66;
% MeshParam.Triangle_Y_from=68;
% MeshParam.Triangle_Y_to=103;
% MeshParam.Fine_Y_to=105;
% MeshParam.Size_Y=120;

% %% mesh DOF parameters
% 
% MeshParam.deltatriangle=0.1;
% MeshParam.deltaboundary=1.0/12.0;

