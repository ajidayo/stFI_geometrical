function MeshParam = MeshParameters_squarefaces_squaresubgrid(MeshMeasurements)

MeshParam.Size_X=MeshMeasurements.XCoord;
MeshParam.Fine_X_from = floor(MeshMeasurements.FineStartAtXCoord)+1;
MeshParam.Fine_X_to   =  ceil(MeshMeasurements.FineEndAtXCoord);

MeshParam.Size_Y=MeshMeasurements.YCoord;
MeshParam.Fine_Y_from = floor(MeshMeasurements.FineStartAtYCoord)+1;
MeshParam.Fine_Y_to   =  ceil(MeshMeasurements.FineEndAtYCoord);

end