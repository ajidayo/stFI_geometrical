function [RefMeshPresetType,MeshMeasurements,LocalUpdateNum] = ParameterPreset(SelectPreset)
switch SelectPreset
    case 1
        RefMeshPresetType = 'FDTD';
        LocalUpdateNum = 1;
        MeshMeasurements.XCoord = 30;
        MeshMeasurements.YCoord = 30;
        MeshMeasurements.ZCoord = 30;
        MeshMeasurements.dx = 1;
        MeshMeasurements.dy = 1;
        MeshMeasurements.dz = 1;
        disp(['Preset: FDTD, RefMeshSize:', ...
            num2str(MeshMeasurements.XCoord), ...
            ' x ', num2str(MeshMeasurements.YCoord), ...
            ' x ', num2str(MeshMeasurements.ZCoord)])
        disp(['Discretization Width:dx, dy, dz = ', ...
            num2str(MeshMeasurements.dx), ...
            ', ', num2str(MeshMeasurements.dy), ...
            ', ', num2str(MeshMeasurements.dz)])
        disp(['LocalUpdateNum = ', num2str(LocalUpdateNum)])
    otherwise 
         warning('Preset undefined.')
end
end