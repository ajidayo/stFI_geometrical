function [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties] ...
    = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum)
switch RefMeshPresetType
    case 'FDTD'
        [sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties] ...
            = FDTD_RefMesh(MeshMeasurements,LocalUpdateNum);
    case 'BeltLike_LocalTimeStep_NoRefinement'
   
    case 'BeltLike_LocalTimeStep_LocalRefinement'
        
    otherwise
        warning('Unexpected RefMeshPresetType. No reference mesh generated.')
end
end