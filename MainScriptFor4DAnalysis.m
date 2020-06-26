clear;

%% Inputs
SelectPreset=1;% Preset {1} available. See ParameterPreset for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] ...
                            = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,UpdNum,RefElemNum] ...
                            = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
%RefImpedance_SpV           = GenerateReferenceImpedancePattern();

%%
[UpdNum,Att]                = Attributes_of_Sp_Elements(sG,sC,sD,UpdNum,RefElemNum);
%[Task,TaskDependanceGraph] = GenerateST_FI_Tasks_4D_ST(sC,Att)
[Task,TaskDependanceGraph]  = GenerateSp_FI_Tasks_4D_ST(sC,Att);
TaskOrder                   = SortTasks(Task,TaskDependanceGraph);
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,UpdNum,RefElemNum);
kappa                       = ComputeKappa_4D_ST(sG,sC,sD,D0,D1,D2,D3,NodePos,UpdNum,RefElemNum,Att);
Kappa_over_Z                = ComputeZ_Matrix(kappa, RefImpedance_SpV).^(-1);
PlaceSources_in_Sp_FI_Region
ConstructTimeMarchingMatrix_4D_ST
SplitTMM_into_FieldsAndSources
ExcludePEC_ST_Planes
InitializeFields
TimeMarch
%PlotMagneticFlux