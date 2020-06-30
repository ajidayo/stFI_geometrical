clear;

%% Inputs
SelectPreset=1;% Preset {none} available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] ...
                            = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem] ...
                            = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
% Usage: NodePos(SpSIdx).Vec
%RefImpedance_SpV           = GenerateReferenceImpedancePattern();
cdt                         = 0.5;
%%
[SpElemProperties,Num_of_Elem.STP] ...
                            = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem);
Task                        = struct;
%[Task,TaskDepGraph] = GenerateST_FI_Tasks_4D_ST;
[Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask] ...
                            = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);
TaskOrder                   = SortTasks(TaskDepGraph,Map_SpElem_to_FirstGlobTask,sC,SpElemProperties);
clearvars TaskDepGraph;
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,UpdNum,Num_of_Elem);
kappa                       = ComputeKappa_4D_ST(sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem);
Kappa_over_Z                = ComputeZ_Matrix(kappa, RefImpedance_SpV).^(-1);
Source                      = PlaceSources_in_Sp_FI_Region;
TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Source);
TMM_Fields                  = ExcludePEC_ST_Planes(TMM_Fields, SpElemProperties);
DoFs_FacesThenEdges         = InitializeFields;
Num_of_Steps                = 100;
Time                        = 0;
DoFs_FacesThenEdges         = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,DoFs_FacesThenEdges,Source);
%PlotMagneticFlux