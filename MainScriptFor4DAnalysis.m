clear;

%% Inputs
SelectPreset=1;% Preset {none} available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] ...
                            = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties] ...
                            = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
                        
% Usage: NodePos(SpSIdx).Vec
% RefImpedance_SpV           = GenerateReferenceImpedancePattern();
% Source =
% Source().UpdNum 
% Source().DualFace_tgt
% Source().WaveformFunctionHandle;
% Source().Area_TargetDualFace;
cdt                         = 0.5;
disp('point1')
%%
[SpElemProperties,Num_of_Elem.STP] ...
    = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem);
disp('point2')
Task                        = struct;
TaskDepGraph = digraph;
Map_SpElem_to_FirstGlobTask = struct;
%[Task,TaskDepGraph] = GenerateST_FI_Tasks_4D_ST;
[Task,TaskDepGraph,SpElemProperties,Map_SpElem_to_FirstGlobTask] ...
                            = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,Num_of_Elem,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);
TaskOrder                   = SortTasks(TaskDepGraph,Map_SpElem_to_FirstGlobTask,sC,SpElemProperties);
clearvars TaskDepGraph;
disp('point3')
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem);
disp('point4')
[kappa,FaceArea]            = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem);
Z                           = ComputeImpedance_for_EachSTPs(RefImpedance_SpV,sC,sD,SpElemProperties,Num_of_Elem);
Kappa_over_Z                = kappa/Z;
disp('point5')
TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Source);
TMM_Fields                  = ExcludePEC_ST_Planes(TMM_Fields, SpElemProperties);
DoFs_FacesThenEdges         = InitializeFields;
Num_of_Steps                = 100;
Time                        = 0;
disp('point6')
DoFs_FacesThenEdges         = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,DoFs_FacesThenEdges,Source);
%PlotMagneticFlux