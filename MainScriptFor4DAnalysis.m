clear;

global SpDIM EPSILON
SpDIM   = 3;
EPSILON = 10^(-7);
%% Inputs
SelectPreset = 1;% Preset = {1} is available. See ParameterPreset for details for each settings.
[RefMeshPresetType,MeshMeasurements,LocalUpdateNum] ...
                            = ParameterPreset(SelectPreset);
[sG,sC,sD,NodePos,Num_of_Elem,SpElemProperties] ...
                            = GenerateReferenceMesh_3D_Sp(RefMeshPresetType,MeshMeasurements,LocalUpdateNum);
% RefImpedance_SpV           = GenerateReferenceImpedancePattern();
RefImpedance_SpV            = ones(Num_of_Elem.SpV,1);
cdt                         = 0.3;
disp(['cdt = ', num2str(cdt)])
disp('point1')
%%
[SpElemProperties,STElemProperties,Num_of_Elem,PrimFacePos] = Properties_of_Sp_Elements(sG,sC,sD,SpElemProperties,Num_of_Elem,NodePos);

disp('point2')
[D0,D1,D2,D3]               = ComputeST_Mesh(sG,sC,sD,SpElemProperties,Num_of_Elem);

Task                        = struct;
TaskDepGraph                = digraph;
Map_SpElem_to_FirstGlobTask = struct;
%[Task,TaskDepGraph] = GenerateST_FI_Tasks_4D_ST;
[Task,TaskDepGraph,SpElemProperties,Map_SpElem_to_FirstGlobTask] ...
                            = GenerateSp_FI_Tasks_4D_ST(sC,sD,SpElemProperties,Num_of_Elem,Task,TaskDepGraph,Map_SpElem_to_FirstGlobTask);
TaskOrder                   = SortTasks(TaskDepGraph,Map_SpElem_to_FirstGlobTask,sC,SpElemProperties);
%clearvars TaskDepGraph;
disp('point3')
disp('point4')
[kappa,FaceArea]            = ComputeKappa_4D_ST(cdt,sG,sC,sD,D0,D1,D2,D3,NodePos,SpElemProperties,Num_of_Elem);
disp('point4.1')
Z                           = ComputeImpedance_for_EachSTPs(RefImpedance_SpV,sC,sD,SpElemProperties,Num_of_Elem);
Kappa_over_Z                = kappa./Z;
disp('point5')

%% 
Source                      = SourceFDTD(MeshMeasurements,SpElemProperties,sC);
TMM                         = ConstructTimeMarchingMatrix_4D_ST(D1,D2,sC,Kappa_over_Z,Source,SpElemProperties,Task,TaskOrder,Num_of_Elem);
[TMM_Fields, TMM_Sources]   = SplitTMM_into_FieldsAndSources(TMM,Source,Num_of_Elem);

%%
%DoFs_FacesThenEdges        = InitializeFields;
DoFs_FacesThenEdges         = zeros(Num_of_Elem.SpP+Num_of_Elem.SpS,1);
%DoFs_FacesThenEdges         = rand(Num_of_Elem.SpP+Num_of_Elem.SpS,1);
Num_of_Steps                = 50;
disp(['Number of Steps = ', num2str(Num_of_Steps)])
Time                        = 0;
disp('point6')
DoFs_FacesThenEdges         = TimeMarch(Num_of_Steps,Time,cdt,TMM_Fields,TMM_Sources,DoFs_FacesThenEdges,Source);
ZConst                      = 5;
PlotMagneticFluxDensity_atZEquals(ZConst,DoFs_FacesThenEdges,FaceArea,PrimFacePos,SpElemProperties,MeshMeasurements)
YConst                      = 5;
PlotMagneticFluxDensity2D_at_YEquals(DoFs_FacesThenEdges,YConst,FaceArea,Num_of_Elem,PrimFacePos,MeshMeasurements)
PlotMagneticFluxDensity3D(DoFs_FacesThenEdges,FaceArea,Num_of_Elem,PrimFacePos,MeshMeasurements)