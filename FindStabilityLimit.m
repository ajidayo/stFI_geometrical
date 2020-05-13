clear;

global EPSILON
global DIM

EPSILON = 10^(-7);
DIM = 2; %number of spatial dimensions

%% Spatial Meshing and Impedance for belt-like subgrid region only with square faces

Img_MeshMeasLocation=imread('MeshMeasurements_square_belt.png');
image(Img_MeshMeasLocation)
% MeshMeasurements=MeshMeasurements_100times100_SquareBelt_belt65to85;

%% test

MeshMeasurements.XCoord=10;
MeshMeasurements.YCoord=10;
MeshMeasurements.FineStartAtYCoord=6;
MeshMeasurements.FineEndAtYCoord=8;

%%
MeshParam = MeshParameters_square_belt(MeshMeasurements);
MeshParam.deltaboundary=1.0/12.0;
UpdateNum_belt=2;
[sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
    = GenerateMesh_square_belt(UpdateNum_belt,MeshParam);

Zinv_p=ones(MeshNum.P,1);
%% Calculate Constitutive Equation

% Future tasks; modify att into nested structures like att.e(e).bound
att = attribute_f_and_e(sC,sG,UpdateNum, MeshNum);


%% Divide_into_induced_subgraphs

% #4: combine both space-time and spatial FI in order to make the size of D small
[subG_bin,subG_sizes,allIdx_stFI,UpdateNum] ...
    = Divide_into_induced_subgraphs(sC,UpdateNum,MeshNum,att);
[Taskorder,task,D,Ctrans] ...
    = Obtain_TaskOrderandIncMat(sC,UpdateNum,allIdx_stFI,subG_bin,subG_sizes,att,first_pIdx,MeshNum);

%% find stability limit of cdt

FuncDomain.from=0;
FuncDomain.to=1;
SearchEpsilon=10^(-4);
StabLim_cdt...
    =BinarySearch_forLogicalFunc...
    (@isstable,FuncDomain,SearchEpsilon,Taskorder,task,Zinv_p,D,Ctrans,sC,sG,UpdateNum,edgevec,subG_bin,first_pIdx,att,MeshNum)
% task: put argumnts for function isstable into a cell to pass them more easily 