clear;

global EPSILON
global DIM

EPSILON=10^(-7);
DIM=2; %number of spatial dimensions


%% Spatial Meshing and Impedance for belt-like subgrid region with triangle faces 
Img_MeshMeasLocation=imread('MeshMeasurements_triangle_belt.png');
image(Img_MeshMeasLocation)
MeshMeasurements=MeshMeasurements_100times100_withTriangle;
MeshParam = MeshParameters_triangle_belt(MeshMeasurements);
MeshParam.deltatriangle=0.1;
MeshParam.deltaboundary=1.0/12.0;
UpdateNum_belt=2;
[sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
    = GenerateMesh_triangular_belt(UpdateNum_belt,MeshParam);

first_i_triangle_scatterer=20;
ImpedanceParam.freespace=1.0;
ImpedanceParam.medium=0.01;
[Zinv_p] ...
    = impedance_triangular(ImpedanceParam,first_i_triangle_scatterer,sC,UpdateNum,first_pIdx,MeshNum,MeshParam);

disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')
gauss_center.x=MeshParam.Size_X/2.0;
gauss_center.y=0.5*(MeshParam.Fine_Y_from-1 ...
    +(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0...
    +(MeshParam.Size_Y-MeshParam.Fine_Y_to));

%% Spatial Meshing and Impedance for belt-like subgrid region only with square faces

Img_MeshMeasLocation=imread('MeshMeasurements_square_belt.png');
image(Img_MeshMeasLocation)
MeshMeasurements=MeshMeasurements_100times100_SquareBelt;
MeshParam = MeshParameters_square_belt(MeshMeasurements);
MeshParam.deltaboundary=1.0/12.0;
UpdateNum_belt=2;
[sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
    = GenerateMesh_square_belt(UpdateNum_belt,MeshParam);

ImpedanceParam.freespace=1.0;
ImpedanceParam.medium=0.01;
Zinv_p = ImpedanceParam.freespace*eyes(MeshNum.P);

disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')
gauss_center.x=MeshParam.Size_X/2.0;
gauss_center.y=0.5*(MeshParam.Fine_Y_from-1 ...
    +(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0...
    +(MeshParam.Size_Y-MeshParam.Fine_Y_to));

%% Calculate Constitutive Equation

% Future tasks; modify att into nested structures like att.e(e).bound
att = attribute_f_and_e(sC,sG,UpdateNum, MeshNum);

cdt=0.41;

% Future tasks; adapt Constitutive to partially non-orthogonal grids:DONE
% but not been tested yet
% Future tasks; adapt Constitutive to subgrid corners
% Future tasks; utilize spatial-FI-like calculation in Constitutive
[kappa,b_area,att,MeshNum] = Constitutive(cdt,sC,sG,UpdateNum,edgevec,first_pIdx,att,MeshNum);
kappaoverZ=kappa.*Zinv_p;
%% calculating initial distribution

GaussParam.Ampl=1;
GaussParam.relaxfact=10;

InitVal ...
    =GaussianDistributBz(GaussParam,tilde_node_position,b_area,MeshNum,gauss_center);

%% Obtain Time-marching Matrix

% #4: combine both space-time and spatial FI in order to make the size of D small
[subG_bin,subG_sizes,allIdx_stFI,UpdateNum] ...
    = Divide_into_induced_subgraphs(sC,UpdateNum,MeshNum,att);

% task: allocate TMM beforehand to reduce overheads
[TMM_Explicit] ...
    = Obtain_TMM_Explicit(kappaoverZ,sC,UpdateNum,allIdx_stFI,subG_bin,subG_sizes,att,first_pIdx,MeshNum);

%% Execute Explicit Calculation

time=0;
variables_f_then_e=[InitVal.f; InitVal.e];

number_of_steps=100

CalPeriod=cdt * number_of_steps;
disp(['Executing Calculation: from ct = ',num2str(time), ' to ct = ',num2str(time+CalPeriod)])

time = time + CalPeriod;
[variables_f_then_e] ...
    =execution_with_TMM_Explicit(TMM_Explicit,variables_f_then_e,number_of_steps);

b_f=variables_f_then_e(1:MeshNum.F);
disp(['plotting Bz at ct =' num2str(time)])
plot_bface_general(b_f,b_area,tilde_node_position,MeshParam,MeshNum)

 
%% Eigenvalue Analysis

% eigenvalues = eigs(TMM_Explicit,20,'largestabs','Tolerance',1e-3);
% plot(eigenvalues)
% IdxUnstabEigVal=find(abs(eigenvalues)>1+EPSILON)

%% Error to Conventional stFI

%b_error=b_reduced-b_obi;
%plot_bface_general(b_error,b_area,tilde_node_position,Size_X,Size_Y,FNum)


%% far-future tasks
% include PML, PEC tasks.