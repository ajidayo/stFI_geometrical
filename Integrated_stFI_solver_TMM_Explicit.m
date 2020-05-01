clear;

global EPSILON
global DIM

EPSILON=10^(-7);
DIM=2; %number of spatial dimensions


%% Spatial Meshing and Impedance for belt-like subgrid region with triangle faces 
Location_MeshMeas=imread('MeshMeasurements.png');
image(Location_MeshMeas)
MeshMeasurements=MeshMeasurements_100times100_withTriangle;
MeshParam = MeshParameters_triangle_belt(MeshMeasurements);
MeshParam.deltatriangle=0.1;
MeshParam.deltaboundary=1.0/12.0;
UpdateNum_obi=2;
[sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
    = GenerateMesh_triangular(UpdateNum_obi,MeshParam);

first_i_triangle_scatterer=20;
ImpedanceParam.freespace=1.0;
ImpedanceParam.medium=0.01;
[impedance_inv_p] ...
    = impedance_triangular(ImpedanceParam,first_i_triangle_scatterer,sC,UpdateNum,first_pIdx,MeshNum,MeshParam);

%%

% Future tasks; modify att into nested structures like att.e(e).bound
att = attribute_f_and_e(sC,sG,UpdateNum, MeshNum);

cdt=0.41;

% Future tasks; adapt Constitutive to partially non-orthogonal grids:DONE
% but not been tested yet
% Future tasks; adapt Constitutive to subgrid corners
% Future tasks; utilize spatial-FI-like calculation in Constitutive
[kappa,b_area,att,MeshNum] = Constitutive(cdt,sC,sG,UpdateNum,edgevec,first_pIdx,att,MeshNum);

%% calculating initial distribution

GaussParam.Ampl=1;
GaussParam.relaxfact=10;

InitVal ...
    =Gaussian_DeadCenter(GaussParam,tilde_node_position,b_area,MeshNum,MeshParam);

kappaoverZ=kappa.*impedance_inv_p;

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