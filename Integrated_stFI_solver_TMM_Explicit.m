clear;

global EPSILON
global DIM
global DISPDEBUGGINGMESSAGE
global DISPCAUTIONMESSAGE
global DOEIGVANALYSIS

EPSILON = 10^(-7)
DIM = 2; %number of spatial dimensions: DO NOT CHANGE unless you rewrited the codes for 3D analysis
DISPCAUTIONMESSAGE =true
DISPDEBUGGINGMESSAGE = true
DOEIGVANALYSIS = true

%% Input (1): Spatial Meshing and Impedance for belt-like subgrid region with triangle faces 

% Img_MeshMeasLocation=imread('MeshMeasurements_triangle_belt.png');
% image(Img_MeshMeasLocation)

% % task; unify two function MeshParameters_hoge... and GenerateMesh_hoge
% MeshMeasurements=MeshMeasurements_100times100_withTriangle;
% MeshParam = MeshParameters_triangle_belt(MeshMeasurements);
% MeshParam.deltatriangle=0.1;
% MeshParam.deltaboundary=1.0/12.0;
% UpdateNum_belt=2;
% [sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
%     = GenerateMesh_triangular_belt(UpdateNum_belt,MeshParam);

% first_i_triangle_scatterer=20;
% ImpedanceParam.freespace=1.0;
% ImpedanceParam.medium=0.01;
% Zinv_p ...
%     = Impedance_TriangleScatterer(ImpedanceParam,first_i_triangle_scatterer,sC,UpdateNum,first_pIdx,MeshNum,MeshParam);

% disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')
% gauss_center.x=MeshParam.Size_X/2.0;
% gauss_center.y=0.5*(MeshParam.Fine_Y_from-1 ...
%     +(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0...
%     +(MeshParam.Size_Y-MeshParam.Fine_Y_to));

%% Input (2): Spatial Meshing and Impedance for belt-like subgrid region only with square faces

% Img_MeshMeasLocation=imread('MeshMeasurements_square_belt.png');
% image(Img_MeshMeasLocation)
% % MeshMeasurements=MeshMeasurements_100times100_SquareBelt_belt65to85;
% 
% % test
% MeshMeasurements.XCoord=10;
% MeshMeasurements.YCoord=10;
% MeshMeasurements.FineStartAtYCoord=6;
% MeshMeasurements.FineEndAtYCoord=8;
% ScattererMeasurements.FromXCoord=6.5;
% ScattererMeasurements.ToXCoord=7.5;
% ScattererMeasurements.FromYCoord=6.5;
% ScattererMeasurements.ToYCoord=7.5;
% GaussParam.Ampl=1;
% GaussParam.relaxfact=1;
% % end test

% % task; unify two function MeshParameters_hoge... and GenerateMesh_hoge
% MeshParam = MeshParameters_square_belt(MeshMeasurements);
% MeshParam.deltaboundary=1.0/12.0;
% UpdateNum_belt=2;
% [sC,sG,UpdateNum,edgevec,first_pIdx,tilde_node_position,MeshNum,MeshParam] ...
%     = GenerateMesh_square_belt(UpdateNum_belt,MeshParam);
% 
% % ScattererMeasurements.FromXCoord=17;
% % ScattererMeasurements.ToXCoord=33;
% % ScattererMeasurements.FromYCoord=67;
% % ScattererMeasurements.ToYCoord=83;
% ImpedanceParam.freespace=1.0;
% ImpedanceParam.medium=1.0;
% %ImpedanceParam.medium=0.01;
% %Zinv_p=ones(MeshNum.P,1);
%  Zinv_p...
%      = Impedance_SquareScatterer_BeltlikeSubgrid(ImpedanceParam,ScattererMeasurements,sC,UpdateNum,first_pIdx,MeshNum,MeshParam,MeshMeasurements);
% disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')
% gauss_center.x=MeshParam.Size_X/2.0;
% gauss_center.y=0.5*(MeshParam.Fine_Y_from-1 ...
%     +(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0...
%     +(MeshParam.Size_Y-MeshParam.Fine_Y_to));

%% Input (3): Spatial Meshing and Impedance for square-like subgrid region only with square faces

Img_MeshMeasLocation=imread('MeshMeasurements_squarefaces_squaresubgrid.png');
image(Img_MeshMeasLocation)

% test 1
% MeshMeasurements.XCoord=100;
% MeshMeasurements.YCoord=100;
% MeshMeasurements.FineStartAtXCoord=15;
% MeshMeasurements.FineEndAtXCoord=35;
% MeshMeasurements.FineStartAtYCoord=15;
% MeshMeasurements.FineEndAtYCoord=35;
% ScattererMeasurements.FromXCoord=20;
% ScattererMeasurements.ToXCoord=30;
% ScattererMeasurements.FromYCoord=20;
% ScattererMeasurements.ToYCoord=30;
% test 1 till here

% test 2
MeshMeasurements.XCoord=10;
MeshMeasurements.YCoord=10;
MeshMeasurements.FineStartAtXCoord=6;
MeshMeasurements.FineEndAtXCoord=9;
MeshMeasurements.FineStartAtYCoord=6;
MeshMeasurements.FineEndAtYCoord=9;
ScattererMeasurements.FromXCoord=7;
ScattererMeasurements.ToXCoord=8;
ScattererMeasurements.FromYCoord=7;
ScattererMeasurements.ToYCoord=8;
% test 2 till here

% task; unify two function MeshParameters_hoge... and GenerateMesh_hoge
MeshParam = MeshParameters_squarefaces_squaresubgrid(MeshMeasurements);
MeshParam.deltaboundary=1.0/12.0;
MeshParam.deltacorner=0; % Not used yet
UpdateNum_subgrid=2;
[sC,sG,UpdateNum,edgevec,first_pIdx,tilde_f,MeshNum,MeshParam] ...
    = GenerateMesh_squarefaces_squaresubgrid(UpdateNum_subgrid,MeshParam);

ImpedanceParam.freespace=1.0;
ImpedanceParam.medium=1.0;
%Zinv_p=ones(MeshNum.P,1);
 Zinv_p...
     = Impedance_SquareScatterer(ImpedanceParam,ScattererMeasurements,sC,UpdateNum,first_pIdx,MeshNum,MeshParam,MeshMeasurements);
disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')

GaussParam.Ampl=1;
GaussParam.relaxfact=1;
gauss_center.x=0.5*MeshMeasurements.XCoord;
gauss_center.y=0.5*MeshMeasurements.YCoord;

%% Calculate Discrete Hodge operator (Constitutive equation, [z]^-1)

% Future tasks; modify att into nested structures like att.e(e).bound
att = attribute_f_and_e(sC,sG,UpdateNum, MeshNum);

cdt=0.55

% Future tasks; utilize spatial-FI-like calculation in Constitutive
[kappa,Area_spatialfaces,att,MeshNum]=Constitutive(cdt,sC,sG,UpdateNum,edgevec,first_pIdx,att,MeshNum);
kappaoverZ=kappa.*Zinv_p;

%% calculating initial distribution of Bz

InitVal ...
    =GaussianDistributBz(GaussParam,tilde_f,Area_spatialfaces,MeshNum,gauss_center);

%% Obtain Time-marching Matrix  for single timestep

% #4: combine both space-time and spatial FI in order to make the size of D small
[subG_bin,subG_sizes,allIdx_stFI,UpdateNum] ...
    = Divide_into_induced_subgraphs(sC,UpdateNum,MeshNum,att);

[Taskorder,task,D,Ctrans] ...
    = Obtain_TaskOrderandIncMat(sC,UpdateNum,allIdx_stFI,subG_bin,subG_sizes,att,first_pIdx,MeshNum);

D_tildeD_Zinv=[D;Ctrans * spdiags(kappaoverZ,0,MeshNum.P,MeshNum.P)];
clearvars D Ctrans

[TMM_Explicit] ...
    = Construct_TMM_Explicit(Taskorder,task,D_tildeD_Zinv,kappaoverZ,sC,UpdateNum,subG_bin,first_pIdx,MeshNum);

%% Execute Explicit Calculation 

time=0;
variables_f_then_e=[InitVal.f; InitVal.e];

number_of_steps=1000000

CalPeriod=cdt * number_of_steps;
disp(['Executing Calculation: from ct = ',num2str(time), ' to ct = ',num2str(time+CalPeriod),' with cdt = ',num2str(cdt)])

time = time + CalPeriod;
[variables_f_then_e] ...
    = execution_with_TMM_Explicit(TMM_Explicit,variables_f_then_e,number_of_steps);

b_f = variables_f_then_e(1:MeshNum.F);
disp(['plotting Bz at ct = ', num2str(time)])
plot_bface_general(b_f,Area_spatialfaces,tilde_f,MeshParam,MeshNum)

%% Eigenvalue Analysis
if DOEIGVANALYSIS ==true
    disp('Calculating Eigenvalues')
    
    eigenvalues = eigs(TMM_Explicit,20,'largestabs');
    figure
    theta = linspace(0,2*pi);
    x = cos(theta);
    y = sin(theta);
    eigv_re=real(eigenvalues);
    eigv_im=imag(eigenvalues);
    plot(eigv_re,eigv_im,'or',x,y,'-b')
    axis equal
    EigValAbs=abs(eigenvalues);
    EigvEpsilon=10^(-12);
    % NoEigvConvergeFlag = ~any(~isnan(EigValAbs));
    % AnyEigvNonConvergeFlag = any(isnan(EigValAbs));
    % NaN entry is blocked by "find" since comparison including NaN always
    % returns false.
    IdxUnstabEigVal=find(EigValAbs>1+EigvEpsilon);
    if  ~any(~isnan(EigValAbs))
        disp('None of the eigenvalues congverged.')
    elseif any(isnan(EigValAbs))
        disp('Some of the eigenvalues did not congverge')
        if size(IdxUnstabEigVal,1)~=0
            disp(['unstable for cdt = ',num2str(cdt),' (EigvEpsilon = ',num2str(EigvEpsilon),')'])
            %disp(EigValAbs)
        else
            disp('(The absolute values of the converged eigenvalues were equal to or less than unity.)')
        end
    else
        if size(IdxUnstabEigVal,1)==0
            disp(['stable for cdt = ',num2str(cdt),' (EigvEpsilon = ',num2str(EigvEpsilon),')'])
            %disp(EigValAbs)
        else
            disp(['unstable for cdt = ',num2str(cdt),' (EigvEpsilon = ',num2str(EigvEpsilon),')'])
            %disp(EigValAbs)
        end
    end
end

%% Error to Conventional stFI

%b_error=b_reduced-b_obi;
%plot_bface_general(b_error,b_area,tilde_node_position,Size_X,Size_Y,FNum)


%% far-future tasks
% include PML, PEC tasks.