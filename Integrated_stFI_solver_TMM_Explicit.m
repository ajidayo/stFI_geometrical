clear;
% Numbers like PNum, FNum, ... should be put together as an structure.
% e.g. MeshNums.FNum=hoge, HodgeGridNum.PNum=hoge, HodgeGridNum.OmegaNum=hoge.
% then it is easy to pass on as an arguement to each functions

% future tasks; Modify Parameters_Mesh into function with arguements such
% as (Size_X,Size_Y,...)

Parameters_Mesh
[sC,sG,denominator,edgevec,first_p,tilde_node_position,MeshNum,MeshParam] ...
    = GenerateMesh_triangular(denominator_obi,MeshParam);

%first_p_for=struct('f',first_p.f,'e',first_p.e);

cdt=0.41

% Future tasks; adapt Constitutive to partially non-orthogonal grids
[kappa,b_area,MeshNum] ...
    =Constitutive(cdt,sC,sG,denominator,edgevec,first_p,tilde_node_position,MeshNum);

first_i_triangle_scatterer=20;
ImpedanceParam.freespace=1.0;
ImpedanceParam.medium=0.01;
[impedance_inv_p] ...
    = impedance_triangular(ImpedanceParam,first_i_triangle_scatterer,sC,denominator,first_p,MeshNum,MeshParam);

A=1;
sigma=10;

InitVal ...
    =Gaussian_DeadCenter_triangle(A,sigma,tilde_node_position,b_area,MeshNum,MeshParam);

kappatimesz=kappa.*impedance_inv_p;
%Zinverse=spdiags(kappatimesz,0,MeshNum.P,MeshNum.P);

%% calculate explicitly using Time-marching Matrix

% #4: combine both space-time and spatial FI in order to make the size of D small
[subG_bin,subG_sizes,allIdx_stFI,denominator,att] ...
    = Divide_into_induced_subgraphs(sC,sG,denominator,MeshNum);

[TMM_Explicit] ...
    = Obtain_TMM_Explicit(kappatimesz,sC,denominator,allIdx_stFI,subG_bin,subG_sizes,att,first_p,MeshNum);

number_of_steps=100

[variables_f_then_e] ...
    =execution_with_TMM_Explicit(TMM_Explicit,InitVal,number_of_steps);

b_f=variables_f_then_e(1:MeshNum.F);

plot_bface_general(b_f,b_area,tilde_node_position,MeshParam,MeshNum)

%%

%b_error=b_reduced-b_obi;
%plot_bface_general(b_error,b_area,tilde_node_position,Size_X,Size_Y,FNum)


%% far-future tasks
% include PML, PEC tasks.