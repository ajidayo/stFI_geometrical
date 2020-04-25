clear;
globals
% Numbers like PNum, FNum, ... should be put together as an structure.
% e.g. MeshNums.FNum=hoge, HodgeGridNum.PNum=hoge, HodgeGridNum.OmegaNum=hoge.
% then it is easy to pass on as an arguement to each functions
Parameters_Mesh
GenerateMesh_triangular

% Future tasks; translate GenerateD and Constitutive into functions 
% Future tasks; store denominators as an structure
%GenerateD_sG_sC_list_3

cdt=0.1;

% Future tasks; adapt Constitutive to partially non-orthogonal grids
Constitutive_directionvector_revised % calculates kapp

Parameters_initial_conditions
initialize_variable

Parameters_scatterer
impedance_triangular
kappatimesz=kappa.*impedance_inverse_p;
Zinverse=spdiags(kappatimesz,0,PNum,PNum);
%D_tildeD_Zinverse=[D;Ctransposed*Zinverse].';


% Future tasks; translate Order into a function 
%Order_undummied_st

number_of_steps=100
% Is it possible to visualize the results in each step in real time? 
%execution_sorted
%plot_bface_general(b_obi,b_area,tilde_node_position,Size_X,Size_Y,FNum)

%% calculate explicitly using Time-marching Matrix

% #1: convert D into a digraph object and then bypass the (intermidiate?) nodes  
% Obtain_TimeMarchingMatrix_Explicit_graph
% #2: delete (intermidiate?) variables directly from D 
% Obtain_TimeMarchingMatrix_Explicit;
% #3: delete (intermidiate?) variables directly from D (thought it fastens up but it didn't)
%Obtain_TimeMarchingMatrix_Explicit_planB;
% #4: combine both space-time and spatial FI in order to make the size of D small
Divide_into_induced_subgraphs;
Obtain_TimeMarchingMatrix_Explicit_subgraphs;

number_of_steps=100
execution_with_TimeMarchingMatrix_Explicit;

plot_bface_general(b_reduced,b_area,tilde_node_position,Size_X,Size_Y,FNum)

%%

b_error=b_reduced-b_obi;
plot_bface_general(b_error,b_area,tilde_node_position,Size_X,Size_Y,FNum)