% clear;

global DIM
global DISPDEBUGGINGMESSAGE
global DISPCAUTIONMESSAGE
global DOEIGVANALYSIS

DIM=2;
DISPDEBUGGINGMESSAGE=true;
DISPCAUTIONMESSAGE=true;
DOEIGVANALYSIS=true;
% number_of_steps=1;
% cdt=0.5;

RANDOMLYINITIALIZE=false;

%% test
% MeshParam.Size_X=100;
% MeshParam.Size_Y=100;
% MeshParam.Fine_X_from=21;
% MeshParam.Fine_X_to=40;
% MeshParam.Fine_Y_from=21;
% MeshParam.Fine_Y_to=40;
% MeshParam.deltaboundary=1.0/12.0;
%% end test

Updated=zeros(MeshParam.Size_X,MeshParam.Size_Y);

F=struct('Ex',[],'Ey',[],'Bz',[],'ExPreserve',[],'EyPreserve',[]);
for j=1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        if MeshParam.Fine_X_from <=i && i <= MeshParam.Fine_X_to ...
                && MeshParam.Fine_Y_from <= j && j <= MeshParam.Fine_Y_to
            F(i,j).Ex=zeros(2,2);
            F(i,j).Ey=zeros(2,2);
            F(i,j).Bz=zeros(2,2);
            if i== MeshParam.Fine_X_from && j== MeshParam.Fine_Y_from
                F(i,j).EyPreserve=0;
                F(i,j).ExPreserve=0;
            elseif i== MeshParam.Fine_X_from
                F(i,j).EyPreserve=0;
            elseif j== MeshParam.Fine_Y_from
                F(i,j).ExPreserve=0;
            end
            if RANDOMLYINITIALIZE==true
                F(i,j).Ex=0.5*rand(2,2);
                F(i,j).Ey=0.5*rand(2,2);
                F(i,j).Bz=0.25*rand(2,2);
            end
        elseif i == MeshParam.Fine_X_to+1 && MeshParam.Fine_Y_from <= j && j <= MeshParam.Fine_Y_to
            F(i,j).Ex=zeros(1,1);
            F(i,j).Ey=zeros(1,2);
            F(i,j).EyPreserve=0;
            F(i,j).Bz=zeros(1,1);
            if RANDOMLYINITIALIZE==true
                F(i,j).Ex=rand(1,1);
                F(i,j).Ey=0.5*rand(1,2);
                F(i,j).Bz=rand(1,1);
            end
        elseif j == MeshParam.Fine_Y_to+1 && MeshParam.Fine_X_from <= i && i <= MeshParam.Fine_X_to
            F(i,j).Ex=zeros(2,1);
            F(i,j).ExPreserve=0;
            F(i,j).Ey=zeros(1,1);
            F(i,j).Bz=zeros(1,1);
            if RANDOMLYINITIALIZE==true
                F(i,j).Ex=0.5*rand(2,1);
                F(i,j).Ey=rand(1,1);
                F(i,j).Bz=rand(1,1);
            end
        else
            F(i,j).Ex=zeros(1,1);
            F(i,j).Ey=zeros(1,1);
            F(i,j).Bz=zeros(1,1);
            if RANDOMLYINITIALIZE==true
                F(i,j).Ex=rand(1,1);
                F(i,j).Ey=rand(1,1);
                F(i,j).Bz=rand(1,1);
            end
        end
    end
end

%%
GaussParam.Ampl=1;
GaussParam.relaxfact=10;
GaussParam.XCenter=0.5*MeshParam.Size_X;
GaussParam.YCenter=0.5*MeshParam.Size_Y;
disp(['A sigma =',num2str(GaussParam.Ampl),', ', num2str(GaussParam.relaxfact)])
disp(['GaussParam.XCenter,  GaussParam.YCenter =',num2str(GaussParam.XCenter), ', ', num2str(GaussParam.YCenter)])

%%
MeshParam.UpdateNum_subgrid=2;
cdt_subgrid=cdt/MeshParam.UpdateNum_subgrid;
delta=MeshParam.deltaboundary;
Area.Bz_next2outercorner=0.5*(1-delta+1-delta+1.0/6.0)/2.0 ...
    +(0.5+delta)*(1-delta+1.0/6.0+1-delta)/2.0 ...
    -0.5*delta*(1-delta);
Area.Bz_outercorner=2*0.5*1.0*(1-delta);
Area.Bz_outsideboundary=2*0.5*(1-delta+1.0/6.0+1-delta)/2.0;
Area.E_outsideboundary=(1.0-MeshParam.deltaboundary)*cdt-cdt^3/12.0;
Area.Bz_insideboundary=0.5*(0.5+delta+0.5+delta-1.0/6.0)*0.5;
Area.Bz_innercorner=2*0.5*(0.5+delta-1.0/6.0)*(0.5+delta);
Area.E_insideboundary_smaller=0.5*(0.5+delta-1.0/6.0+0.5+delta-1.0/6.0+cdt^2/6.0)*cdt_subgrid;
Area.E_insideboundary_larger =0.5*(0.5+delta        +0.5+delta        +cdt^2/6.0)*cdt_subgrid;
Area.E_innercorner=0.5*( ...
    ((sqrt(10))^(-1))*(1.5+3*MeshParam.deltaboundary+1.0/6.0)...
    +((sqrt(10))^(-1))*(1.5+3*MeshParam.deltaboundary+cdt^2/2.0+1.0/6.0)...
    )*cdt_subgrid;

if ~RANDOMLYINITIALIZE
    [F] = Initialize_Conventional_stFI_Explicit_w_GaussianDistribution(F,GaussParam,MeshParam,Area);
end


%%
for timestep=1:number_of_steps
    timestep
   [F,Updated] = Update_Conventional_stFI_Explicit(F,cdt,MeshParam,Area,Updated);
end

[B_mesh_Conventional] = Plot_Bz_ConventionalstFI(F,MeshParam,Area);

