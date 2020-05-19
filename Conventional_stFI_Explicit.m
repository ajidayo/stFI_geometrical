clear;

global DIM
global DISPDEBUGGINGMESSAGE
global DISPCAUTIONMESSAGE
global DOEIGVANALYSIS

DIM=2;
DISPDEBUGGINGMESSAGE=true;
DISPCAUTIONMESSAGE=true;
DOEIGVANALYSIS=true;
number_of_steps=100;
cdt=0.30;

%% test
MeshParam.Size_X=10;
MeshParam.Size_Y=10;
MeshParam.Fine_X_from=6;
MeshParam.Fine_X_to=8;
MeshParam.Fine_Y_from=6;
MeshParam.Fine_Y_to=8;
MeshParam.deltaboundary=1.0/12.0;
%% end test

F=struct('Ex',[],'Ey',[],'Bz',[]);
for j=1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        if MeshParam.Fine_X_from <=i && i <= MeshParam.Fine_X_to ...
                && MeshParam.Fine_Y_from <= j && j <= MeshParam.Fine_Y_to
            F(i,j).Ex=zeros(2,2);
            F(i,j).Ey=zeros(2,2);
            F(i,j).Bz=zeros(2,2);
        elseif i == MeshParam.Fine_X_to+1 && MeshParam.Fine_Y_from <= j && j <= MeshParam.Fine_Y_to
            F(i,j).Ex=zeros(1,1);
            F(i,j).Ey=zeros(1,2);
            F(i,j).Bz=zeros(1,1);        
        elseif j == MeshParam.Fine_Y_to+1 && MeshParam.Fine_X_from <= i && i <= MeshParam.Fine_X_to 
            F(i,j).Ex=zeros(2,1);
            F(i,j).Ey=zeros(1,1);
            F(i,j).Bz=zeros(1,1);
        else
            F(i,j).Ex=zeros(1,1);
            F(i,j).Ey=zeros(1,1);
            F(i,j).Bz=zeros(1,1);
        end
    end
end

%%
GaussParam.Ampl=1;
GaussParam.relaxfact=1;
GaussParam.Xcenter=0.5*MeshMeasurements.XCoord;
GaussParam.Ycenter=0.5*MeshMeasurements.YCoord;
disp(['A sigma =',num2str(GaussParam.Ampl),', ', num2str(GaussParam.relaxfact)])
disp(['gauss_center.x,  gauss_center.y =',num2str(gauss_center.x), ', ', num2str(gauss_center.y)])

%%
MeshParam.UpdateNum_subgrid=2;
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

%%
for timestep=1:number_of_steps
   F=Update_Conventional_stFI_Explicit(F,cdt,MeshParam,Area);
end
