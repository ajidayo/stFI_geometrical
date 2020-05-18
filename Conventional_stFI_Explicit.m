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

%% test
MeshParam.Size_X=10;
MeshParam.Size_Y=10;
MeshParam.Fine_X_from=6;
MeshParam.Fine_X_to=8;
MeshParam.Fine_Y_from=6;
MeshParam.Fine_Y_to=8;
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
            F(i,j).Ex=zeros(1,2);
            F(i,j).Ey=zeros(1,1);
            F(i,j).Bz=zeros(1,1);
        end
    end
end

for timestep=1:number_of_steps
   F=Update_Conventional_stFI_Explicit(F,MeshParam);
end
