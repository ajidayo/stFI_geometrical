function [sC,sG,denominator,edgevec,first_p,tilde_node_position,MeshNum,MeshParam]...
    = GenerateMesh_triangular(denominator_obi,MeshParam)

disp('GenerateMesh_triangular:CALLED')

global DIM

%global MeshNum

disp('Initializing Spatial Mesh Information ')

MeshParam.coarse_FNum_former=MeshParam.Size_X*(MeshParam.Fine_Y_from-1);
MeshParam.fine_FNum_former=MeshParam.Size_X*2*(MeshParam.Triangle_Y_from-MeshParam.Fine_Y_from);
MeshParam.triangle_FNum=MeshParam.Size_X*4*(MeshParam.Triangle_Y_to-MeshParam.Triangle_Y_from+1);
MeshParam.fine_FNum_latter=MeshParam.Size_X*2*(MeshParam.Fine_Y_to-MeshParam.Triangle_Y_to);
MeshParam.coarse_FNum_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to);
FNum=MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter+MeshParam.coarse_FNum_latter;
MeshNum.F=FNum;

MeshParam.coarse_ENum_X_former=MeshParam.coarse_FNum_former;
MeshParam.fine_ENum_X_former=MeshParam.fine_FNum_former;
MeshParam.triangle_ENum_XandDiag=MeshParam.triangle_FNum;
MeshParam.fine_ENum_X_latter=MeshParam.Size_X*2*(MeshParam.Fine_Y_to-MeshParam.Triangle_Y_to+1);
MeshParam.coarse_ENum_X_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to-1);
MeshParam.ENum_XandDiag=MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter+MeshParam.coarse_ENum_X_latter;

MeshParam.coarse_ENum_Y_former=MeshParam.coarse_FNum_former;
MeshParam.fine_ENum_Y_former=MeshParam.fine_FNum_former;
MeshParam.triangle_ENum_Y=MeshParam.Size_X*2*(MeshParam.Triangle_Y_to-MeshParam.Triangle_Y_from+1);
MeshParam.fine_ENum_Y_latter=MeshParam.fine_FNum_latter;
MeshParam.coarse_ENum_Y_latter=MeshParam.coarse_FNum_latter;
MeshParam.ENum_Y=MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.coarse_ENum_Y_latter;
ENum=MeshParam.ENum_XandDiag+MeshParam.ENum_Y;
MeshNum.E=ENum;


MeshParam.coarse_NNum_former=MeshParam.Size_X*(MeshParam.Fine_Y_from-1);
MeshParam.fine_NNum_former=MeshParam.Size_X*2*(MeshParam.Triangle_Y_from-MeshParam.Fine_Y_from);
MeshParam.triangle_NNum=MeshParam.Size_X*2*(MeshParam.Triangle_Y_to-MeshParam.Triangle_Y_from+1);
MeshParam.fine_NNum_latter=MeshParam.Size_X*2*(MeshParam.Fine_Y_to-MeshParam.Triangle_Y_to+1);
MeshParam.coarse_NNum_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to-1);
NNum=MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter+MeshParam.coarse_NNum_latter;
MeshNum.N=NNum;

%% check
disp('Check; [FNum, ENum, NNum]=')
disp([MeshNum.F MeshNum.E MeshNum.N])


%% Allocate Arrays
sC=sparse(FNum,ENum);
sG=sparse(ENum,NNum);
denominator_node=zeros(NNum,1);
denominator_edge=zeros(ENum,1);
denominator_face=zeros(FNum,1);
edge_vector=zeros(ENum,DIM);
tilde_edge_vector=zeros(ENum,DIM);
first_p_for_f=zeros(FNum,1);
first_p_for_e=zeros(ENum,1);
%first_Omega_for_f=zeros(FNum,1);
%first_Omega_for_e=zeros(ENum,1);

tilde_node_position=zeros(FNum,DIM);


%% initialize sC

% coarce-quadrilateral zone 1: main
for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X-1
        f=MeshParam.Size_X*(j-1)+i;
        e=f;
        sC(f,e)=1;
        e=f+MeshParam.Size_X;
        sC(f,e)=-1;
        e=f+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=f+1+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
%        disp([f e])
%        pause(0.05)
    end
    i=MeshParam.Size_X;
    f=MeshParam.Size_X*(j-1)+i;
    e=f;
    sC(f,e)=1;
    e=f+MeshParam.Size_X;
    sC(f,e)=-1;
    e=f+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=f-(MeshParam.Size_X-1)+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% coarce-quadrilateral zone 1: the row right outside the boundary
%j=MeshParam.Fine_Y_from-1;
for i=1:MeshParam.Size_X-1
    f=i+MeshParam.coarse_FNum_former-MeshParam.Size_X;
    e=i+MeshParam.coarse_ENum_X_former-MeshParam.Size_X;
    sC(f,e)=1;
    e=2*i-1 +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=2*i   +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=f+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=f+1+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end
i=MeshParam.Size_X;
f=i+MeshParam.coarse_FNum_former-MeshParam.Size_X;
e=i+MeshParam.coarse_ENum_X_former-MeshParam.Size_X;
sC(f,e)=1;
e=2*i-1 +MeshParam.coarse_FNum_former;
sC(f,e)=-1;
e=2*i   +MeshParam.coarse_FNum_former;
sC(f,e)=-1;
e=f+MeshParam.ENum_XandDiag;
sC(f,e)=-1;
e=f-(MeshParam.Size_X-1)+MeshParam.ENum_XandDiag;
sC(f,e)=1;

% fine-quadrilateral zone 1:main
for j=MeshParam.Fine_Y_from:MeshParam.Triangle_Y_from-2
    for i=1:MeshParam.Size_X-1
        f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
            +MeshParam.coarse_FNum_former;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=1;
        e=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from) + 2*i-1 ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 + 1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
        
        f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
            +MeshParam.coarse_FNum_former;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=1;
        e=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from) + 2*i ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i + 1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
    end
    i=MeshParam.Size_X;
    f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_FNum_former;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 +1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_FNum_former;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) +1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% fine-quadrilateral zone 1: the row next to the triangular zone
j=MeshParam.Triangle_Y_from-1;
for i=1:MeshParam.Size_X-1
    f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_FNum_former;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e= 4*i-2 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_FNum_former;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e=4*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i + 1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end
i=MeshParam.Size_X;
f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
    +MeshParam.coarse_FNum_former;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
    +MeshParam.coarse_ENum_X_former;
sC(f,e)=1;
e= 4*i-2 ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=1;

f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
    +MeshParam.coarse_FNum_former;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
    +MeshParam.coarse_ENum_X_former;
sC(f,e)=1;
e=4*i-1 ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 1 ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=1;

% triangular zone: odd numbered row
for j=MeshParam.Triangle_Y_from:2:MeshParam.Triangle_Y_to-1
    for i=1:MeshParam.Size_X-1
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
    end
    i=MeshParam.Size_X;
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1-2*MeshParam.Size_X ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% triangular zone: even numbered row
for j=MeshParam.Triangle_Y_from+1:2:MeshParam.Triangle_Y_to-2
    for i=1:MeshParam.Size_X-1
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;        

        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;        
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=1;
        e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        sC(f,e)=1;  

    end
    i=MeshParam.Size_X;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1-2*MeshParam.Size_X...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% triangular zone: even numbered row, next to the boundary 
j=MeshParam.Triangle_Y_to;
for i=1:MeshParam.Size_X-1
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e= 2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=-1;   
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e= 2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=1;
    e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
end
i=MeshParam.Size_X;

f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=1;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i-1 ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=-1;

f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=1;
e= 2*i-1 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=1;

f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=1;
e= 2*i ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=-1;

f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=1;
e=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
sC(f,e)=-1;
e=2*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+2*i+1-2*MeshParam.Size_X...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
sC(f,e)=1;


% fine-quadrilateral zone 2: main
for j=MeshParam.Triangle_Y_to+1:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X-1
        f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;        
        sC(f,e)=1;
        e=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 + 1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
       
        f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        sC(f,e)=1;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i + 2*MeshParam.Size_X ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i + 1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
    end
    i=MeshParam.Size_X;
    
    f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=1;
    e=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 + 1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
    
    f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=1;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i + 2*MeshParam.Size_X ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% coarce-quadrilateral zone 2 : the row right outside the boundary
%j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Size_X-1
    f=i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
    e=2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag ... 
        +MeshParam.fine_ENum_X_latter-2*MeshParam.Size_X;
    sC(f,e)=1;
    e=2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag ... 
        +MeshParam.fine_ENum_X_latter-2*MeshParam.Size_X;
    sC(f,e)=1;
    e=i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
    sC(f,e)=-1;
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=i+1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end
i=MeshParam.Size_X;
f=i ...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
e=2*i-1 ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag ...
    +MeshParam.fine_ENum_X_latter-2*MeshParam.Size_X;
sC(f,e)=1;
e=2*i ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag ...
    +MeshParam.fine_ENum_X_latter-2*MeshParam.Size_X;
sC(f,e)=1;
e=i ...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
sC(f,e)=-1;
e=i ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
sC(f,e)=-1;
e=i-(MeshParam.Size_X-1) ...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
sC(f,e)=1;

% coarce-quadrilateral zone 2: main
for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y-1
    for i=1:MeshParam.Size_X-1
        f=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
        e=MeshParam.Size_X*(j-2-MeshParam.Fine_Y_to) +i...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
        sC(f,e)=1;
        e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to) +i...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
        sC(f,e)=-1;
        e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
        sC(f,e)=-1;
        e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i +1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
        sC(f,e)=1;
    end
    i=MeshParam.Size_X;
    
    f=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
    e=MeshParam.Size_X*(j-2-MeshParam.Fine_Y_to) +i...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
    sC(f,e)=1;
    e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to) +i...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
    sC(f,e)=-1;
    e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i-(MeshParam.Size_X-1)...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end

% coarce-quadrilateral zone 2: the last row
j=MeshParam.Size_Y;
for i=1:MeshParam.Size_X-1
    f=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
    e=MeshParam.Size_X*(j-2-MeshParam.Fine_Y_to) +i...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
    sC(f,e)=1;
    e=i;
    sC(f,e)=-1;
    e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=-1;
    e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i+1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    sC(f,e)=1;
end
i=MeshParam.Size_X;

f=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
    +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
e=MeshParam.Size_X*(j-2-MeshParam.Fine_Y_to) +i...
    +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
sC(f,e)=1;
e=i;
sC(f,e)=-1;
e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
sC(f,e)=-1;
e=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i-(MeshParam.Size_X-1)...
    +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
sC(f,e)=1;

%% initialize sG

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X-1
        e=MeshParam.Size_X*(j-1)+i;
        n=MeshParam.Size_X*(j-1)+i;
        sG(e,n)=-1;
        n=MeshParam.Size_X*(j-1)+i+1;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;
    e=MeshParam.Size_X*(j-1)+i;
    n=MeshParam.Size_X*(j-1)+i;
    sG(e,n)=-1;
    n=MeshParam.Size_X*(j-1)+1;
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_from:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X-1
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i+1 ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+1 ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
end

for j=MeshParam.Triangle_Y_from:2:MeshParam.Triangle_Y_to-1
    for i=1:MeshParam.Size_X-1
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
end

for j=MeshParam.Triangle_Y_from+1:2:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X-1
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
    sG(e,n)=1;
end

for j=MeshParam.Triangle_Y_to+1:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Size_X-1
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=1;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
    sG(e,n)=1;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
    sG(e,n)=-1;
    n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X-1
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
        n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
        sG(e,n)=-1;
        n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i+1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
        sG(e,n)=1;
    end
    i=MeshParam.Size_X;     
    e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
    n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
    sG(e,n)=-1;
    n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
    sG(e,n)=1;
end

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-1)+i ...
            +MeshParam.ENum_XandDiag;
        n=MeshParam.Size_X*(j-1)+i;
        sG(e,n)=-1;
        n=MeshParam.Size_X*(j  )+i;
        sG(e,n)=1;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=1:MeshParam.Size_X
    e=MeshParam.Size_X*(j-1)+i ...
        +MeshParam.ENum_XandDiag;
    n=MeshParam.Size_X*(j-1)+i;
    sG(e,n)=-1;
    n=2*i-1 ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_from:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
    end
end

for j=MeshParam.Triangle_Y_from:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former;
        sG(e,n)=1;
        
    end
end

for j=MeshParam.Triangle_Y_to+1:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=1;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        n=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=-1;
        n=2*MeshParam.Size_X*(j+1-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum;
        sG(e,n)=1;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Size_X
    e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    n=2*i-1 ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum ...
        +MeshParam.fine_NNum_latter-2*MeshParam.Size_X;
    sG(e,n)=-1;
    n=MeshParam.Size_X*(j+1-MeshParam.Fine_Y_to-2)+i ...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y-1
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
        n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
        sG(e,n)=-1;
        n=MeshParam.Size_X*(j+1-MeshParam.Fine_Y_to-2)+i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
        sG(e,n)=1;
    end
end

j=MeshParam.Size_Y;
for i=1:MeshParam.Size_X
    e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    n=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i...
        +MeshParam.coarse_NNum_former+MeshParam.fine_NNum_former+MeshParam.triangle_NNum+MeshParam.fine_NNum_latter;
    sG(e,n)=-1;
    n=i;
    sG(e,n)=1;
end

%%

[row_sC,col_sC]=find(sC);
[row_sG,col_sG]=find(sG);

%% initializing edge_vector

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-1)+i;
        edge_vector(e,1)=1.0;
        edge_vector(e,2)=0.0;
    end
end

j=MeshParam.Fine_Y_from;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=1.0/6.0;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=-1.0/6.0;
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
               
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
    end
end

j=MeshParam.Triangle_Y_from;
for i=1:MeshParam.Size_X
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=0.5;
    
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=-0.5*MeshParam.deltatriangle/(0.5-MeshParam.deltatriangle);  

    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=0.5*MeshParam.deltatriangle/(0.5-MeshParam.deltatriangle);
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=-0.5;  
end

for j=MeshParam.Triangle_Y_from+2:2:MeshParam.Triangle_Y_to-1
    for i=1:MeshParam.Size_X
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.5;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=-0.5;
    end
end

for j=MeshParam.Triangle_Y_from+1:2:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;

        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=-0.5;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.5;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
    end
end

j=MeshParam.Triangle_Y_to+1;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=0.5*MeshParam.deltatriangle/(0.5-MeshParam.deltatriangle);
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=-0.5*MeshParam.deltatriangle/(0.5-MeshParam.deltatriangle);
end

for j=MeshParam.Triangle_Y_to+2:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        edge_vector(e,1)=0.5;
        edge_vector(e,2)=0.0;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=-1.0/6.0;
    
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    edge_vector(e,1)=0.5;
    edge_vector(e,2)=1.0/6.0;
end


for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
        edge_vector(e,1)=1.0;
        edge_vector(e,2)=0.0;
    end
end

% Y edges

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-1)+i ...
            +MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=1.0;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=1:MeshParam.Size_X
    e=MeshParam.Size_X*(j-1)+i ...
        +MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=1.0-MeshParam.deltaboundary;
end

j=MeshParam.Fine_Y_from;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=0.5+MeshParam.deltaboundary;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=0.5-1.0/6.0+MeshParam.deltaboundary;
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;
        
    end    
end

% adjustment for row j=MeshParam.Triangle_Y_from-1
j=MeshParam.Triangle_Y_from-1;
for i=1:MeshParam.Size_X    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
    edge_vector(e,2)=edge_vector(e,2)-0.5*MeshParam.deltatriangle*(0.5-MeshParam.deltatriangle);
end

for j=MeshParam.Triangle_Y_from:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;  
    end
end

% adjustment for row j=MeshParam.Triangle_Y_from;
j=MeshParam.Triangle_Y_from;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    edge_vector(e,2)=edge_vector(e,2)+0.5*MeshParam.deltatriangle*(0.5-MeshParam.deltatriangle);
end

% adjustment for row j=MeshParam.Triangle_Y_to;
j=MeshParam.Triangle_Y_to;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
    edge_vector(e,2)=edge_vector(e,2)+0.5*MeshParam.deltatriangle*(0.5-MeshParam.deltatriangle);
end

j=MeshParam.Triangle_Y_to+1;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=0.5;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=0.5-0.5*MeshParam.deltatriangle*(0.5-MeshParam.deltatriangle);
end

for j=MeshParam.Triangle_Y_to+2:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        edge_vector(e,1)=0.0;
        edge_vector(e,2)=0.5;
    end
end

% adjustment for row j=MeshParam.Fine_Y_to-1
j=MeshParam.Fine_Y_to;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    edge_vector(e,2)=edge_vector(e,2)+MeshParam.deltaboundary;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
    edge_vector(e,2)=edge_vector(e,2)-1.0/6.0+MeshParam.deltaboundary;
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Size_X
    e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=1.0-MeshParam.deltaboundary;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
       edge_vector(e,1)=0.0;
       edge_vector(e,2)=1.0;
    end
end

%% initializing tilde_edge_vector

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-1)+i;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=1.0;
    end
end

j=MeshParam.Fine_Y_from;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former;
    tilde_edge_vector(e,1)=-0.25;
    tilde_edge_vector(e,2)=0.75;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
        +MeshParam.coarse_ENum_X_former;
    tilde_edge_vector(e,1)=0.25;
    tilde_edge_vector(e,2)=0.75;
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5;
    end
end

j=MeshParam.Triangle_Y_from;
for i=1:MeshParam.Size_X
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    tilde_edge_vector(e,1)=-2*MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;  

    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    tilde_edge_vector(e,1)=MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=0.5-MeshParam.deltatriangle;

    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    tilde_edge_vector(e,1)=-MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=0.5-MeshParam.deltatriangle;  
    
    e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
    tilde_edge_vector(e,1)=2*MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;  
end

for j=MeshParam.Triangle_Y_from+2:2:MeshParam.Triangle_Y_to-1
    for i=1:MeshParam.Size_X
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=-2*MeshParam.deltatriangle;
        tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5-2*MeshParam.deltatriangle;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5-2*MeshParam.deltatriangle;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=2*MeshParam.deltatriangle;
        tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;
    end
end

for j=MeshParam.Triangle_Y_from+1:2:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5-2*MeshParam.deltatriangle;

        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=2*MeshParam.deltatriangle;
        tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;
        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=-2*MeshParam.deltatriangle;
        tilde_edge_vector(e,2)=2*MeshParam.deltatriangle;

        
        e=4*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5-2*MeshParam.deltatriangle;
    end
end


j=MeshParam.Triangle_Y_to+1;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    tilde_edge_vector(e,1)=-MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=0.5-MeshParam.deltatriangle;
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    tilde_edge_vector(e,1)=MeshParam.deltatriangle;
    tilde_edge_vector(e,2)=0.5-MeshParam.deltatriangle;
end

for j=MeshParam.Triangle_Y_to+2:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5;
        
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=0.5;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Size_X
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    tilde_edge_vector(e,1)=0.25;
    tilde_edge_vector(e,2)=0.75;
    
    
    e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
        +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag;
    tilde_edge_vector(e,1)=-0.25;
    tilde_edge_vector(e,2)=0.75;
end


for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)+i ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X_former+MeshParam.triangle_ENum_XandDiag+MeshParam.fine_ENum_X_latter;
        tilde_edge_vector(e,1)=0.0;
        tilde_edge_vector(e,2)=1.0;
    end
end

% Y edges

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-1)+i ...
            +MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-1.0;
        tilde_edge_vector(e,2)=0.0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-0.5;
        tilde_edge_vector(e,2)=0.0;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-0.5;
        tilde_edge_vector(e,2)=0.0;
        
    end    
end

for j=MeshParam.Triangle_Y_from:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-(0.5-2*MeshParam.deltatriangle);
        tilde_edge_vector(e,2)=0.0;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_from)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-(0.5-2*MeshParam.deltatriangle);
        tilde_edge_vector(e,2)=0.0;  
    end
end

for j=MeshParam.Triangle_Y_to+1:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i-1 ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-0.5;
        tilde_edge_vector(e,2)=0.0;
        
        e=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1)+2*i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.ENum_XandDiag;
        tilde_edge_vector(e,1)=-0.5;
        tilde_edge_vector(e,2)=0.0;
    end
end

for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1)+i ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y_former+MeshParam.triangle_ENum_Y+MeshParam.fine_ENum_Y_latter+MeshParam.ENum_XandDiag;
       tilde_edge_vector(e,1)=-1.0;
       tilde_edge_vector(e,2)=0.0;
    end
end

%% initializing tilde_node_position

% coarce-quadrilateral zone 1: main
for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        f=MeshParam.Size_X*(j-1)+i;
        tilde_node_position(f,1)=-0.5+i;
        tilde_node_position(f,2)=-0.5+j;
    end
end

% fine-quadrilateral zone 1:main
for j=MeshParam.Fine_Y_from:MeshParam.Triangle_Y_from-1
    for i=1:MeshParam.Size_X
        f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i-1 ...
            +MeshParam.coarse_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1);
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from);        

        f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from) + 2*i ...
            +MeshParam.coarse_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i);
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from);
    end
end

% triangular zone: odd numbered row
for j=MeshParam.Triangle_Y_from:2:MeshParam.Triangle_Y_to-1
    for i=1:MeshParam.Size_X
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1)-MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)+MeshParam.deltatriangle;

        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1)+MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)-MeshParam.deltatriangle;

        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i)-MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)-MeshParam.deltatriangle;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i)+MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)+MeshParam.deltatriangle;
    end
end

% triangular zone: even numbered row
for j=MeshParam.Triangle_Y_from+1:2:MeshParam.Triangle_Y_to
    for i=1:MeshParam.Size_X
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1)-MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)-MeshParam.deltatriangle;

        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1)+MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)+MeshParam.deltatriangle;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i)-MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)+MeshParam.deltatriangle;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        tilde_node_position(f,1)=-0.25+0.5*(2*i)+MeshParam.deltatriangle;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from)-MeshParam.deltatriangle;
    end
end

% fine-quadrilateral zone 2: main
for j=MeshParam.Triangle_Y_to+1:MeshParam.Fine_Y_to
    for i=1:MeshParam.Size_X
        f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
        tilde_node_position(f,1)=-0.25+0.5*(2*i-1);
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from);

        
        f=2*MeshParam.Size_X*(j-MeshParam.Triangle_Y_to-1) + 2*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum;
        tilde_node_position(f,1)=-0.25+0.5*(2*i);
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+0.25+0.5*(j-MeshParam.Fine_Y_from);
    end
end

% coarce-quadrilateral zone 2: main
for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        f=MeshParam.Size_X*(j-1-MeshParam.Fine_Y_to)+i...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter;
        tilde_node_position(f,1)=-0.5+i;
        tilde_node_position(f,2)=MeshParam.Fine_Y_from-1+(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0+0.5+(j-1-MeshParam.Fine_Y_to);
    end
end


%% initializing denominators
for f=1:MeshParam.coarse_FNum_former
    denominator_face(f)=1;
end
for f=MeshParam.coarse_FNum_former+1:MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter
    denominator_face(f)=denominator_obi;
end
for f=MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former+MeshParam.triangle_FNum+MeshParam.fine_FNum_latter+1:FNum
    denominator_face(f)=1;
end


% genrate denominator_edge,denominator_node
% can be prosessed in parallel
for e=1:ENum
    denominator_edge(e)=0;
end
for kk=1:size(row_sC,1)
    f=row_sC(kk);
    e=col_sC(kk);
    %disp(e)
    %i=0;
    if denominator_edge(e)<denominator_face(f)
        denominator_edge(e)=denominator_face(f);
        %i=i+1;
    end
end
for kk=1:size(row_sG,1)
    e=row_sG(kk);
    n=col_sG(kk);
    if denominator_node(n)<denominator_edge(e)
        denominator_node(n)=denominator_edge(e);
    end
end
% end denominator_edge, denominator_node

%% generate st faces
%disp('generating st faces')
p=1;
%Omega=1;
for f=1:FNum
    first_p_for_f(f)=p;
    %first_Omega_for_f(f)=Omega;
    %l(p)=0;
    p=p+denominator_face(f)+1;
    %Omega=Omega+denominator_face(f);
end
%OmegaNum=Omega-1;
%Omega=1;
for e=1:ENum
    first_p_for_e(e)=p;
    %first_Omega_for_e(e)=Omega;
    %l(p)=0;
    p=p+denominator_edge(e)+1;
    %Omega=Omega+denominator_edge(e);
end
%OmegaDualNum=Omega-1;
PNum=p-1;

MeshNum.P=PNum;

%% END
%sC,sG
denominator.n=denominator_node;
denominator.e=denominator_edge;
denominator.f=denominator_face;

for e=1:MeshNum.E
    edgevec.prim(e).vec=edge_vector(e,:).';
end
for e=1:MeshNum.E
    edgevec.dual(e).vec=tilde_edge_vector(e,:).';
end
first_p.f=first_p_for_f;
first_p.e=first_p_for_e;
%first_omega.f=first_Omega_for_f;
%first_omega.e=first_Omega_for_e;

%tilde_node_position;
 
disp('GenerateMesh_triangular:STOPPED')

end