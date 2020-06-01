function [sC,sG,UpdateNum,edgevec,first_pIdx,tilde_f,MeshNum,MeshParam] ...
    = GenerateMesh_squarefaces_squaresubgrid(MeshParam)

global DIM
global DISPDEBUGGINGMESSAGE 
  
if DISPDEBUGGINGMESSAGE 
    disp('GenerateMesh_square_belt:CALLED')
    disp('Initializing Spatial Mesh Information ')
end
%global MeshNum

MeshParam.FacePerRow_Subgrid=...
    4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1);

MeshParam.XEdgePerRow_Subgrid=...
    4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1);

MeshParam.YEdgePerRow_Subgrid=...
    4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)+2 ...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)-1;

MeshParam.NodePerRow_Subgrid=...
    4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)+2 ...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)-1;
    

MeshParam.coarse_FNum_former=MeshParam.Size_X*(MeshParam.Fine_Y_from-1);
MeshParam.fine_FNum=MeshParam.FacePerRow_Subgrid*(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1);
MeshParam.coarse_FNum_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to);
MeshNum.F=MeshParam.coarse_FNum_former+MeshParam.fine_FNum+MeshParam.coarse_FNum_latter;
    
MeshParam.coarse_ENum_X_former=MeshParam.coarse_FNum_former;
MeshParam.fine_ENum_X=MeshParam.XEdgePerRow_Subgrid*(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)...
    +2*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1);
MeshParam.coarse_ENum_X_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to-1);
MeshParam.ENum_X=MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X+MeshParam.coarse_ENum_X_latter;

MeshParam.coarse_ENum_Y_former=MeshParam.coarse_FNum_former;
MeshParam.fine_ENum_Y=MeshParam.YEdgePerRow_Subgrid*(MeshParam.Fine_Y_to+1-MeshParam.Fine_Y_from);
MeshParam.coarse_ENum_Y_latter=MeshParam.coarse_FNum_latter;
MeshParam.ENum_Y=MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y+MeshParam.coarse_ENum_Y_latter;
MeshNum.E=MeshParam.ENum_X+MeshParam.ENum_Y;

MeshParam.coarse_NNum_former=MeshParam.Size_X*(MeshParam.Fine_Y_from-1);
MeshParam.fine_NNum=MeshParam.NodePerRow_Subgrid*(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)...
    +2*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)+1 ...
    +MeshParam.Size_X-(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)-1;
MeshParam.coarse_NNum_latter=MeshParam.Size_X*(MeshParam.Size_Y-MeshParam.Fine_Y_to-1);
MeshNum.N=MeshParam.coarse_NNum_former+MeshParam.fine_NNum+MeshParam.coarse_NNum_latter;

%% check
if DISPDEBUGGINGMESSAGE
    disp('Check; [MeshNum.F, MeshNum.E, MeshNum.N]=')
    disp([MeshNum.F MeshNum.E MeshNum.N])
end
%% Allocate Arrays
sC=sparse(MeshNum.F,MeshNum.E);
sG=sparse(MeshNum.E,MeshNum.N);
UpdateNum.n=zeros(MeshNum.N,1);
UpdateNum.e=zeros(MeshNum.E,1);
UpdateNum.f=zeros(MeshNum.F,1);

Dummy=cell(MeshNum.E,1);
for e=1:MeshNum.E
    Dummy{e}=zeros(DIM,1);
end
edgevec.prim=struct('vec',Dummy);
edgevec.dual=struct('vec',Dummy);
clearvars Dummy

first_pIdx.f=zeros(MeshNum.F,1);
first_pIdx.e=zeros(MeshNum.E,1);

Dummy=cell(MeshNum.F,1);
for f=1:MeshNum.F
    Dummy{f}=zeros(DIM,1);
end
tilde_f=struct('position',Dummy);

%% initialize sC

% coarce-quadrilateral zone 1: main
for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        f=MeshParam.Size_X*(j-1)+i;
        e=f;
        sC(f,e)=1;
        e=f+MeshParam.Size_X;
        sC(f,e)=-1;
        e=f+MeshParam.ENum_X;
        sC(f,e)=-1;
        if i==MeshParam.Size_X
            e=f-(MeshParam.Size_X-1)+MeshParam.ENum_X;
        else
            e=f+1+MeshParam.ENum_X;
        end
        sC(f,e)=1;
    end
end

% coarce-quadrilateral zone 1: the row right outside the boundary
for i=1:MeshParam.Fine_X_from-1
    f=i+MeshParam.coarse_FNum_former-MeshParam.Size_X;
    e=f;
    sC(f,e)=1;
    e=f+MeshParam.Size_X;
    sC(f,e)=-1;
    e=f+MeshParam.ENum_X;
    sC(f,e)=-1;
    e=f+1+MeshParam.ENum_X;
    sC(f,e)=1;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    f=i+MeshParam.coarse_FNum_former-MeshParam.Size_X;
    e=i+MeshParam.coarse_ENum_X_former-MeshParam.Size_X;
    sC(f,e)=1;
    e=4*(i-MeshParam.Fine_X_from+1)-3 +MeshParam.Fine_X_from-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=4*(i-MeshParam.Fine_X_from+1)-1 +MeshParam.Fine_X_from-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=f+MeshParam.ENum_X;
    sC(f,e)=-1;
    e=f+1+MeshParam.ENum_X;
    sC(f,e)=1;
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    f=i+MeshParam.coarse_FNum_former-MeshParam.Size_X;
    e=f;
    sC(f,e)=1;
    e= i-MeshParam.Fine_X_to ...
        +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=f+MeshParam.ENum_X;
    sC(f,e)=-1;
    if i==MeshParam.Size_X
        e=f-(MeshParam.Size_X-1)+MeshParam.ENum_X;
    else
        e=f+1+MeshParam.ENum_X;
    end
    sC(f,e)=1;
end

% fine-quadrilateral zone:main
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-2
        f=i  +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        e=i  +MeshParam.XEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=1;
        e=i  +MeshParam.XEdgePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=-1;
        e=i  +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
        sC(f,e)=-1;
        e=i+1+MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
        sC(f,e)=1;
    end
end
i=MeshParam.Fine_X_from-1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    f=i  +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_FNum_former;
    e=i  +MeshParam.XEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e=i  +MeshParam.XEdgePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=-1;
    e=i  +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    sC(f,e)=-1;
    e=i+1+MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    sC(f,e)=1;
    e=i+2+MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    sC(f,e)=1;
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localfaceIdx=1:4
            f=localfaceIdx...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...      
                +MeshParam.coarse_FNum_former;
            e=localfaceIdx...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...      
                +MeshParam.coarse_ENum_X_former;
            sC(f,e)=1;
            if localfaceIdx==1 || localfaceIdx==3
                e=localfaceIdx+1 ...
                    +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                    +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_ENum_X_former;
                sC(f,e)=-1;
            elseif j==MeshParam.Fine_Y_to
                e=localfaceIdx/2 ...
                    +2*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                    +MeshParam.XEdgePerRow_Subgrid*(MeshParam.Fine_Y_to+1-MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_ENum_X_former;
                sC(f,e)=-1;
            else    
                e=localfaceIdx-1 ...
                    +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                    +MeshParam.XEdgePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_ENum_X_former;
                sC(f,e)=-1;
            end
            e=localfaceIdx ...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
            sC(f,e)=-1;
            e=localfaceIdx+2 ...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
            sC(f,e)=1;            
        end
    end
end
i=MeshParam.Fine_X_to+1; 
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    f=i-MeshParam.Fine_X_to ...
        +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
        +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_FNum_former;
    e=i-MeshParam.Fine_X_to ...
        +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
        +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    if j==MeshParam.Fine_Y_to
        e= i -MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +2*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.XEdgePerRow_Subgrid *(j+1-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
    else
        e= i -MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.XEdgePerRow_Subgrid *(j+1-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
    end
    sC(f,e)=-1;
    e=1 +i-(MeshParam.Fine_X_to+1) ...
        +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
        +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    sC(f,e)=-1;
    e=2 +i-(MeshParam.Fine_X_to+1) ...
        +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
        +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    sC(f,e)=-1;
    if i==MeshParam.Size_X
        e=1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    else
        e= i+1-(MeshParam.Fine_X_to+1) ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +2 ...
            +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
    end
    sC(f,e)=1;
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        f=i-MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        e=i-MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        sC(f,e)=1;
        if j==MeshParam.Fine_Y_to
            e=i  -MeshParam.Fine_X_to ...
                +MeshParam.Fine_X_from-1 +2*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
                +MeshParam.XEdgePerRow_Subgrid *(MeshParam.Fine_Y_to+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_X_former;
        else
            e=i  -MeshParam.Fine_X_to ...
                +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
                +MeshParam.XEdgePerRow_Subgrid *(j+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_X_former;
        end
        sC(f,e)=-1;
        e=i -(MeshParam.Fine_X_to+1) ...
            +MeshParam.Fine_X_from-1 ...
            +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +2 ...
            +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
        sC(f,e)=-1;
     %   disp(e)
        if i==MeshParam.Size_X
            e=1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
        else
            e=i+1-(MeshParam.Fine_X_to+1) ...
                +MeshParam.Fine_X_from-1 ...
                +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
                +2 ...
                +MeshParam.YEdgePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.ENum_X;
        end
        sC(f,e)=1;
    end
end

% coarce-quadrilateral zone: the row right outside the boundary
j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    f=i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
    e=i ...
        +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    if j==MeshParam.Size_Y
        e=i;
    else
        e=i ...
            +MeshParam.Size_X*(j+1-(MeshParam.Fine_Y_to+2)) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
    end
    sC(f,e)=-1;
    e=i ...
        +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    sC(f,e)=-1;
    e=i+1 ...
        +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    sC(f,e)=1;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    f=i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
    e=1 ...
        +2*(i -MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    e=2 ...
        +2*(i -MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    if j==MeshParam.Size_Y
        e=i;
    else
        e=i ...
            +MeshParam.Size_X*(j+1-(MeshParam.Fine_Y_to+2)) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
    end
    sC(f,e)=-1;
    e=i ...
        +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    sC(f,e)=-1;
    e=i+1 ...
        +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    sC(f,e)=1;
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    f=i ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
    e=i-MeshParam.Fine_X_to...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    sC(f,e)=1;
    if j==MeshParam.Size_Y
        e=i;
    else
        e=i ...
            +MeshParam.Size_X*(j+1-(MeshParam.Fine_Y_to+2)) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
    end
    sC(f,e)=-1;
    e=i ...
        +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    sC(f,e)=-1;
    if i==MeshParam.Size_X
        e=1 ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
    else
        e=i+1 ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
    end
    sC(f,e)=1;
end

% coarce-quadrilateral zone 2: main
for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        f=i ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1))...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
        e=i ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+2))...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
        sC(f,e)=1;
        if j==MeshParam.Size_Y
            e=i;
        else
            e=i ...
                +MeshParam.Size_X*(j+1-(MeshParam.Fine_Y_to+2)) ...
                +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
        end
        sC(f,e)=-1;
        e=i ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
        sC(f,e)=-1;
        if i==MeshParam.Size_X
            e=1 ...
                +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
                +MeshParam.ENum_X;
        else
            e=i+1 ...
                +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1)) ...
                +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
                +MeshParam.ENum_X;
        end
        sC(f,e)=1;
    end
end

%% initialize sG

% X edges

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1);
        n=i+MeshParam.Size_X*(j-1);
        sG(e,n)=-1;
        if i==MeshParam.Size_X
            n=1  +MeshParam.Size_X*(j-1);
        else
            n=i+1+MeshParam.Size_X*(j-1);
        end
        sG(e,n)=1;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        n=i...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=i+1 ...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;        
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localfaceIdx=1:4
            e=localfaceIdx ...
                +4*(i-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_ENum_X_former;        
            n=localfaceIdx ...
                +4*(i-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;        
            sG(e,n)=-1;
            n=localfaceIdx+2 ...
                +4*(i-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;        
            sG(e,n)=1;
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        e=i  -MeshParam.Fine_X_to ...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;
        if i==MeshParam.Fine_X_to+1
            n=1  ...
                +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;
        else
            n=i  -(MeshParam.Fine_X_to+1) ...
                +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +2 ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;
        end
        sG(e,n)=-1;
        if i==MeshParam.Size_X
            n=1   ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;
        else
            n=i+1-(MeshParam.Fine_X_to+1) ...
                +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +2 ...
                +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_NNum_former;
        end
        sG(e,n)=1;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    e=i ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    n=i...
        +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    n=i+1 ...
        +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +2*(i-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;
        n=localedgeIdx ...
            +2*(i-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=localedgeIdx+1 ...
            +2*(i-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
    end
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    e=i  -MeshParam.Fine_X_to ...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_ENum_X_former;
    n=i  -MeshParam.Fine_X_to ...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    if i==MeshParam.Size_X
        n=1   ...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
    else
        n=i+1-MeshParam.Fine_X_to ...
            +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_NNum_former;
    end
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
        n=i...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
        sG(e,n)=-1;
       if i==MeshParam.Size_X
         n=1 ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
       else
        n=i+1 ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2)...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
       end
        sG(e,n)=1;
    end
end

%% Y edges

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1) ...
            +MeshParam.ENum_X;
        n=i+MeshParam.Size_X*(j-1);
        sG(e,n)=-1;
        n=i+MeshParam.Size_X*(j  );
        sG(e,n)=1;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=1:MeshParam.Fine_X_from-1
        e=i+MeshParam.Size_X*(j-1) ...
            +MeshParam.ENum_X;
        n=i+MeshParam.Size_X*(j-1);
        sG(e,n)=-1;
        n=i+MeshParam.Size_X*(j  );
        sG(e,n)=1;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to+1
    e=i+MeshParam.Size_X*(j-1) ...
        +MeshParam.ENum_X;
    n=i+MeshParam.Size_X*(j-1);
    sG(e,n)=-1;
    n=1+4*(i-MeshParam.Fine_X_from) ...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
end
for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
    e=i+MeshParam.Size_X*(j-1) ...
        +MeshParam.ENum_X;
    n=i+MeshParam.Size_X*(j-1);
    sG(e,n)=-1;
    n=i-(MeshParam.Fine_X_to+1)...
        +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
        +2 ...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        n=i...
            +MeshParam.NodePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        n=i  ...
            +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=1;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localedgeIdx=1:4
            e=localedgeIdx ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former ...
                +MeshParam.ENum_X;
            n=localedgeIdx ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former;
            sG(e,n)=-1;
            if localedgeIdx==1 || localedgeIdx==3
                n=localedgeIdx+1 ...
                    +4*(i-MeshParam.Fine_X_from) ...
                    +MeshParam.Fine_X_from-1 ...
                    +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_NNum_former;
            elseif j==MeshParam.Fine_Y_to
                n=localedgeIdx/2 ...
                    +2*(i-MeshParam.Fine_X_from) ...
                    +MeshParam.Fine_X_from-1 ...
                    +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_NNum_former;
            else
                n=localedgeIdx-1 ...
                    +4*(i-MeshParam.Fine_X_from) ...
                    +MeshParam.Fine_X_from-1 ...
                    +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                    +MeshParam.coarse_NNum_former;
            end
            sG(e,n)=1;
        end
    end
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +4*(i-MeshParam.Fine_X_from) ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
      %  disp(e)
        n=localedgeIdx ...
            +4*(i-MeshParam.Fine_X_from) ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        if localedgeIdx==1
            n=localedgeIdx+1 ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former; 
        elseif j==MeshParam.Fine_Y_to
            n=1 ...
                +2*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former;
        else
            n=localedgeIdx/2 ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former;
        end
        sG(e,n)=1;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        e=i-(MeshParam.Fine_X_to+1)...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
            +2 ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        n=i-(MeshParam.Fine_X_to+1)...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
            +2 ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_NNum_former;
        sG(e,n)=-1;
        if j==MeshParam.Fine_Y_to
            n=i-MeshParam.Fine_X_to ...
                +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former;
        else
            n=i-(MeshParam.Fine_X_to+1)...
                +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
                +2 ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.NodePerRow_Subgrid*(j+1-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_NNum_former;
        end
        sG(e,n)=1;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    n=i...
        +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    if j==MeshParam.Size_Y
        n=i;
    else
        n=i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
    end
    sG(e,n)=1;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to+1
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    n=1 ...
        +2*(i-MeshParam.Fine_X_from) ...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    if j==MeshParam.Size_Y
        n=i;
    else
        n=i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
    end
    sG(e,n)=1;
end
for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    n=i-MeshParam.Fine_X_to...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.NodePerRow_Subgrid*(j  -MeshParam.Fine_Y_from)...
        +MeshParam.coarse_NNum_former;
    sG(e,n)=-1;
    if j==MeshParam.Size_Y
        n=i;
    else
        n=i ...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
    end
    sG(e,n)=1;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
        n=i...
            +MeshParam.Size_X*(j  -MeshParam.Fine_Y_to-2)...
            +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
        sG(e,n)=-1;
        if j~=MeshParam.Size_Y
            n=i ...
                +MeshParam.Size_X*(j+1-MeshParam.Fine_Y_to-2)...
                +MeshParam.coarse_NNum_former+MeshParam.fine_NNum;
        else
            n=i;
        end
        sG(e,n)=1;
    end
end

%%

[row_sC,col_sC]=find(sC);
[row_sG,col_sG]=find(sG);

%% initializing edge_vector

% X edges

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1);
        edgevec.prim(e).vec(1)=1.0;
        edgevec.prim(e).vec(2)=0.0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        edgevec.prim(e).vec(1)=1.0;
        edgevec.prim(e).vec(2)=0.0;
        if i==MeshParam.Fine_X_from-1 
            if j==MeshParam.Fine_Y_from
                edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
                edgevec.prim(e).vec(2)=   -MeshParam.deltaboundary;
            else
                edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
            end
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localedgeIdx=1:4
            e=localedgeIdx ...
                +4*(i-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_ENum_X_former;         
            edgevec.prim(e).vec(1)=0.5;
            edgevec.prim(e).vec(2)=0.0;
            if j==MeshParam.Fine_Y_from
                if i==MeshParam.Fine_X_from
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)=0.5    +MeshParam.deltaboundary;
                        edgevec.prim(e).vec(2)=1.0/6.0;
                    elseif localedgeIdx == 2
                        edgevec.prim(e).vec(1)=0.5+MeshParam.deltaboundary-1.0/6.0;
                        edgevec.prim(e).vec(2)=0.0;
                    elseif localedgeIdx == 3
                        edgevec.prim(e).vec(1)= 0.5;
                        edgevec.prim(e).vec(2)=-1.0/6.0;
                    end
                elseif i==MeshParam.Fine_X_to
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 0.5;
                        edgevec.prim(e).vec(2)= 1.0/6.0;
                    elseif localedgeIdx == 3
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary;
                        edgevec.prim(e).vec(2)=-1.0/6.0;
                    elseif localedgeIdx == 4
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.0;
                    end
                else
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 0.5;
                        edgevec.prim(e).vec(2)= 1.0/6.0;
                    elseif localedgeIdx == 3
                        edgevec.prim(e).vec(1)= 0.5;
                        edgevec.prim(e).vec(2)=-1.0/6.0;
                    end
                end
            else
                if i==MeshParam.Fine_X_from
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary;
                        edgevec.prim(e).vec(2)= 0.0;
                    elseif localedgeIdx ==2
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.0;
                    end
                elseif i==MeshParam.Fine_X_to
                    if localedgeIdx == 3
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary;
                        edgevec.prim(e).vec(2)= 0.0;
                    elseif localedgeIdx == 4
                        edgevec.prim(e).vec(1)= 0.5+MeshParam.deltaboundary-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.0;
                    end
                end
            end
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        e=i  -MeshParam.Fine_X_to ...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;
        edgevec.prim(e).vec(1)=1.0;
        edgevec.prim(e).vec(2)=0.0;            
        if i==MeshParam.Fine_X_to+1
            if j==MeshParam.Fine_Y_from
                edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
                edgevec.prim(e).vec(2)=   +MeshParam.deltaboundary;
            else
                edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
                edgevec.prim(e).vec(2)=0.0;
            end
        end
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    e=i ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    edgevec.prim(e).vec(1)=1.0;
    edgevec.prim(e).vec(2)=0.0;
    if i==MeshParam.Fine_X_from-1
        edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=    MeshParam.deltaboundary;
    end
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +2*(i-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;       
        edgevec.prim(e).vec(1)=0.5;
        edgevec.prim(e).vec(2)=0.0;
        if i==MeshParam.Fine_X_from && localedgeIdx ==1
            edgevec.prim(e).vec(1)= 0.5    +MeshParam.deltaboundary;
            edgevec.prim(e).vec(2)=-1.0/6.0;
        elseif i==MeshParam.Fine_X_to && localedgeIdx ==2
            edgevec.prim(e).vec(1)= 0.5    +MeshParam.deltaboundary;
            edgevec.prim(e).vec(2)= 1.0/6.0;
        elseif localedgeIdx ==1
            edgevec.prim(e).vec(1)= 0.5;
            edgevec.prim(e).vec(2)=-1.0/6.0;
        elseif localedgeIdx ==2
            edgevec.prim(e).vec(1)= 0.5;
            edgevec.prim(e).vec(2)= 1.0/6.0;
        end
    end
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    e=i  -MeshParam.Fine_X_to ...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_ENum_X_former;
    %disp(e)
    edgevec.prim(e).vec(1)=1.0;
    edgevec.prim(e).vec(2)=0.0;
    if i==MeshParam.Fine_X_to+1
        edgevec.prim(e).vec(1)=1.0-MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=   -MeshParam.deltaboundary;
    end
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
        edgevec.prim(e).vec(1)=1.0;
        edgevec.prim(e).vec(2)=0.0;
    end
end

%%

% Y edges

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1) ...
            +MeshParam.ENum_X;
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=1:MeshParam.Fine_X_from-1
    e=i+MeshParam.Size_X*(j-1) ...
        +MeshParam.ENum_X;
    edgevec.prim(e).vec(1)=0.0;
    edgevec.prim(e).vec(2)=1.0;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    e=i+MeshParam.Size_X*(j-1) ...
        +MeshParam.ENum_X;
    edgevec.prim(e).vec(1)=0.0;
    edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    if i==MeshParam.Fine_X_from
        edgevec.prim(e).vec(1)=-MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    else
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    end
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    e=i+MeshParam.Size_X*(j-1) ...
        +MeshParam.ENum_X;
    edgevec.prim(e).vec(1)=0.0;
    edgevec.prim(e).vec(2)=1.0;
    if i==MeshParam.Fine_X_to+1
        edgevec.prim(e).vec(1)=    MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localedgeIdx=1:4
            e=localedgeIdx ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former ...
                +MeshParam.ENum_X;
            edgevec.prim(e).vec(1)=0.0;
            edgevec.prim(e).vec(2)=0.5;
            if i==MeshParam.Fine_X_from
                if j==MeshParam.Fine_Y_from
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5    +MeshParam.deltaboundary;
                    elseif localedgeIdx == 2
                        edgevec.prim(e).vec(1)=-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5;
                    elseif localedgeIdx == 3
                        edgevec.prim(e).vec(1)=0.0;
                        edgevec.prim(e).vec(2)=0.5+MeshParam.deltaboundary-1.0/6.0;
                    end
                elseif j==MeshParam.Fine_Y_to
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5;
                    elseif localedgeIdx == 2
                        edgevec.prim(e).vec(1)=-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5    +MeshParam.deltaboundary;
                    elseif localedgeIdx == 4
                        edgevec.prim(e).vec(1)=0.0;
                        edgevec.prim(e).vec(2)=0.5+MeshParam.deltaboundary-1.0/6.0;
                    end
                else
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5;
                    elseif  localedgeIdx == 2
                        edgevec.prim(e).vec(1)=-1.0/6.0;
                        edgevec.prim(e).vec(2)= 0.5;
                    end
                end
            else
                if j==MeshParam.Fine_Y_from
                    if localedgeIdx ==1
                        edgevec.prim(e).vec(1)= 0.0;
                        edgevec.prim(e).vec(2)= 0.5+MeshParam.deltaboundary;
                    elseif localedgeIdx ==3
                        edgevec.prim(e).vec(1)= 0.0;
                        edgevec.prim(e).vec(2)= 0.5+MeshParam.deltaboundary-1.0/6.0;
                    end
                elseif j==MeshParam.Fine_Y_to
                    if localedgeIdx == 3
                        edgevec.prim(e).vec(1)= 0.0;
                        edgevec.prim(e).vec(2)= 0.5+MeshParam.deltaboundary;
                    elseif localedgeIdx == 4
                        edgevec.prim(e).vec(1)= 0.0;
                        edgevec.prim(e).vec(2)= 0.5+MeshParam.deltaboundary-1.0/6.0;
                    end
                end
            end
        end
    end
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +4*(i-MeshParam.Fine_X_from) ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        if j==MeshParam.Fine_Y_from && localedgeIdx==1
                edgevec.prim(e).vec(1)=-1.0/6.0;
                edgevec.prim(e).vec(2)= 0.5 +MeshParam.deltaboundary;
        elseif j==MeshParam.Fine_Y_to && localedgeIdx==2
                edgevec.prim(e).vec(1)= 1.0/6.0;
                edgevec.prim(e).vec(2)= 0.5    +MeshParam.deltaboundary;
        else
            if localedgeIdx==1
                edgevec.prim(e).vec(1)=-1.0/6.0;
                edgevec.prim(e).vec(2)= 0.5;
            elseif localedgeIdx==2
                edgevec.prim(e).vec(1)= 1.0/6.0;
                edgevec.prim(e).vec(2)= 0.5;
            end
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        e=i-(MeshParam.Fine_X_to+1)...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
            +2 ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0;
    end
end

%j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    edgevec.prim(e).vec(1)=0.0;
    edgevec.prim(e).vec(2)=1.0;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to+1
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    if i==MeshParam.Fine_X_from
        edgevec.prim(e).vec(1)=    MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    elseif i==MeshParam.Fine_X_to+1
        edgevec.prim(e).vec(1)=   -MeshParam.deltaboundary;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    else
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0-MeshParam.deltaboundary;
    end
end
for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
    e=i ...
        +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
        +MeshParam.ENum_X;
    edgevec.prim(e).vec(1)=0.0;
    edgevec.prim(e).vec(2)=1.0;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
        edgevec.prim(e).vec(1)=0.0;
        edgevec.prim(e).vec(2)=1.0;
    end
end



%% initializing edgevec.dual

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1);
        edgevec.dual(e).vec(1)=0.0;
        edgevec.dual(e).vec(2)=1.0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_X_former;
        edgevec.dual(e).vec(1)=0.0;
        edgevec.dual(e).vec(2)=1.0;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localfaceIdx=1:4
            e=localfaceIdx ...
                +4*(i-MeshParam.Fine_X_from)...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_ENum_X_former;
            edgevec.dual(e).vec(1)=0.0;
            edgevec.dual(e).vec(2)=0.5;
            if j==MeshParam.Fine_Y_from
                if localfaceIdx ==1
                    edgevec.dual(e).vec(1)=-0.25;
                    edgevec.dual(e).vec(2)= 0.75;
                elseif localfaceIdx ==3
                    edgevec.dual(e).vec(1)= 0.25;
                    edgevec.dual(e).vec(2)= 0.75;
                end
            end
        end
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        e=i  -MeshParam.Fine_X_to ...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;
        edgevec.dual(e).vec(1)=0.0;
        edgevec.dual(e).vec(2)=1.0;
    end
end

j=MeshParam.Fine_Y_to+1;
for i=1:MeshParam.Fine_X_from-1
    e=i ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
        +MeshParam.coarse_ENum_X_former;
    edgevec.dual(e).vec(1)=0.0;
    edgevec.dual(e).vec(2)=1.0;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +2*(i-MeshParam.Fine_X_from)...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
            +MeshParam.coarse_ENum_X_former;
        %disp(e)
        if localedgeIdx ==1
            edgevec.dual(e).vec(1)= 0.25;
            edgevec.dual(e).vec(2)= 0.75;
        elseif localedgeIdx ==2
            edgevec.dual(e).vec(1)=-0.25;
            edgevec.dual(e).vec(2)= 0.75;
        end
    end
end
for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    e=i  -MeshParam.Fine_X_to ...
        +2*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from)...
        +MeshParam.Fine_X_from-1 ...
        +MeshParam.XEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
        +MeshParam.coarse_ENum_X_former;
    edgevec.dual(e).vec(1)=0.0;
    edgevec.dual(e).vec(2)=1.0;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-2) ...
            +MeshParam.coarse_ENum_X_former+MeshParam.fine_ENum_X;
        edgevec.dual(e).vec(1)=0.0;
        edgevec.dual(e).vec(2)=1.0;
    end
end

%% Y edges

for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        e=i+MeshParam.Size_X*(j-1) ...
            +MeshParam.ENum_X;
        edgevec.dual(e).vec(1)=-1.0;
        edgevec.dual(e).vec(2)= 0.0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        e=i...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        edgevec.dual(e).vec(1)=-1.0;
        edgevec.dual(e).vec(2)= 0.0;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localedgeIdx=1:4
            e=localedgeIdx ...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_ENum_Y_former ...
                +MeshParam.ENum_X;
            edgevec.dual(e).vec(1)=-0.5;
            edgevec.dual(e).vec(2)= 0.0;
            if i==MeshParam.Fine_X_from
                if localedgeIdx==1
                    edgevec.dual(e).vec(1)=-0.75;
                    edgevec.dual(e).vec(2)= 0.25;
                elseif localedgeIdx==2
                    edgevec.dual(e).vec(1)=-0.75;
                    edgevec.dual(e).vec(2)=-0.25;
                end
            end
        end
    end
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for localedgeIdx=1:2
        e=localedgeIdx ...
            +4*(i-MeshParam.Fine_X_from) ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        if localedgeIdx==1
            edgevec.dual(e).vec(1)=-0.75;
            edgevec.dual(e).vec(2)=-0.25;
        elseif localedgeIdx==2
            edgevec.dual(e).vec(1)=-0.75;
            edgevec.dual(e).vec(2)= 0.25;
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        e=i-(MeshParam.Fine_X_to+1)...
            +4*(MeshParam.Fine_X_to+1-MeshParam.Fine_X_from) ...
            +2 ...
            +MeshParam.Fine_X_from-1 ...
            +MeshParam.YEdgePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_ENum_Y_former ...
            +MeshParam.ENum_X;
        edgevec.dual(e).vec(1)=-1.0;
        edgevec.dual(e).vec(2)= 0.0;
    end
end

for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        e=i ...
            +MeshParam.Size_X*(j-MeshParam.Fine_Y_to-1) ...
            +MeshParam.coarse_ENum_Y_former+MeshParam.fine_ENum_Y...
            +MeshParam.ENum_X;
        edgevec.dual(e).vec(1)=-1.0;
        edgevec.dual(e).vec(2)= 0.0;
    end
end


%% initializing tilde_node_position


% coarce-quadrilateral zone 1: main
for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        f=MeshParam.Size_X*(j-1)+i;
        tilde_f(f).position(1)=-0.5+i;
        tilde_f(f).position(2)=-0.5+j;
    end
end

% fine-quadrilateral zone:main
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        f=i  +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        tilde_f(f).position(1)=-0.5+i;
        tilde_f(f).position(2)=-0.5+j;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localfaceIdx=1:4
            f=localfaceIdx...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_FNum_former;
        tilde_f(f).position(1)=-0.5+i-0.25+0.5*ceil((localfaceIdx-2)/2);
        tilde_f(f).position(2)=-0.5+j+0.25-0.5*mod(localfaceIdx,2);
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        f=i-MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        tilde_f(f).position(1)=-0.5+i;
        tilde_f(f).position(2)=-0.5+j;
    end
end


% coarce-quadrilateral zone 2: main
for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        f=i ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1))...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
        tilde_f(f).position(1)=-0.5+i;
        tilde_f(f).position(2)=-0.5+j;
    end
end

%% initializing UpdateNum


% coarce-quadrilateral zone 1: main
for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        f=MeshParam.Size_X*(j-1)+i;
        UpdateNum.f(f)=1;
    end
end

% fine-quadrilateral zone:main
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
        f=i  +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        UpdateNum.f(f)=1;
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localfaceIdx=1:4
            f=localfaceIdx...
                +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
                +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
                +MeshParam.coarse_FNum_former;
        UpdateNum.f(f)=MeshParam.UpdateNum_subgrid;
        end
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        f=i-MeshParam.Fine_X_to ...
            +MeshParam.Fine_X_from-1 +4*(MeshParam.Fine_X_to-MeshParam.Fine_X_from+1) ...
            +MeshParam.FacePerRow_Subgrid *(j  -MeshParam.Fine_Y_from) ...
            +MeshParam.coarse_FNum_former;
        UpdateNum.f(f)=1;
    end
end

% coarce-quadrilateral zone 2: main
for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        f=i ...
            +MeshParam.Size_X*(j  -(MeshParam.Fine_Y_to+1))...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum;
        UpdateNum.f(f)=1;
    end
end

% genrate UpdateNum.e,UpdateNum.n
% can be prosessed in parallel

for e=1:MeshNum.E
    UpdateNum.e(e)=0;
end
for kk=1:size(row_sC,1)
    f=row_sC(kk);
    e=col_sC(kk);
    %disp(e)
    %i=0;
    if UpdateNum.e(e)<UpdateNum.f(f)
        UpdateNum.e(e)=UpdateNum.f(f);
        %i=i+1;
    end
end
for kk=1:size(row_sG,1)
    e=row_sG(kk);
    n=col_sG(kk);
    if UpdateNum.n(n)<UpdateNum.e(e)
        UpdateNum.n(n)=UpdateNum.e(e);
    end
end
% end UpdateNum.e, UpdateNum.n

%% generate st faces
%disp('generating st faces')
p=1;
%Omega=1;
for f=1:MeshNum.F
    first_pIdx.f(f)=p;
    %first_Omega_for_f(f)=Omega;
    %l(p)=0;
    p=p+UpdateNum.f(f)+1;
    %Omega=Omega+UpdateNum.f(f);
end
%OmegaNum=Omega-1;
%Omega=1;
for e=1:MeshNum.E
    first_pIdx.e(e)=p;
    %first_Omega_for_e(e)=Omega;
    %l(p)=0;
    p=p+UpdateNum.e(e)+1;
    %Omega=Omega+UpdateNum.e(e);
end
%OmegaDualNum=Omega-1;
PNum=p-1;

MeshNum.P=PNum;

%% END
if DISPDEBUGGINGMESSAGE 
    disp('GenerateMesh_square_belt:ENDED')
end



end