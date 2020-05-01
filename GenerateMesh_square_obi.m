PNumCheck=((Fine_Y_from-1-0)*2+(Fine_Y_to-Fine_Y_from+1)*3+(Size_Y-Fine_Y_to)*2)*Size_X ...
 +((Fine_Y_from-1)*2+(Fine_Y_to-Fine_Y_from+1+1)*3+(Size_Y-Fine_Y_to-1)*2)*Size_X ...
 +((Fine_Y_from-1)*2+(Fine_Y_to-Fine_Y_from+1)*3+(Size_Y-Fine_Y_to)*2)*Size_X;
%disp(PNumCheck)
OmegaNumCheck=((Fine_Y_from-1-0)*1+(Fine_Y_to-Fine_Y_from+1)*2+(Size_Y-Fine_Y_to)*1)*Size_X;
OmegaDualNumCheck=((Fine_Y_from-1)*1+(Fine_Y_to-Fine_Y_from+1+1)*2+(Size_Y-Fine_Y_to-1)*1)*Size_X ...
 +((Fine_Y_from-1)*1+(Fine_Y_to-Fine_Y_from+1)*2+(Size_Y-Fine_Y_to)*1)*Size_X;

%% Allocation

% OUTPUTS
%sC=spalloc(FNum,ENum,ENum*2);
%sG=spalloc(ENum,NNum,ENum*2);
%sC=sparse(i,j,v,FNum,ENum);
%sG=sparse(i,j,v,ENum,NNum);
sC=sparse(FNum,ENum);
sG=sparse(ENum,NNum);
denominator_node=zeros(NNum,1);
denominator_edge=zeros(ENum,1);
denominator_face=zeros(FNum,1);
edge_vector=zeros(ENum,dim);
tilde_edge_vector=zeros(ENum,dim);
first_p_for_f=zeros(FNum,1);
first_p_for_e=zeros(ENum,1);
first_Omega_for_f=zeros(FNum,1);
first_Omega_for_e=zeros(ENum,1);
%TEMPS
%x=0;y=0;
%COUNTERS
%n=0;e=0;f=0;p=0;i=0;j=0;


%% initializing Incidence Matrices sC and sG
disp('initializing sC')
% initialize sC
for f=1:Size_X*(Size_Y-1)
    e=f;
    sC(f,e)=1;
    e=f+Size_X;
    sC(f,e)=-1;
end
for i=1:Size_X
    f=Size_X*(Size_Y-1)+i;
    e=f;
    sC(f,e)=1;
    e=f-Size_X*(Size_Y-1);
    sC(f,e)=-1;
end
for i=1:Size_X-1
    for j=1:Size_Y
        f=i+(j-1)*Size_X;
        e=f+Size_X*Size_Y;
        sC(f,e)=-1;
        e=f+Size_X*Size_Y+1;
        sC(f,e)=1;
    end
end
for j=1:Size_Y
    f=j*Size_X;
    e=f+Size_X*Size_Y;
    sC(f,e)=-1;
    e=f+Size_X*Size_Y+1-Size_X;
    sC(f,e)=1;
end
%END=0

% initialize sG
disp('initializing sG')
for e=1:MESHSIZE_FACES-1
    n=e;
    sG(e,n)=-1;
    n=e+1;
    sG(e,n)=1;
end
for j=1:Size_Y-1
    e=Size_X*j;
    n=e+1;
    sG(e,n)=0;
    n=e+1-Size_X;
    sG(e,n)=1;
end
e=MESHSIZE_FACES;
n=e;
sG(e,n)=-1;
n=e+1-Size_X;
sG(e,n)=1;
for i=1:Size_X
    for j=1:Size_Y-1
        e=i+(j-1)*Size_X+Size_X*Size_Y;
        n=e-Size_X*Size_Y;
        sG(e,n)=-1;
        n=e-Size_X*Size_Y+Size_X;
        sG(e,n)=1;
    end
end
for i=1:Size_X
    e=Size_X*(Size_Y-1)+i+Size_X*Size_Y;
    n=e-Size_X*Size_Y;
    sG(e,n)=-1;
    n=i;
    sG(e,n)=1;
end
%%
[row_sC,col_sC]=find(sC);
[row_sG,col_sG]=find(sG);
%% initializing edge_vector, tilde_edge_vector
disp('initializing edge vectors and tilde edge vectors')
for e=1:Size_X*Size_Y
    edge_vector(e,1)=1.0;
    edge_vector(e,2)=0.0;
end
for e=Size_X*Size_Y+1:Size_X*Size_Y*2
    edge_vector(e,1)=0.0;
    edge_vector(e,2)=1.0;
end
for e=1:Size_X*Size_Y
    tilde_edge_vector(e,1)=0.0;
    tilde_edge_vector(e,2)=1.0;
end
for e=Size_X*Size_Y+1:Size_X*Size_Y*2
    tilde_edge_vector(e,1)=-1.0;
    tilde_edge_vector(e,2)=0.0;
end

%% initializing denominators
disp('initializing denominators')
%initialize denominator_face
for f=1:Size_X*(Fine_Y_from-1)
    denominator_face(f)=1;
end
for f=Size_X*(Fine_Y_from-1)+1:Size_X*Fine_Y_to
    %denominator_face(f)=2;
    %denominator_face(f)=1;
    denominator_face(f)=denominator_obi;
end
for f=Size_X*Fine_Y_to+1:Size_X*Size_Y
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
disp('generating st faces')
p=1;
Omega=1;
for f=1:FNum
    first_p_for_f(f)=p;
    first_Omega_for_f(f)=Omega;
    %l(p)=0;
    p=p+denominator_face(f)+1;
    Omega=Omega+denominator_face(f);
end
OmegaNum=Omega-1
Omega=1;
for e=1:ENum
    first_p_for_e(e)=p;
    first_Omega_for_e(e)=Omega;
    %l(p)=0;
    p=p+denominator_edge(e)+1;
    Omega=Omega+denominator_edge(e);
end
OmegaDualNum=Omega-1
PNum=p-1

%% Allocate PNum size matrices here
l=ones(PNum,1); %1:uncalculated variables 0:given variables
variables_initializer=zeros(PNum,1);
%variable_looper_for_f=zeros(FNum,2);
%variable_looper_for_e=zeros(ENum,2);

%% Initializing l(p)
p=1;
for f=1:FNum
    l(p)=0;
    p=p+1;
    for i=1:denominator_face(f)
        l(p)=1;
        p=p+1;
    end
end
for e=1:ENum
    l(p)=0;
    p=p+1;
    for i=1:denominator_edge(e)
        l(p)=1;
        p=p+1;
    end
end
%% END
disp('GenerateMesh:END')