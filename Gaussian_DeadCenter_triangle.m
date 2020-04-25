function [InitVal] ...
    =Gaussian_DeadCenter_triangle(A,sigma,tilde_node_position,b_area,MeshNum,MeshParam)
%% gaussian distribution parameters


% dead center

gauss_center.x=MeshParam.Size_X/2.0;
gauss_center.y=0.5*(MeshParam.Fine_Y_from-1 ...
    +(MeshParam.Fine_Y_to-MeshParam.Fine_Y_from+1)/2.0...
    +(MeshParam.Size_Y-MeshParam.Fine_Y_to));

disp('Initial conditions: Gaussian Distribution of Bz, centered at the Dead center of the mesh')
disp('A sigma =')
disp([A sigma])
disp('x_o y_o =')
disp([gauss_center.x gauss_center.y])


%% calculate initial values
InitBz=zeros(MeshNum.F,1);
for f=1:MeshNum.F
    x=tilde_node_position(f,1);
    y=tilde_node_position(f,2);
    InitBz(f)=A*exp(-((x-gauss_center.x)^2+(y-gauss_center.y)^2)/sigma) ...
        *b_area(f);
end
InitVal.f=InitBz;
InitVal.e=sparse(MeshNum.E,1);
disp('Gaussian_DeadCenter_triangle; ENDED')
end
