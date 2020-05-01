function [InitVal] ...
    =GaussianDistributBz(GaussParam,tilde_node_position,b_area,MeshNum,gauss_center)
disp('GaussianDistributBz; CALLED')
%% gaussian distribution parameters

disp('A sigma =')
disp([GaussParam.Ampl GaussParam.relaxfact])
disp('x_o y_o =')
disp([gauss_center.x gauss_center.y])


%% calculate initial values
InitBz=zeros(MeshNum.F,1);
for f=1:MeshNum.F
    x=tilde_node_position(f,1);
    y=tilde_node_position(f,2);
    InitBz(f)=GaussParam.Ampl*exp(-((x-gauss_center.x)^2+(y-gauss_center.y)^2)/GaussParam.relaxfact) ...
        *b_area(f);
end
InitVal.f=InitBz;
InitVal.e=zeros(MeshNum.E,1);
disp('GaussianDistributBz; ENDED')
end
