function [InitVal] ...
    =GaussianDistributBz(GaussParam,tilde_f,b_area,MeshNum,gauss_center)

global DISPDEBUGGINGMESSAGE
if DISPDEBUGGINGMESSAGE
    disp('GaussianDistributBz; CALLED')
end
%% gaussian distribution parameters

disp(['A sigma =',num2str(GaussParam.Ampl),', ', num2str(GaussParam.relaxfact)])
disp(['gauss_center.x,  gauss_center.y =',num2str(gauss_center.x), ', ', num2str(gauss_center.y)])


%% calculate initial values
InitBz=zeros(MeshNum.F,1);
for f=1:MeshNum.F
    x=tilde_f(f).position(1);
    y=tilde_f(f).position(2);
    InitBz(f)=GaussParam.Ampl*exp(-((x-gauss_center.x)^2+(y-gauss_center.y)^2)/GaussParam.relaxfact) ...
        *b_area(f);
end
InitVal.f=InitBz;
InitVal.e=zeros(MeshNum.E,1);
if DISPDEBUGGINGMESSAGE
    disp('GaussianDistributBz; ENDED')
end

end
