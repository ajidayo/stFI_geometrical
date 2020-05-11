function [Zinv_p]...
    = Impedance_SquareScatterer(ImpedanceParam,ScattererMeasurements,sC,UpdateNum,first_pIdx,MeshNum,MeshParam,MeshMeasurements)
%% allocation
impedance_inverse_f=zeros(MeshNum.F,1);
impedance_inverse_e=zeros(MeshNum.E,1);
Zinv_p=zeros(MeshNum.P,1);

%%

disp('initializing impedance')

%% for f faces

for f=1:MeshNum.F
    impedance_inverse_f(f)=1.0/ImpedanceParam.freespace;
end

ScatFrom_i = floor(ScattererMeasurements.FromXCoord);
ScatTo_i   =  ceil(ScattererMeasurements.ToXCoord);
ScatFrom_j = floor(ScattererMeasurements.FromXCoord);
ScatTo_j   =  ceil(ScattererMeasurements.ToXCoord);


% medium
for j=ScatFrom_j:ScatTo_j
    for i=ScatFrom_i:ScatTo_i
        for localfaceIdx=1:4
            f=localfaceIdx...
                +4*(i-MeshParam.Fine_X_from) ...
                +MeshParam.Fine_X_from-1 ...
                +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from)...
                +MeshParam.coarse_FNum_former;
            impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        end
    end
end


%% for e faces (dummy)
for e=1:MeshNum.E
    row_sC_e=find(sC(:,e));
    for ff=1:size(row_sC_e,1)
        f=row_sC_e(ff);
        impedance_inverse_e(e)=impedance_inverse_e(e)+0.5*impedance_inverse_f(f);
    end
end

%% expand to p

for f=1:MeshNum.F
    for p = first_pIdx.f(f):first_pIdx.f(f)+UpdateNum.f(f)
        Zinv_p(p)=impedance_inverse_f(f);
    end
end

for e=1:MeshNum.E
    for p = first_pIdx.e(e):first_pIdx.e(e)+UpdateNum.e(e)
        Zinv_p(p)=impedance_inverse_e(e);
    end
end