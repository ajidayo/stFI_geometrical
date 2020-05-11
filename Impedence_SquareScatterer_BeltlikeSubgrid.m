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

ScatFrom_i=floor(2*ScattererMeasurements.FromXCoord);
ScatTo_i=ceil(2*ScattererMeasurements.ToXCoord);
ScatFrom_j=2*((1/2)*floor(2*ScattererMeasurements.FromYCoord)...
    -floor(MeshMeasurements.FineStartAtYCoord))...
    +floor(MeshMeasurements.FineStartAtYCoord);
ScatFrom_j=round(ScatFrom_j);
ScatTo_j=2*((1/2)*ceil(2*ScattererMeasurements.ToYCoord)...
    -floor(MeshMeasurements.FineStartAtYCoord))...
    +floor(MeshMeasurements.FineStartAtYCoord);
ScatTo_j=round(ScatTo_j);
% medium
for j=ScatFrom_j:ScatTo_j
    for i=ScatFrom_i:ScatTo_i
        f=2*MeshParam.Size_X*(j-MeshParam.Fine_Y_from)+i...
            +MeshParam.coarse_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;       
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
