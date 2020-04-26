function [impedance_inv_p] ...
    = impedance_triangular(ImpedanceParam,first_i_triangle_scatterer,sC,denominator,first_pIdx,MeshNum,MeshParam)


%% allocation
impedance_inverse_f=zeros(MeshNum.F,1);
impedance_inverse_e=zeros(MeshNum.E,1);
impedance_inv_p=zeros(MeshNum.P,1);

%%

disp('initializing impedance')

%% for f faces

for f=1:MeshNum.F
    impedance_inverse_f(f)=1.0/ImpedanceParam.freespace;
end

% medium
% triangular zone: odd numbered row
for j=MeshParam.Triangle_Y_from+2:2:MeshParam.Triangle_Y_to-3
    for i=first_i_triangle_scatterer:first_i_triangle_scatterer+(j-(MeshParam.Triangle_Y_from+2))/2.0-1
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
    end    
        i=first_i_triangle_scatterer+(j-(MeshParam.Triangle_Y_from+2))/2.0;
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
end

% triangular zone: even numbered row
for j=MeshParam.Triangle_Y_from+3:2:MeshParam.Triangle_Y_to-2
    for i=first_i_triangle_scatterer:first_i_triangle_scatterer+(j-(MeshParam.Triangle_Y_from+3))/2.0-1
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        
        f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i ...
            +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
        impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
        
    end
    i=first_i_triangle_scatterer+(j-(MeshParam.Triangle_Y_from+3))/2.0;
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-3 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-2 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
    
    f=4*MeshParam.Size_X*(j  -MeshParam.Triangle_Y_from)+4*i-1 ...
        +MeshParam.coarse_FNum_former+MeshParam.fine_FNum_former;
    impedance_inverse_f(f)=1.0/ImpedanceParam.medium;
end



%% for e faces (dummy)
for e=1:MeshNum.E
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1)
        f=row_sC_e(ff);
        impedance_inverse_e(e)=impedance_inverse_e(e)+0.5*impedance_inverse_f(f);
    end
end

%% expand to p

first_pIdx_f=first_pIdx.f;
first_pIdx_e=first_pIdx.e;
denominator_face=denominator.f;
denominator_edge=denominator.e;
for f=1:MeshNum.F
    for p = first_pIdx_f(f):first_pIdx_f(f)+denominator_face(f)
        impedance_inv_p(p)=impedance_inverse_f(f);
    end
end

for e=1:MeshNum.E
    for p = first_pIdx_e(e):first_pIdx_e(e)+denominator_edge(e)
        impedance_inv_p(p)=impedance_inverse_e(e);
    end
end

disp('impedance_triangular;ENDED')
end