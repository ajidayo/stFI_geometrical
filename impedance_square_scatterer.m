%% allocation
impedance_inverse_f=zeros(FNum,1);
impedance_inverse_e=zeros(ENum,1);
impedance_inverse_p=zeros(PNum,1);

%%

disp('initializing impedance')

%% for f faces

for f=1:FNum
    impedance_inverse_f(f)=1.0/impedance_freespace;
end

% medium
for j=first_j_square_scatterer:last_j_square_scatterer
    for i=first_i_square_scatterer:last_i_square_scatterer
        f=Size_X*(j-1)+i;
        impedance_inverse_f(f)=1.0/impedance_medium;       
    end
end


%% for e faces (dummy)
for e=1:ENum
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1)
        f=row_sC_e(ff);
        impedance_inverse_e(e)=impedance_inverse_e(e)+0.5*impedance_inverse_f(f);
    end
end

%% expand to p

for f=1:FNum
    for p = first_p_for_f(f):first_p_for_f(f)+denominator_face(f)
        impedance_inverse_p(p)=impedance_inverse_f(f);
    end
end

for e=1:ENum
    for p = first_p_for_e(e):first_p_for_e(e)+denominator_edge(e)
        impedance_inverse_p(p)=impedance_inverse_e(e);
    end
end