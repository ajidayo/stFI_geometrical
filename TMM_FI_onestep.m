function [TMM_onestep] = ...
    TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappaoverZ,subG_bin,first_pIdx,MeshNum)
% Arguement subSGIdx is unneccesary when passing f_lst_SG and e_lst_SG as
% arguements.
% When passing subG_bin.f and subG_bin.e as arguements, subSGIdx is neccesary.

Conserve_pp = true(MeshNum.P,1);
Map_p_next_f = logical(sparse(MeshNum.P,MeshNum.F));
Map_f_p_pres = logical(sparse(MeshNum.F,MeshNum.P));
Map_p_next_e = logical(sparse(MeshNum.P,MeshNum.E));
Map_e_p_pres = logical(sparse(MeshNum.E,MeshNum.P));

Map_f_inc_e_p_pres=logical(sparse(MeshNum.F,MeshNum.P));

% f_lst_SG is the index of "equal-to-SG_tgt" elements 
f_lst_SG = find(subG_bin.f==SG_tgt);
e_lst_SG = find(subG_bin.e==SG_tgt);

for ff=1:size(f_lst_SG,1)
    f=f_lst_SG(ff);
    p_next=first_pIdx.f(f)+Timesection_tgt;
    p_pres=first_pIdx.f(f)+Timesection_tgt-1;
    Map_p_next_f(p_next,f)=true;
    Map_f_p_pres(f,p_pres)=true;
    Conserve_pp(p_next)=false;
end
for ee=1:size(e_lst_SG,1)
    e=e_lst_SG(ee);
    p_next=first_pIdx.e(e)+Timesection_tgt;
    p_pres=first_pIdx.e(e)+Timesection_tgt-1;
    Map_p_next_e(p_next,e)=true;
    Map_e_p_pres(e,p_pres)=true;
    Conserve_pp(p_next)=false;
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1)%faces incident to e
        f=row_sC_e(ff);
        p_pres=first_pIdx.f(f)+Timesection_tgt-1;
        Map_f_inc_e_p_pres(f,p_pres)=true;
    end
end

z=spdiags(kappaoverZ.^(-1),0,MeshNum.P,MeshNum.P);
zinv=spdiags(kappaoverZ,0,MeshNum.P,MeshNum.P);

%%%%%%% BUG HERE; Map_f_p_pres misses fs incident to e s.t. att_bound_SG_e
%%%%%%% == true
TMM_onestep = Map_p_next_f*(Map_f_p_pres - sC*Map_p_next_e.'...
    *z*Map_p_next_e*( Map_e_p_pres - sC.'*Map_f_inc_e_p_pres )*zinv)...
    +z*Map_p_next_e*( Map_e_p_pres - sC.'*Map_f_inc_e_p_pres )*zinv;
% the minus for the dual grid is because [z] includes the sign-flipping
%%% history
% TMM_onestep = Map_p_next_f*(Map_f_p_pres ...
%     - spdiags(limtoSG_f,0,MeshNum.F,MeshNum.F)*sC*Map_p_next_e.'...
%     *Z*Map_p_next_e...
%     *( Map_e_p_pres - spdiags(limtoSG_e,0,MeshNum.E,MeshNum.E)*sC.'*Map_f_p_pres )...
%     *Zinv)...
%     +Z*Map_p_next_e...
%     *( Map_e_p_pres - spdiags(limtoSG_e,0,MeshNum.E,MeshNum.E)*sC.'*Map_f_p_pres )*...
%     Zinv;
%%% history
TMM_onestep = TMM_onestep + spdiags(Conserve_pp,0,MeshNum.P,MeshNum.P);

%[a,b,c]=find(TMM_onestep,20)

end