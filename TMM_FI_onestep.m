function [TMM_onestep] = ...
    TMM_FI_onestep(SG_tgt,Timesection_tgt,sC,kappatimesz,subG_f_bin,subG_e_bin,first_p)
% Arguement subSGIdx is unneccesary when passing f_lst_SG and e_lst_SG as
% arguements.
% When passing subG_f_bin and subG_e_bin as arguements, subSGIdx is neccesary.

first_p_for_f=first_p.f;
first_p_for_e=first_p.e;

PNum=size(kappatimesz,1);
FNum=size(sC,1);
ENum=size(sC,2);

Conserve_pp = true(PNum,1);
Map_p_next_f = logical(sparse(PNum,FNum));
Map_f_p_pres = logical(sparse(FNum,PNum));
Map_p_next_e = logical(sparse(PNum,ENum));
Map_e_p_pres = logical(sparse(ENum,PNum));

Map_f_inc_e_p_pres=logical(sparse(FNum,PNum));

% % limtoSG_f is the "equal-to-SG_tgt" pattern of subG_f_bin  
% limtoSG_f = subG_f_bin==SG_tgt;
% limtoSG_e = subG_e_bin==SG_tgt;

% f_lst_SG is the index of "equal-to-SG_tgt" elements 
f_lst_SG = find(subG_f_bin==SG_tgt);
e_lst_SG = find(subG_e_bin==SG_tgt);

for ff=1:size(f_lst_SG,1)
    f=f_lst_SG(ff);
    p_next=first_p_for_f(f)+Timesection_tgt;
    p_pres=first_p_for_f(f)+Timesection_tgt-1;
    Map_p_next_f(p_next,f)=true;
    Map_f_p_pres(f,p_pres)=true;
    Conserve_pp(p_next)=false;
end
for ee=1:size(e_lst_SG,1)
    e=e_lst_SG(ee);
    p_next=first_p_for_e(e)+Timesection_tgt;
    p_pres=first_p_for_e(e)+Timesection_tgt-1;
    Map_p_next_e(p_next,e)=true;
    Map_e_p_pres(e,p_pres)=true;
    Conserve_pp(p_next)=false;
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1)%faces incident to e
        f=row_sC_e(ff);
        p_pres=first_p_for_f(f)+Timesection_tgt-1;
        Map_f_inc_e_p_pres(f,p_pres)=true;
    end
end

Z=spdiags(kappatimesz.^(-1),0,PNum,PNum);
Zinv=spdiags(kappatimesz,0,PNum,PNum);

%%%%%%% BUG HERE; Map_f_p_pres misses fs incident to e s.t. att_bound_SG_e
%%%%%%% == true
TMM_onestep = Map_p_next_f*(Map_f_p_pres - sC*Map_p_next_e.'...
    *Z*Map_p_next_e*( Map_e_p_pres - sC.'*Map_f_inc_e_p_pres )*Zinv)...
    +Z*Map_p_next_e*( Map_e_p_pres - sC.'*Map_f_inc_e_p_pres )*Zinv;
% the minus for the dual grid is because [z] includes the sign-flipping
%%% history
% TMM_onestep = Map_p_next_f*(Map_f_p_pres ...
%     - spdiags(limtoSG_f,0,FNum,FNum)*sC*Map_p_next_e.'...
%     *Z*Map_p_next_e...
%     *( Map_e_p_pres - spdiags(limtoSG_e,0,ENum,ENum)*sC.'*Map_f_p_pres )...
%     *Zinv)...
%     +Z*Map_p_next_e...
%     *( Map_e_p_pres - spdiags(limtoSG_e,0,ENum,ENum)*sC.'*Map_f_p_pres )*...
%     Zinv;
%%% history
TMM_onestep = TMM_onestep + spdiags(Conserve_pp,0,PNum,PNum);

%[a,b,c]=find(TMM_onestep,20)

end