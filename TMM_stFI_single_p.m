function [TMM_onestep] = ...
    TMM_stFI_single_p(p_tgt, omega_tgt, D_tildeD_Zinv)
PNum=size(D_tildeD_Zinv,2);

TMM_onestep=speye(PNum);

TMM_onestep(p_tgt,:) = ...
    -(1.0/D_tildeD_Zinv(omega_tgt,p_tgt))*D_tildeD_Zinv(omega_tgt,:);
TMM_onestep(p_tgt,p_tgt)=0;

end