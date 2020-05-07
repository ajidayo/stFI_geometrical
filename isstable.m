function [Stability] = isstable(cdt,Taskorder,task,Zinv_p,D,Ctrans,sC,sG,UpdateNum,edgevec,subG_bin,first_pIdx,att,MeshNum)
[kappa,~,~,MeshNum] = Constitutive(cdt,sC,sG,UpdateNum,edgevec,first_pIdx,att,MeshNum);

kappaoverZ=kappa.*Zinv_p;
D_tildeD_Zinv=[D;Ctrans * spdiags(kappaoverZ,0,MeshNum.P,MeshNum.P)];
[TMM_Explicit] ...
    = Construct_TMM_Explicit(Taskorder,task,D_tildeD_Zinv,kappaoverZ,sC,UpdateNum,subG_bin,first_pIdx,MeshNum);

EigvEpsilon=10^(-3);
eigenvalues = eigs(TMM_Explicit,20,'largestabs');
% eigenvalues = eigs(TMM_Explicit,1,'largestabs','SubspaceDimension', 50,'MaxIterations',500);

%% somehow works
IdxUnstabEigVal=find(abs(eigenvalues)>1+EigvEpsilon);
if size(IdxUnstabEigVal,1)==0
    Stability = true;
else
    Stability = false;
end

%% somehow does not work
% if isnan(Abs_LargestEigV)
%     disp('"eigs" did not converge')
%     Stability = false;
% elseif abs(Abs_LargestEigV)>1+EigvEpsilon
%     Stability = true;
% elseif abs(Abs_LargestEigV)<1+EigvEpsilon
%     Stability = false;
% else 
%     disp('error in function "isstable"')
%     pause
% end

%% end

end