% for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
%     for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to-MeshParam.Fine_X_from:MeshParam.Fine_X_to
%         for localfaceIdx=1:4
%             f=localfaceIdx...
%                 +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
%                 +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
%                 +MeshParam.coarse_FNum_former;
%             f
%             Area_spatialfaces(f)
%         end
%     end
% end
% for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1-MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to-1
%     for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1-MeshParam.Fine_X_from-1:MeshParam.Fine_X_to-1
%         for localfaceIdx=1:4
%             f=localfaceIdx...
%                 +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
%                 +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
%                 +MeshParam.coarse_FNum_former;
%             kappa(first_pIdx.f(f))
%         end
%     end
% end
% for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
%     for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to-MeshParam.Fine_X_from:MeshParam.Fine_X_to
%         for localfaceIdx=1:4
%             f=localfaceIdx...
%                 +4*(i-MeshParam.Fine_X_from)+MeshParam.Fine_X_from-1 ...
%                 +MeshParam.FacePerRow_Subgrid*(j-MeshParam.Fine_Y_from) ...
%                 +MeshParam.coarse_FNum_former;
%             f
%             kappa(first_pIdx.f(f)+1)
%         end
%     end
% end
Epsilon_SymmetryTest=10^(-14);
Bz_Proposed_diagonallysymmetricity_error=zeros(MeshParam.Size_X,MeshParam.Size_Y);
for j=1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        Temp=B_mesh_Proposed(i,j)-B_mesh_Proposed(j,i);
        UnSymFlag=Temp>Epsilon_SymmetryTest;
        if UnSymFlag
            disp(['unsymmetric found: i=',num2str(i),', j=',num2str(j)])
        end
        Bz_Proposed_diagonallysymmetricity_error(i,j)=UnSymFlag;
    end
end
UnSymFlag=any(any(find(Bz_Proposed_diagonallysymmetricity_error>EPSILON)));



