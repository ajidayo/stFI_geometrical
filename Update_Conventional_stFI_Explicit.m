function [F] = Update_Conventional_stFI_Explicit(F,Area_faces,Length_edges,Length_dualedges,MeshParam)

UpdateNum_subgrid=2;

[F]=Update_E_CoarseRegion(F,MeshParam,Area_faces,Length_edges,Length_dualedges);
[F]=Update_E_RegionBoundary_FirstTimeSec(F,MeshParam,Area_faces,Length_edges,Length_dualedges);
[F]=Update_E_FineRegion(F,MeshParam,Area_faces,Length_edges,Length_dualedges);
[F]=Update_Bz_FineRegion(F,MeshParam);
for TimeSec=2:UpdateNum_subgrid   
    [F]=Update_E_RegionBoundary_NotFirstTimeSec(F,MeshParam,Area_faces,Length_edges,Length_dualedges);
    [F]=Update_E_FineRegion(F,MeshParam,Area_faces,Length_edges,Length_dualedges);
    [F]=Update_Bz_FineRegion(F,MeshParam);
end
[F]=Update_Bz_CoarseRegion(F,MeshParam);

end

%%
function [F]=Update_E_CoarseRegion(F,MeshParam,Area_faces,Length_edges,Length_dualedges)
% update Ex
for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
        F(i,j).Ex=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ex ...
            +((cdt/(UpdateNum_Subgrid*1.0))*F(i,j).Bz ...
            - (cdt/(UpdateNum_Subgrid*1.0))*F(i,j-1).Bz)...
            );
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Fine_X_from-1
        F(i,j).Ex=F(i,j).Ex+F(i,j).Bz-F(i,j-1).Bz;
    end
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
        F(i,j).Ex=F(i,j).Ex+F(i,j).Bz-F(i,j-1).Bz;
    end
end
for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
    end
end
% update Ey
for j=1:MeshParam.Fine_Y_from-1
    for i=1:MeshParam.Size_X
    end
end
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=1:MeshParam.Fine_X_from-1
    end
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
    end
end
for j=MeshParam.Fine_Y_to+1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
    end
end
end

%% 
function [F]=Update_E_RegionBoundary_FirstTimeSec(F,MeshParam,Area_faces,Length_edges,Length_dualedges)
% update Ex

j=MeshParam.Fine_Y_from;
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
end
j=MeshParam.Fine_Y_to+1;
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
end
% update Ey
i=MeshParam.Fine_X_from;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end

end

%%
function [F]=Update_E_RegionBoundary_NotFirstTimeSec(F,MeshParam,Area_faces,Length_edges,Length_dualedges)
% update Ex
j=MeshParam.Fine_Y_from;
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
end
j=MeshParam.Fine_Y_to+1;
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
end
% update Ey
i=MeshParam.Fine_X_from;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end
end

%%
function [F]=Update_E_FineRegion(F,MeshParam,Area_faces,Length_edges,Length_dualedges)

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to
    end
end

end

%%
function [F]=Update_Bz_CoarseRegion(F,MeshParam)
for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
    end
end
for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Fine_X_from-1
    end
    for i=MeshParam.Fine_X_to+1:MeshParam.Size_X
    end
end
for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
    end
end
end

%%
function [F]=Update_Bz_RegionBoundary(F,MeshParam)
j=MeshParam.Fine_Y_from-1;
for i=MeshParam.Fine_X_from-1:MeshParam.Fine_X_to+1
end
j=MeshParam.Fine_Y_to+1;
for i=MeshParam.Fine_X_from-1:MeshParam.Fine_X_to+1
end
i=MeshParam.Fine_X_from-1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end
i=MeshParam.Fine_X_to+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
end

end

%%
function [F]=Update_Bz_FineRegion(F,MeshParam)
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    end
end
end

