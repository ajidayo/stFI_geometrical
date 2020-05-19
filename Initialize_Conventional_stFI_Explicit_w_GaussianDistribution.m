function [F] = Initialize_Conventional_stFI_Explicit_w_GaussianDistribution(F,GaussParam)

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        F(i,j).Bz=0;
    end
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Fine_X_from-2
        F(i,j).Bz=0;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_from-1)-(MeshParam.Fine_X_to+1):MeshParam.Fine_X_to+1
    F(i,j).Bz=0;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    F(i,j).Bz=0;
end

i=MeshParam.Fine_X_from-1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    F(i,j).Bz=0;
end
i=MeshParam.Fine_X_from+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    F(i,j).Bz=0;
end
j=MeshParam.Fine_Y_to+1;
for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_from-1)-(MeshParam.Fine_X_to+1):MeshParam.Fine_X_to+1
    F(i,j).Bz=0;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    F(i,j).Bz=0;
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        F(i,j).Bz=0;
    end
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        F(i,j).Bz=0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localjIdx=1:2
            for localiIdx=1:2
                F(i,j).Bz(localiIdx,localjIdx)=0;
            end
        end
    end
end

end