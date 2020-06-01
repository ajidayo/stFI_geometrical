function [B_mesh] = Plot_Bz_ConventionalstFI(F,MeshParam,Area)

B_mesh=zeros(MeshParam.Size_X,MeshParam.Size_Y);
area_squareoid=zeros(MeshParam.Size_X,MeshParam.Size_Y);

for j=1:MeshParam.Size_Y
    for i=1:MeshParam.Size_X

        if MeshParam.Fine_X_from < i && i < MeshParam.Fine_X_to ...
                && MeshParam.Fine_Y_from < j && j < MeshParam.Fine_Y_to
            for localjIdx=1:2
                for localiIdx=1:2
                    B_mesh(i,j)=B_mesh(i,j)+F(i,j).Bz(localiIdx,localjIdx);
                end
            end
            area_squareoid(i,j)=1.0;
        elseif ((MeshParam.Fine_X_from < i && i < MeshParam.Fine_X_to) && ...
                (j==MeshParam.Fine_Y_from || j==MeshParam.Fine_Y_to)) || ...
                ((MeshParam.Fine_Y_from < j && j < MeshParam.Fine_Y_to) && ...
                (i==MeshParam.Fine_X_from || i==MeshParam.Fine_X_to))
            for localjIdx=1:2
                for localiIdx=1:2
                    B_mesh(i,j)=B_mesh(i,j)+F(i,j).Bz(localiIdx,localjIdx);
                end
            end
            area_squareoid(i,j)=0.5+2*Area.Bz_insideboundary;
        elseif ((i==MeshParam.Fine_X_from || i == MeshParam.Fine_X_to) && ...
                (j==MeshParam.Fine_Y_from || j == MeshParam.Fine_Y_to))
            for localjIdx=1:2
                for localiIdx=1:2
                    B_mesh(i,j)=B_mesh(i,j)+F(i,j).Bz(localiIdx,localjIdx);
                end
            end
            area_squareoid(i,j)=0.25+Area.Bz_innercorner+2*Area.Bz_insideboundary;
        elseif ((i == MeshParam.Fine_X_from-1 || i == MeshParam.Fine_X_to+1)...
                && MeshParam.Fine_Y_from+1 <= j && j <= MeshParam.Fine_Y_to-1) ...
                || ((j == MeshParam.Fine_Y_from-1 || j == MeshParam.Fine_Y_to+1)...
                && MeshParam.Fine_X_from+1 <= i && i <= MeshParam.Fine_X_to-1)
            B_mesh(i,j)=F(i,j).Bz;
            area_squareoid(i,j)=Area.Bz_outsideboundary;
        elseif ((i==MeshParam.Fine_X_from || i == MeshParam.Fine_X_to) && ...
                (j==MeshParam.Fine_Y_from-1 || j == MeshParam.Fine_Y_to+1)) || ...
                ((j==MeshParam.Fine_Y_from || j == MeshParam.Fine_Y_to) && ...
                (i==MeshParam.Fine_X_from-1 || i == MeshParam.Fine_X_to+1))
            B_mesh(i,j)=F(i,j).Bz;
            area_squareoid(i,j)=Area.Bz_next2outercorner;
        elseif ((i==MeshParam.Fine_X_from-1 || i == MeshParam.Fine_X_to+1) && ...
                (j==MeshParam.Fine_Y_from-1 || j == MeshParam.Fine_Y_to+1))
            B_mesh(i,j)=F(i,j).Bz;
            area_squareoid(i,j)=Area.Bz_outercorner;
        else
            B_mesh(i,j)=F(i,j).Bz;
            area_squareoid(i,j)=1.0;
        end
    end
end

B_mesh=B_mesh./area_squareoid;

figure('name','Bz, Calculated by New Method')
xa = gca;
mesh(B_mesh.')
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
xlim([0+0.5 100+0.5])
ylim([0+0.5 100+0.5])
xticks([0+0.5 25+0.5 50+0.5 75+0.5 100+0.5])
yticks([0+0.5 25+0.5 50+0.5 75+0.5 100+0.5])
xa.FontSize = 20;
xticklabels([0 25 50 75 100])
yticklabels([0 25 50 75 100])
zlabel('B_{z}','FontSize',30)
pbaspect([1 1 1])
%color = colorbar('southoutside');
color = colorbar;
color.Label.String = 'B_{z}';

