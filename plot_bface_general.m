function [B_mesh_Proposed,Area_squareoid]=plot_bface_general(b_f,Area_spatialfaces,tilde_f,MeshParam,MeshNum)
% title should be included in the argument


B_mesh_Proposed=zeros(MeshParam.Size_X,MeshParam.Size_Y);
Area_squareoid=zeros(MeshParam.Size_X,MeshParam.Size_Y);

for f=1:MeshNum.F
    i = ceil(tilde_f(f).position(1));
    j = ceil(tilde_f(f).position(2));
    
    B_mesh_Proposed(i,j)=B_mesh_Proposed(i,j)+b_f(f);
    Area_squareoid(i,j)=Area_squareoid(i,j)+Area_spatialfaces(f);
end

B_mesh_Proposed=B_mesh_Proposed./Area_squareoid;

%%

figure('name','Bz, Calculated by New Method')
xa = gca;
mesh(B_mesh_Proposed.')
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

end
