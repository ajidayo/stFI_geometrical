function plot_bface_general(b_f,b_area,tilde_node_position,MeshParam,MeshNum)
% title should be included in the argument


B_mesh=zeros(MeshParam.Size_X,MeshParam.Size_Y);
area_squareoid=zeros(MeshParam.Size_X,MeshParam.Size_Y);

for f=1:MeshNum.F
    i = ceil(tilde_node_position(f,1));
    j = ceil(tilde_node_position(f,2));
    
    B_mesh(i,j)=B_mesh(i,j)+b_f(f);
    area_squareoid(i,j)=area_squareoid(i,j)+b_area(f);
end

B_mesh=B_mesh./area_squareoid;

%%

% figure('name','Bz, Calculated by New Method')
% ax = gca;
% contour(Bz_obi)
% xlabel('x','FontSize',30)
% ylabel('y','FontSize',30)
% xticks([0 25 50 75 100])
% yticks([0 25 50 75 100])
% ax.FontSize = 20;

%%

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

end
