function PlotMagneticFluxDensity3D(DoFs_FacesThenEdges,FaceArea,Num_of_Elem,SpElemProperties,MeshMeasurements)
global SpDIM EPSILON

dx = MeshMeasurements.dx;
dy = MeshMeasurements.dy;
dz = MeshMeasurements.dz;
XSize = MeshMeasurements.XCoord/MeshMeasurements.dx;
YSize = MeshMeasurements.YCoord/MeshMeasurements.dy;
ZSize = MeshMeasurements.ZCoord/MeshMeasurements.dz;

% Disp_dx = dx;
% Disp_dy = dy;
% Disp_dz = dz;

B_SquareoidMesh = zeros(XSize+1,YSize+1,ZSize+1,SpDIM);
Area_Squareoid  = ones(XSize+1,YSize+1,ZSize+1,SpDIM);
for SpPIdx = 1:Num_of_Elem.SpP
    xIdx = dx*round(SpElemProperties.SpP.Position.x(SpPIdx)/dx);
    yIdx = dy*round(SpElemProperties.SpP.Position.y(SpPIdx)/dy);
    zIdx = dz*round(SpElemProperties.SpP.Position.z(SpPIdx)/dz);
    if abs(xIdx - round(xIdx)) < EPSILON
        dimIdx = 1
        xIdx = round(xIdx)+1;
        yIdx =  ceil(yIdx+EPSILON);
        zIdx =  ceil(zIdx+EPSILON);
    elseif abs(yIdx - round(yIdx)) < EPSILON
        dimIdx = 2
        xIdx =  ceil(xIdx+EPSILON);
        yIdx = round(yIdx)+1;
        zIdx =  ceil(zIdx+EPSILON);
    elseif abs(zIdx - round(zIdx)) < EPSILON
        dimIdx = 3
        xIdx =  ceil(xIdx+EPSILON);
        yIdx =  ceil(yIdx+EPSILON);
        zIdx = round(zIdx)+1;
    else
        disp("hoge")
        continue;
    end
    B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) ...
        = B_SquareoidMesh(xIdx,yIdx,zIdx,dimIdx) + DoFs_FacesThenEdges(SpPIdx);
    Area_Squareoid(xIdx,yIdx,zIdx,dimIdx)  ...
        =  Area_Squareoid(xIdx,yIdx,zIdx,dimIdx) +       FaceArea.Prim(SpPIdx);
end

B_SquareoidMesh=B_SquareoidMesh./Area_Squareoid;

PlottingVectorIdx = 0;
for zIdx = 1:ZSize 
    for yIdx = 1:YSize
        for xIdx = 1:XSize
            PlottingVectorIdx = PlottingVectorIdx+1;
            VecPosX(PlottingVectorIdx) = (xIdx-0.5)*dx;
            VecPosY(PlottingVectorIdx) = (yIdx-0.5)*dy;
            VecPosZ(PlottingVectorIdx) = (zIdx-0.5)*dz;
            B_Vol_X(PlottingVectorIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,1)...
                +B_SquareoidMesh(xIdx+1,yIdx  ,zIdx  ,1));
            B_Vol_Y(PlottingVectorIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,2)...
                +B_SquareoidMesh(xIdx  ,yIdx+1,zIdx  ,2));
            B_Vol_Z(PlottingVectorIdx) = ...
                (1/6)*...
                (B_SquareoidMesh(xIdx  ,yIdx  ,zIdx  ,3)...
                +B_SquareoidMesh(xIdx  ,yIdx  ,zIdx+1,3));
        end
    end
end

%%

figure('name','B, Calculated by New Method')
xa = gca;
scale = 100;
quiver3(VecPosX,VecPosY,VecPosZ,B_Vol_X,B_Vol_Y,B_Vol_Z,scale)
xlabel('x','FontSize',30)
ylabel('y','FontSize',30)
zlabel('z','FontSize',30)
xlim([0 XSize*dx])
ylim([0 YSize*dy])
ylim([0 ZSize*dz])
xticks([0 dx*XSize/4 dx*XSize/2 dx*XSize*3/4 dx*XSize])
yticks([0 dy*YSize/4 dy*YSize/2 dy*YSize*3/4 dy*YSize])
zticks([0 dz*ZSize/4 dz*ZSize/2 dz*ZSize*3/4 dz*ZSize])
xa.FontSize = 20;
xticks([0 dx*XSize/4 dx*XSize/2 dx*XSize*3/4 dx*XSize])
yticks([0 dy*YSize/4 dy*YSize/2 dy*YSize*3/4 dy*YSize])
zticks([0 dz*ZSize/4 dz*ZSize/2 dz*ZSize*3/4 dz*ZSize])
pbaspect([1 1 1])
%color = colorbar('southoutside');
color = colorbar;
color.Label.String = 'B';


end