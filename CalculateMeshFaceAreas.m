function [MeshFaceAreas]=CalculateMeshFaceAreas(MeshParam,cdt)
cdt_subgrid=cdt/MeshParam.UpdateNum_subgrid;
delta=MeshParam.deltaboundary;
MeshFaceAreas.Bz_next2outercorner=0.5*(1-delta+1-delta+1.0/6.0)/2.0 ...
    +(0.5+delta)*(1-delta+1.0/6.0+1-delta)/2.0 ...
    -0.5*delta*(1-delta);
MeshFaceAreas.Bz_outercorner=2*0.5*1.0*(1-delta);
MeshFaceAreas.Bz_outsideboundary=2*0.5*(1-delta+1.0/6.0+1-delta)/2.0;
MeshFaceAreas.E_outsideboundary=(1.0-MeshParam.deltaboundary)*cdt-cdt^3/12.0;
MeshFaceAreas.Bz_insideboundary_Initial  = 0.5*(0.5+delta+0.5+delta-1.0/6.0)*0.5;
MeshFaceAreas.Bz_insideboundary_HemiStep = MeshFaceAreas.Bz_insideboundary_Initial+cdt^2/12.0;
MeshFaceAreas.Bz_innercorner_Initial  = 2*0.5*(0.5+delta-1.0/6.0          )*(0.5+delta          );
MeshFaceAreas.Bz_innercorner_HemiStep = 2*0.5*(0.5+delta-1.0/6.0+cdt^2/6.0)*(0.5+delta+cdt^2/6.0);
MeshFaceAreas.E_insideboundary_smaller=0.5*(0.5+delta-1.0/6.0+0.5+delta-1.0/6.0+cdt^2/6.0)*cdt_subgrid;
MeshFaceAreas.E_insideboundary_larger =0.5*(0.5+delta        +0.5+delta        +cdt^2/6.0)*cdt_subgrid;
MeshFaceAreas.E_innercorner=0.5*( ...
    ((sqrt(10))^(-1))*(1.5+3*MeshParam.deltaboundary+1.0/6.0)...
    +((sqrt(10))^(-1))*(1.5+3*MeshParam.deltaboundary+cdt^2/2.0+1.0/6.0)...
    )*cdt_subgrid;
end