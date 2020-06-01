function [F] = Initialize_Conventional_stFI_Explicit_w_GaussianDistribution(F,GaussParam,MeshParam,Area)

for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
        *1.0;
    end
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Fine_X_from-2
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
        *1.0;
    end
end

for j=MeshParam.Fine_Y_from-1:(MeshParam.Fine_Y_to+1)-(MeshParam.Fine_Y_from-1):MeshParam.Fine_Y_to+1
    for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_to+1)-(MeshParam.Fine_X_from-1):MeshParam.Fine_X_to+1
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_outercorner;
    end
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to-MeshParam.Fine_X_from:MeshParam.Fine_X_to
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_next2outercorner;
    end
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_outsideboundary;
    end
end
for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_to+1)-(MeshParam.Fine_X_from-1):MeshParam.Fine_X_to+1
    for j=MeshParam.Fine_Y_from-1:(MeshParam.Fine_Y_to+1)-(MeshParam.Fine_Y_from-1):MeshParam.Fine_Y_to+1
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_outercorner;
    end
    for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_next2outercorner;
    end
    for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_outsideboundary;
    end
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
        *1.0;
    end
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        x=i-0.5;
        y=j-0.5;
        F(i,j).Bz=GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
        *1.0;
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
        for localjIdx=1:2
            for localiIdx=1:2
                x=i-0.75+0.5*(localiIdx-1);
                y=j-0.75+0.5*(localjIdx-1);
                F(i,j).Bz(localiIdx,localjIdx)=...
                    GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
                    *0.25;
            end
        end
        if j==MeshParam.Fine_Y_from
            localjIdx=1;
        else % j==MeshParam.Fine_Y_to
            localjIdx=2;
        end
        for localiIdx=1:2
            x=i-0.75+0.5*(localiIdx-1);
            y=j-0.75+0.5*(localjIdx-1);
            F(i,j).Bz(localiIdx,localjIdx)=...
                GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
                *Area.Bz_insideboundary;
        end
    end
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to-MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localjIdx=1:2
            for localiIdx=1:2
                x=i-0.75+0.5*(localiIdx-1);
                y=j-0.75+0.5*(localjIdx-1);
                F(i,j).Bz(localiIdx,localjIdx)=...
                    GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
                    *0.25;
            end
        end
        if i==MeshParam.Fine_X_from
            localiIdx=1;
        else
            localiIdx=2;
        end
        for localjIdx=1:2
            x=i-0.75+0.5*(localiIdx-1);
            y=j-0.75+0.5*(localjIdx-1);
            F(i,j).Bz(localiIdx,localjIdx)=...
                GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
                *Area.Bz_insideboundary;
        end
    end
end

for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to-MeshParam.Fine_X_from:MeshParam.Fine_X_to
        if j==MeshParam.Fine_Y_from
            localjIdx=1;
        else
            localjIdx=2;
        end
        if i==MeshParam.Fine_X_from
            localiIdx=1;
        else
            localiIdx=2;
        end       
        x=i-0.75+0.5*(localiIdx-1);
        y=j-0.75+0.5*(localjIdx-1);
        F(i,j).Bz(localiIdx,localjIdx)=...
            GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_innercorner;
        if i==MeshParam.Fine_X_from
            localiIdx=2;
        else
            localiIdx=1;
        end
        x=i-0.75+0.5*(localiIdx-1);
        y=j-0.75+0.5*(localjIdx-1);
        F(i,j).Bz(localiIdx,localjIdx)=...
            GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_insideboundary;
        if j==MeshParam.Fine_Y_from
            localjIdx=2;
        else
            localjIdx=1;
        end
        x=i-0.75+0.5*(localiIdx-1);
        y=j-0.75+0.5*(localjIdx-1);
        F(i,j).Bz(localiIdx,localjIdx)=...
            GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *0.25;
        if i==MeshParam.Fine_Y_from
            localiIdx=1;
        else
            localiIdx=2;
        end
        x=i-0.75+0.5*(localiIdx-1);
        y=j-0.75+0.5*(localjIdx-1);
        F(i,j).Bz(localiIdx,localjIdx)=...
            GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
            *Area.Bz_insideboundary;        
    end
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
        for localjIdx=1:2
            for localiIdx=1:2
                x=i-0.75+0.5*(localiIdx-1);
                y=j-0.75+0.5*(localjIdx-1);
                F(i,j).Bz(localiIdx,localjIdx)=...
                    GaussParam.Ampl*exp(-((x-GaussParam.XCenter)^2+(y-GaussParam.XCenter)^2)/GaussParam.relaxfact) ...
                    *0.25;
            end
        end
    end
end

end