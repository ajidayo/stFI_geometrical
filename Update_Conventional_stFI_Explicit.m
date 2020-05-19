function [F] = Update_Conventional_stFI_Explicit(F,cdt,MeshParam,Area)

cdt_subgrid=cdt/MeshParam.UpdateNum_subgrid;

[F]=Update_E_CoarseRegion(F,cdt,MeshParam,Area);
[F]=Update_E_RegionBoundary_FirstTimeSec(F,cdt,cdt_subgrid,MeshParam,Area);
[F]=Update_E_FineRegion(F,cdt_subgrid,MeshParam);
[F]=Update_Bz_FineRegion(F,MeshParam);
for TimeSec=2:MeshParam.UpdateNum_subgrid
    [F]=Update_E_RegionBoundary_NotFirstTimeSec(F,cdt_subgrid,MeshParam,Area);
    [F]=Update_E_FineRegion(F,cdt_subgrid,MeshParam);
    [F]=Update_Bz_FineRegion(F,MeshParam);
end
[F]=Update_Bz_CoarseRegion(F,MeshParam);
    
end

%%
function [F]=Update_E_CoarseRegion(F,cdt,MeshParam,Area)

%% update Ex
j=1;
for i=1:MeshParam.Size_X
    F(i,j).Ex=((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/1.0)*F(i,j).Bz ...
        -(cdt/1.0)*F(i,MeshParam.Size_Y).Bz...
        );
end
for j=2:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X
        F(i,j).Ex=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ex ...
            +(cdt/1.0)*F(i,j).Bz ...
            -(cdt/1.0)*F(i,j-1).Bz...
            );
    end
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+2
    for i=1:MeshParam.Fine_X_from-2
        F(i,j).Ex=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ex ...
            +(cdt/1.0)*F(i,j).Bz ...
            -(cdt/1.0)*F(i,j-1).Bz...
            );
    end
end
j=MeshParam.Fine_Y_from-1;
for i=MeshParam.Fine_X_from-1:...
        (MeshParam.Fine_X_to+1)-(MeshParam.Fine_X_from-1):...
        MeshParam.Fine_X_to+1
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/Area.Bz_outercorner)*F(i,j  ).Bz ...
        -(cdt/1.0                )*F(i,j-1).Bz...
        );
end
for i=MeshParam.Fine_X_from:...
        MeshParam.Fine_X_to-MeshParam.Fine_X_from:...
        MeshParam.Fine_X_to
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/Area.Bz_next2outercorner)*F(i,j  ).Bz ...
        -(cdt/1.0                     )*F(i,j-1).Bz...
        );
end
for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/Area.Bz_outsideboundary)*F(i,j  ).Bz ...
        -(cdt/1.0                    )*F(i,j-1).Bz...
        );
end
j=MeshParam.Fine_Y_to+2;
for i=MeshParam.Fine_X_from-1:...
        (MeshParam.Fine_X_to+1)-(MeshParam.Fine_X_from-1):...
        MeshParam.Fine_X_to+1
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/1.0                )*F(i,j  ).Bz ...
        -(cdt/Area.Bz_outercorner)*F(i,j-1).Bz...
        );
end
for i=MeshParam.Fine_X_from:...
        MeshParam.Fine_X_to-MeshParam.Fine_X_from:...
        MeshParam.Fine_X_to
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/1.0                     )*F(i,j  ).Bz ...
        -(cdt/Area.Bz_next2outercorner)*F(i,j-1).Bz...
        );
end
for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    F(i,j).Ex=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ex ...
        +(cdt/1.0                    )*F(i,j  ).Bz ...
        -(cdt/Area.Bz_outsideboundary)*F(i,j-1).Bz...
        );
end

for i=MeshParam.Fine_X_from-1:...
        (MeshParam.Fine_X_to+1)-(MeshParam.Fine_X_from-1):...
        MeshParam.Fine_X_to+1
    j=MeshParam.Fine_Y_from;
    F(i,j).Ex=...
        (Area.E_outsideboundary/1.0)*( ...
        (1.0/Area.E_outsideboundary)*F(i,j).Ex ...
        +(cdt/Area.Bz_next2outercorner      )*F(i,j  ).Bz ...
        -(cdt/Area.Bz_outercorner)*F(i,j-1).Bz...
        );
    j=MeshParam.Fine_Y_from+1;
    F(i,j).Ex=...
        (Area.E_outsideboundary/1.0)*( ...
        (1.0/Area.E_outsideboundary)*F(i,j).Ex ...
        +(cdt/Area.Bz_outsideboundary )*F(i,j  ).Bz ...
        -(cdt/Area.Bz_next2outercorner)*F(i,j-1).Bz...
        );
    for j=MeshParam.Fine_Y_from+2:MeshParam.Fine_Y_to-2
        F(i,j).Ex=...
            (Area.E_outsideboundary/1.0)*( ...
            (1.0/Area.E_outsideboundary)*F(i,j).Ex ...
            +(cdt/Area.Bz_outsideboundary)*F(i,j  ).Bz ...
            -(cdt/Area.Bz_outsideboundary)*F(i,j-1).Bz...
            );
    end
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+2   
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X
        F(i,j).Ex=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ex ...
            + (cdt/1.0)*F(i,j).Bz ...
            - (cdt/1.0)*F(i,j-1).Bz...
            );
    end
end
for j=MeshParam.Fine_Y_to+3:MeshParam.Size_Y
    for i=1:MeshParam.Size_X
        F(i,j).Ex=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ex ...
            + (cdt/1.0)*F(i,j).Bz ...
            - (cdt/1.0)*F(i,j-1).Bz...
            );
    end
end

%% update Ey
for j=1:MeshParam.Size_Y
    i=1;
    F(i,j).Ey=((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/1.0)*F(i  ,j).Bz ...
        +(cdt/1.0)*F(MeshParam.Size_X,j).Bz...
        );
    for i=2:MeshParam.Fine_X_from-2
        F(i,j).Ey=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ey ...
            -(cdt/1.0)*F(i  ,j).Bz ...
            +(cdt/1.0)*F(i-1,j).Bz...
            );
    end
end

for j=1:MeshParam.Fine_Y_from-2
    for i=MeshParam.Fine_X_from-1:MeshParam.Fine_X_to+2
        F(i,j).Ey=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ey ...
            -(cdt/1.0)*F(i  ,j).Bz ...
            +(cdt/1.0)*F(i-1,j).Bz...
            );
    end
end
i=MeshParam.Fine_X_from-1;
for j=MeshParam.Fine_Y_from-1:...
        (MeshParam.Fine_Y_to+1)-(MeshParam.Fine_Y_from-1):...
        MeshParam.Fine_Y_to+1
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/Area.Bz_outercorner)*F(i  ,j).Bz ...
        +(cdt/1.0                )*F(i-1,j).Bz...
        );
end
for j=MeshParam.Fine_Y_from:...
        MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:...
        MeshParam.Fine_Y_to
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/Area.Bz_next2outercorner)*F(i  ,j).Bz ...
        +(cdt/1.0                     )*F(i-1,j).Bz...
        );
end
for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/Area.Bz_outsideboundary)*F(i  ,j).Bz ...
        +(cdt/1.0                    )*F(i-1,j).Bz...
        );
end
i=MeshParam.Fine_X_to+2;
for j=MeshParam.Fine_Y_from-1:...
        (MeshParam.Fine_Y_to+1)-(MeshParam.Fine_Y_from-1):...
        MeshParam.Fine_Y_to+1
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/1.0                )*F(i  ,j).Bz ...
        +(cdt/Area.Bz_outercorner)*F(i-1,j).Bz...
        );
end
for j=MeshParam.Fine_Y_from:...
        MeshParam.Fine_Y_to-MeshParam.Fine_Y_from:...
        MeshParam.Fine_Y_to
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/1.0                     )*F(i  ,j).Bz ...
        +(cdt/Area.Bz_next2outercorner)*F(i-1,j).Bz...
        );
end
for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    F(i,j).Ey=...
        ((1.0*cdt)/1.0)*( ...
        (1.0/(1.0*cdt))*F(i,j).Ey ...
        -(cdt/1.0                    )*F(i  ,j).Bz ...
        +(cdt/Area.Bz_outsideboundary)*F(i-1,j).Bz...
        );
end

for j=MeshParam.Fine_Y_from-1:...
        (MeshParam.Fine_Y_to+1)-(MeshParam.Fine_Y_from-1):...
        MeshParam.Fine_Y_to+1
    i=MeshParam.Fine_X_from;
    F(i,j).Ey=...
        (Area.E_outsideboundary/1.0)*( ...
        (1.0/Area.E_outsideboundary)*F(i,j).Ey ...
        -(cdt/Area.Bz_next2outercorner)*F(i  ,j).Bz ...
        +(cdt/Area.Bz_outercorner     )*F(i-1,j).Bz...
        );
    i=MeshParam.Fine_X_from+1;
    F(i,j).Ey=...
        (Area.E_outsideboundary/1.0)*( ...
        (1.0/Area.E_outsideboundary)*F(i,j).Ey ...
        -(cdt/Area.Bz_outsideboundary )*F(i  ,j).Bz ...
        +(cdt/Area.Bz_next2outercorner)*F(i-1,j).Bz...
        );
    for i=MeshParam.Fine_X_from+2:MeshParam.Fine_X_to-2
        F(i,j).Ey=...
            (Area.E_outsideboundary/1.0)*( ...
            (1.0/Area.E_outsideboundary)*F(i,j).Ey ...
            -(cdt/Area.Bz_outsideboundary)*F(i  ,j).Bz ...
            +(cdt/Area.Bz_outsideboundary)*F(i-1,j).Bz...
            );
    end
end

for i=MeshParam.Fine_X_from-1:MeshParam.Fine_X_to+2   
    for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y
        F(i,j).Ey=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ey ...
            - (cdt/1.0)*F(i  ,j).Bz ...
            + (cdt/1.0)*F(i-1,j).Bz...
            );
    end
end
for i=MeshParam.Fine_X_to+3:MeshParam.Size_X
    for j=1:MeshParam.Size_Y
        F(i,j).Ey=((1.0*cdt)/1.0)*( ...
            (1.0/(1.0*cdt))*F(i,j).Ey ...
            - (cdt/1.0)*F(i  ,j).Bz ...
            + (cdt/1.0)*F(i-1,j).Bz...
            );
    end
end
%% end 
end

%% 
function [F]=Update_E_RegionBoundary_FirstTimeSec(F,cdt,cdt_subgrid,MeshParam,Area)

%% update Ex
j=MeshParam.Fine_Y_from;

i=MeshParam.Fine_X_from;
localjIdx=1;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i,j-1).Bz...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i,j-1).Bz...
    );
localjIdx=2;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_innercorner    )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );

i=MeshParam.Fine_X_to;
localjIdx=1;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i,j-1).Bz...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i,j-1).Bz...
    );
localjIdx=2;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_innercorner    )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );

for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    for localiIdx=1:2
        localjIdx=1;
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx) ...
            - (cdt        /Area.Bz_outsideboundary)*F(i,j-1).Bz...
            );
        localjIdx=2;
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
            - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
            );
    end
end

j=MeshParam.Fine_Y_to+1;

i=MeshParam.Fine_X_from;
localiIdx=1;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i,j  ).Bz ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i,j-1).Bz(localiIdx,2        )...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i,j  ).Bz ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j-1).Bz(localiIdx,2        ) ...
);

i=MeshParam.Fine_X_to;
localiIdx=1;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i,j  ).Bz ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j-1).Bz(localiIdx,2        ) ...
);
localiIdx=2;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i,j  ).Bz ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i,j-1).Bz(localiIdx,2        )...
    );

for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    for localiIdx=1:2
        F(i,j).Ex(localiIdx,1)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1) ...
            + (cdt        /Area.Bz_outsideboundary)*F(i,j  ).Bz...
            - (cdt_subgrid/Area.Bz_insideboundary )*F(i,j-1).Bz(localiIdx,2) ...
            );
    end
end

%% update Ey
i=MeshParam.Fine_X_from;

j=MeshParam.Fine_Y_from;
localiIdx=1;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i-1,j).Bz...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i-1,j).Bz...
    );
localiIdx=2;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner    )*F(i-1,j).Bz(localiIdx-1,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );

j=MeshParam.Fine_Y_to;
localiIdx=1;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i-1,j).Bz...
    );
localjIdx=2;
F(i,j).Exy(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    + (cdt        /Area.Bz_next2outercorner)*F(i-1,j).Bz...
    );
localiIdx=2;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner    )*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for localjIdx=1:2
        localiIdx=1;
        F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
            - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx,localjIdx) ...
            + (cdt        /Area.Bz_outsideboundary)*F(i-1,j).Bz...
            );
        localiIdx=2;
        F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
            - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
            );
    end
end

i=MeshParam.Fine_X_to+1;

j=MeshParam.Fine_Y_from;
localjIdx=1;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i  ,j).Bz ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i-1,j).Bz(2,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i  ,j).Bz ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i-1,j).Bz(2,localjIdx) ...
);

j=MeshParam.Fine_Y_to;
localjIdx=1;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i  ,j).Bz ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i-1,j).Bz(2,localjIdx) ...
);
localjIdx=2;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    - (cdt        /Area.Bz_next2outercorner)*F(i  ,j).Bz ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i-1,j).Bz(2,localjIdx)...
    );

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for localjIdx=1:2
        F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
            - (cdt        /Area.Bz_outsideboundary)*F(i  ,j).Bz...
            + (cdt_subgrid/Area.Bz_insideboundary )*F(i-1,j).Bz(2,localjIdx) ...
            );
    end
end
end

%%
function [F]=Update_E_RegionBoundary_NotFirstTimeSec(F,cdt_subgrid,MeshParam,Area)
%% update Ex
j=MeshParam.Fine_Y_from;
i=MeshParam.Fine_X_from;
localjIdx=1;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    );
localjIdx=2;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_innercorner    )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );

i=MeshParam.Fine_X_to;
localjIdx=1;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i,j  ).Bz(localiIdx,localjIdx) ...
    );
localjIdx=2;
localiIdx=1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_innercorner    )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );

for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    for localiIdx=1:2
        localjIdx=1;
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary )*F(i,j  ).Bz(localiIdx,localjIdx) ...
            );
        localjIdx=2;
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
            - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
            );
    end
end

i=MeshParam.Fine_X_to;
for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    localjIdx=1;
    localiIdx =1;
    F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
        (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
        + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
        - (cdt_subgrid/0.25                  )*F(i,j-1).Bz(localiIdx,localjIdx+1)...
        );
    localiIdx =2;
    F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_larger/0.5)*( ...
        (0.5/Area.E_insideboundary_larger)*F(i,j).Ex(localiIdx,localjIdx) ...
        + (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
        - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j-1).Bz(localiIdx,localjIdx+1)...
        );
    
    localjIdx=2;
    localiIdx =1;
    F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
        (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
        + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
        - (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
        );
    localiIdx =2;
    F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
        (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
        + (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
        - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
        );
end

j=MeshParam.Fine_Y_to;
for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    localjIdx=1;
    for localiIdx =1:2
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
            - (cdt_subgrid/0.25                  )*F(i,j-1).Bz(localiIdx,localjIdx+1)...
            );    
    end
    localjIdx=2;
    for localiIdx =1:2                
        F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
            - (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
            );
    end
end

i=MeshParam.Fine_X_to;
j=MeshParam.Fine_Y_to;
localjIdx = 1;
localiIdx = 1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/0.25                  )*F(i,j-1).Bz(localiIdx,localjIdx+1)...
    );
localiIdx = 2;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_larger/0.5)*( ...
    (0.5/Area.E_insideboundary_larger)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j-1).Bz(localiIdx,localjIdx+1)...
    );
localjIdx=2;
localiIdx =1;
F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/0.25                  )*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );
localiIdx =2;
F(i,j).Ex(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ex(localiIdx,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner   )*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
    - (cdt_subgrid/Area.Bz_insideboundary)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
    );

j=MeshParam.Fine_Y_to+1;
i=MeshParam.Fine_X_from;
localiIdx=1;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i,j-1).Bz(localiIdx,2        )...
    );
localiIdx=2;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j-1).Bz(localiIdx,2        ) ...
);

i=MeshParam.Fine_X_to;
localiIdx=1;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i,j-1).Bz(localiIdx,2        ) ...
);
localiIdx=2;
F(i,j).Ex(localiIdx,1        )=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1        ) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i,j-1).Bz(localiIdx,2        )...
    );

for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
    for localiIdx=1:2
        F(i,j).Ex(localiIdx,1)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,1) ...
            - (cdt_subgrid/Area.Bz_insideboundary )*F(i,j-1).Bz(localiIdx,2) ...
            );
    end
end

%% update Ey
i=MeshParam.Fine_X_from;

j=MeshParam.Fine_Y_from;
localiIdx=1;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    );
localiIdx=2;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner    )*F(i-1,j).Bz(localiIdx-1,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );

j=MeshParam.Fine_Y_to;
localiIdx=1;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary  )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    );
localjIdx=2;
F(i,j).Exy(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_innercorner     )*F(i  ,j).Bz(localiIdx,localjIdx) ...
    );
localiIdx=2;
localjIdx=1;
F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
    (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(localiIdx,localjIdx)=(Area.E_insideboundary_smaller/0.5)*( ...
    (0.5/Area.E_insideboundary_smaller)*F(i,j).Ey(localiIdx,localjIdx) ...
    - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner    )*F(i  ,j).Bz(localiIdx-1,localjIdx)...
    );

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for localjIdx=1:2
        localiIdx=1;
        F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
            - (cdt_subgrid/Area.Bz_insideboundary )*F(i  ,j).Bz(localiIdx,localjIdx) ...
            );
        localiIdx=2;
        F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
            (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
            - (cdt_subgrid/0.25                  )*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
            );
    end
end

i=MeshParam.Fine_X_to+1;

j=MeshParam.Fine_Y_from;
localjIdx=1;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i-1,j).Bz(2,localjIdx)...
    );
localjIdx=2;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i-1,j).Bz(2,localjIdx) ...
);

j=MeshParam.Fine_Y_to;
localjIdx=1;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    + (cdt_subgrid/Area.Bz_insideboundary  )*F(i-1,j).Bz(2,localjIdx) ...
);
localjIdx=2;
F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
    (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
    + (cdt_subgrid/Area.Bz_innercorner     )*F(i-1,j).Bz(2,localjIdx)...
    );

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for localjIdx=1:2
        F(i,j).Ey(1,localjIdx)=((0.5*cdt_subgrid)/0.75)*( ...
            (0.75/(0.5*cdt_subgrid))*F(i,j).Ey(1,localjIdx) ...
            + (cdt_subgrid/Area.Bz_insideboundary )*F(i-1,j).Bz(2,localjIdx) ...
            );
    end
end
end

%%
function [F]=Update_E_FineRegion(F,cdt_subgrid,MeshParam)

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
        localjIdx=1;
        for localiIdx =1:2
            F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
                (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
                + (cdt_subgrid/0.25)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
                - (cdt_subgrid/0.25)*F(i,j-1).Bz(localiIdx,localjIdx+1)...
                );
        end
        localjIdx=2;
        for localiIdx =1:2
            F(i,j).Ex(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
                (0.5/(0.5*cdt_subgrid))*F(i,j).Ex(localiIdx,localjIdx) ...
                + (cdt_subgrid/0.25)*F(i,j  ).Bz(localiIdx,localjIdx  ) ...
                - (cdt_subgrid/0.25)*F(i,j  ).Bz(localiIdx,localjIdx-1)...
                );
        end     
    end
end

for j=MeshParam.Fine_Y_from+1:MeshParam.Fine_Y_to-1
    for i=MeshParam.Fine_X_from+1:MeshParam.Fine_X_to-1
        localiIdx=1;
        for localjIdx =1:2
            F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
                (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
                - (cdt_subgrid/0.25)*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
                + (cdt_subgrid/0.25)*F(i-1,j).Bz(localiIdx+1,localjIdx)...
                );
        end
        localiIdx=2;
        for localjIdx =1:2
            F(i,j).Ey(localiIdx,localjIdx)=((0.5*cdt_subgrid)/0.5)*( ...
                (0.5/(0.5*cdt_subgrid))*F(i,j).Ey(localiIdx,localjIdx) ...
                - (cdt_subgrid/0.25)*F(i  ,j).Bz(localiIdx  ,localjIdx) ...
                + (cdt_subgrid/0.25)*F(i  ,j).Bz(localiIdx-1,localjIdx)...
                );
        end
    end
end

end


%%
function [F]=Update_Bz_CoarseRegion(F,MeshParam)
for j=1:MeshParam.Fine_Y_from-2
    for i=1:MeshParam.Size_X-1
        F(i,j).Bz=F(i,j).Bz ...
            +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
            -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
    end
    i=MeshParam.Size_X;
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(1  ,j  ).Ey;
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=1:MeshParam.Fine_X_from-2
        F(i,j).Bz=F(i,j).Bz ...
            +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
            -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
    end
end

j=MeshParam.Fine_Y_from-1;
for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_from-1)-(MeshParam.Fine_X_to+1):MeshParam.Fine_X_to+1
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex(1,1)-F(i  ,j+1).Ex(2,1) ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
end

i=MeshParam.Fine_X_from-1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey(1,1)+F(i+1,j  ).Ey(1,2);
end
i=MeshParam.Fine_X_from+1;
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey(2,1)-F(i  ,j  ).Ey(2,2)+F(i+1,j  ).Ey;
end
j=MeshParam.Fine_Y_to+1;
for i=MeshParam.Fine_X_from-1:(MeshParam.Fine_X_from-1)-(MeshParam.Fine_X_to+1):MeshParam.Fine_X_to+1
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
end
for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex(1,1)+F(i  ,j  ).Ex(2,1)-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
end

for j=MeshParam.Fine_Y_from-1:MeshParam.Fine_Y_to+1
    for i=MeshParam.Fine_X_to+2:MeshParam.Size_X-1
        F(i,j).Bz=F(i,j).Bz ...
            +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
            -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
    end
    i=MeshParam.Size_X;
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(1  ,j  ).Ey;
end

for j=MeshParam.Fine_Y_to+2:MeshParam.Size_Y-1
    for i=1:MeshParam.Size_X-1
        F(i,j).Bz=F(i,j).Bz ...
            +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
            -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
    end
    i=1:MeshParam.Size_X;
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,j+1).Ex ...
        -F(i  ,j  ).Ey+F(1  ,j  ).Ey;
end
j=MeshParam.Size_Y;
for i=1:MeshParam.Size_X-1
    F(i,j).Bz=F(i,j).Bz ...
        +F(i  ,j  ).Ex-F(i  ,1  ).Ex ...
        -F(i  ,j  ).Ey+F(i+1,j  ).Ey;
end
i=MeshParam.Size_X;
F(i,j).Bz=F(i,j).Bz ...
    +F(i  ,j  ).Ex-F(i  ,1  ).Ex ...
    -F(i  ,j  ).Ey+F(1  ,j  ).Ey;

end
%%

function [F]=Update_Bz_FineRegion(F,MeshParam)
for j=MeshParam.Fine_Y_from:MeshParam.Fine_Y_to
    for i=MeshParam.Fine_X_from:MeshParam.Fine_X_to
        for localjIdx=1:2
            for localiIdx=1:2
                F(i,j).Bz(localiIdx,localjIdx)=F(i,j).Bz(localiIdx,localjIdx) ...
                    +F(i  ,j            ).Ex(localiIdx,localjIdx                    )...
                    -F(i  ,j+localjIdx-1).Ex(localiIdx,localjIdx+sign(1.5-localjIdx)) ...
                    -F(i            ,j  ).Ey(localiIdx                    ,localjIdx)...
                    +F(i+localiIdx-1,j  ).Ey(localiIdx+sign(1.5-localiIdx),localjIdx);
            end
        end
    end
end
end

