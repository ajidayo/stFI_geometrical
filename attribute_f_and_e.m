function [att] ...
    =attribute_f_and_e(sC,sG,UpdateNum,MeshNum)

global EPSILON

% allocate attribute vector
att.e.boundaryedge=logical(sparse(MeshNum.E,1));
att.f.inc_bounde=logical(sparse(MeshNum.F,1));
att.e.adj_bounde=logical(sparse(MeshNum.E,1));
att.e.bound_of_SG=logical(sparse(MeshNum.E,1));

% make adjacent matrix
AdjE=logical(sG*sG.');

%disp('checkpoint alpha')

for e_target=1:MeshNum.E    
    %%%%%% CAPSULIZE FROM %%%%%%
    sC_e_target=sC(:,e_target);
    row_sC_e_target=find(sC_e_target);
    % check boundary-attribute
    check_bound=0.0;
    for ff=1:size(row_sC_e_target,1)
        f=row_sC_e_target(ff);
        check_bound=check_bound+sC(f,e_target)*UpdateNum.f(f);
    end
    if abs(check_bound)<EPSILON % not a boundary
         continue;
    end
    %%%%%% CAPSULIZE TO %%%%%%    
    
    %disp(e_target)
    
    % if e_target is a dt-varying boundary
    att.e.boundaryedge(e_target)=true;
    % for all f incident to e_target
    for ff=1:size(row_sC_e_target,1)
        f=row_sC_e_target(ff);
        att.f.inc_bounde(f)=true;
    end
end % for e_target

att_bound_e=att.e.boundaryedge;
allIdx_bound_e=find(att_bound_e==true);

for ee=1:size(allIdx_bound_e,1)
    e_bound=allIdx_bound_e(ee);
    AdjE_e_bound=AdjE(e_bound,:);
    col_AdjE_e_bound=find(AdjE_e_bound);
    for eee=1:size(col_AdjE_e_bound,2)
        e_target=col_AdjE_e_bound(eee);
        if att.e.boundaryedge(e_target)==false
            att.e.adj_bounde(e_target)=true;
        end
    end
end % for e_target

clearvars AdjE

%disp('checkpoint bravo')

att_inc_bound_f=att.f.inc_bounde;
allIdx_stFI.f=find(att_inc_bound_f==true);
allIdx_stFI_f=allIdx_stFI.f;

for ff=1:size(allIdx_stFI_f,1)
    %disp(ff)
    f_target=allIdx_stFI_f(ff);
    sC_f_target=sC(f_target,:);
    col_sC_f_target=find(sC_f_target);
    for ee=1:size(col_sC_f_target,2)
        e=col_sC_f_target(ee);
        if att.e.boundaryedge(e)==false && att.e.adj_bounde(e)==false
            %disp(e)
            att.e.bound_of_SG(e)=true;
        end
    end
end

end