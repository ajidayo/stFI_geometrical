function [kappa,b_area,MeshNum] ...
    =Constitutive(cdt,sC,sG,denominator,edgevec,first_pIdx,att,MeshNum)
disp('Constitutive: CALLED')

%global EPSILON
global DIM

exception=0;% dummy
%% Parameters
%cdt=0.3;
%% Allocation
%endpoints_tilde_e=zeros(2,1);

DirVec_NodeDelta=zeros(DIM,1); %(xy)

first_stn_n=zeros(MeshNum.N,1);
kappa=zeros(MeshNum.P,1);
LorenzInvariant_dS=sparse(MeshNum.P,1);

b_area=zeros(MeshNum.F,1);
%% generate st nodes

stn=1;
for n=1:MeshNum.N
    first_stn_n(n)=stn;
    stn=stn+denominator.n(n)+1;
end
stNNum=stn-1;
MeshNum.stn=stNNum;

Dummy=cell(MeshNum.stn,1);
for stn=1:MeshNum.stn
    Dummy{stn}=sparse(DIM,1);
end
stnInfo=struct('NodeDeltaVec',Dummy);
clearvars Dummy


%% calculate kappa for boundary edges
disp('calculating kappa for boundary edges')
for e_tar=1:MeshNum.E
    if att.bound_e(e_tar)==false
         continue;
    end

    row_sC_e_tar=find(sC(:,e_tar));
  
    % check exception-attribute
    if exception==1 %exception (corner)
    else %non-exceptional dt-boundary edge
        denominator_e_tar=denominator.e(e_tar);
        p_tar=first_pIdx.e(e_tar);        
        endpoints_e = endpoint(e_tar,sG);

        e_connec_to_ep1=e_connec2ntar_shareincfwith_e_tar(endpoints_e(1),e_tar,sC,sG,edgevec);
        e_connec_to_ep2=e_connec2ntar_shareincfwith_e_tar(endpoints_e(2),e_tar,sC,sG,edgevec);
        DirVec_e_tar=DirectionVec(e_tar,edgevec);

        % correct directions of vec_e1c, vec_e2c: the direction is the same
        % as the direction of edgevec.dual.(the direction of g, or tilde e)
        for e_con=1:2
            e_connec_to_ep1(e_con).DirVec=e_connec_to_ep1(e_con).DirVec...
                *sign(DirVec_e_tar.'*e_connec_to_ep1(e_con).DirVec);
            e_connec_to_ep2(e_con).DirVec=e_connec_to_ep2(e_con).DirVec...
                *sign(DirVec_e_tar.'*e_connec_to_ep2(e_con).DirVec);
        end
        
        % check that the two edges are parallel: if parallel then the direction of the two edges is the shift direction for boundary nodes

        Parallel=ParaCheck_e_cons(e_connec_to_ep1,e_connec_to_ep2);
        if Parallel.is==true
            DirVec_NodeDelta=e_connec_to_ep1(Parallel.e_con_ep1_Idx).DirVec;
        else
            disp(['error; Cannot define DirVec_NodeDelta for e_tar = ', num2str(e_tar)])
        end  
        costheta=DirVec_NodeDelta.'*DirVec_e_tar;
        sintheta=sqrt(1-costheta^2);
        
        % find endpoints of tilde_e from sC
        endpoints_tilde_e=find(row_sC_e_tar);
%         ep_tilde=0;
%         for ff=1:size(row_sC_e_tar,1)
%             ep_tilde=ep_tilde+1;
%             f=row_sC_e_tar(ff);
%             endpoints_tilde_e(ep_tilde)=f;
%         ends
        % end find endpoints of tilde_e from sC
        
        for i=1:denominator_e_tar
            p_tar=p_tar+1;
            timing_p_target=(i-0.5)/denominator_e_tar;
            Length_e_tar=sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec);
            Area_p_tar=Length_e_tar*cdt/denominator_e_tar;
            Length_tilde_e=sqrt(edgevec.dual(e_tar).vec.'*edgevec.dual(e_tar).vec);
            kappa(p_tar) = Length_tilde_e/Area_p_tar;
            
            stn_FutureEdge_ep1=first_stn_n(endpoints_e(1))+i;
            stn_FutureEdge_ep2=first_stn_n(endpoints_e(2))+i;
            
            df1=denominator.f(endpoints_tilde_e(1));
            % calculate difference_time
            for j=1:df1
                % can be written without if
                if (1.0/df1)*(j-1)< timing_p_target && timing_p_target < (1.0/df1)*j
                    timing_ep1tilde=(1.0/df1)*(j-0.5);
                    break;
                end
            end
            df2=denominator.f(endpoints_tilde_e(1));
            for j=1:df2
                if (1.0d0/df2)*(j-1)< timing_p_target && timing_p_target < (1.0d0/df2)*j
                    timing_ep2tilde=(1.0/df2)*(j-0.5);
                    break;
                end
            end
            
            DeltaTime_etilde=cdt*(...
                sC(endpoints_tilde_e(1),e_tar)*timing_ep1tilde+sC(endpoints_tilde_e(2),e_tar)*timing_ep2tilde);
            DeltaArea=DeltaTime_etilde/kappa(p_tar); % N.B. delta_area has sign
            delta_position=DeltaArea/(Length_e_tar*sintheta);          
            stnInfo(stn_FutureEdge_ep1).NodeDeltaVec=...
                stnInfo(stn_FutureEdge_ep1).NodeDeltaVec+delta_position*DirVec_NodeDelta;
            stnInfo(stn_FutureEdge_ep2).NodeDeltaVec=...
                stnInfo(stn_FutureEdge_ep2).NodeDeltaVec+delta_position*DirVec_NodeDelta;
            
            LorenzInvariant_dS(p_tar)=-Length_tilde_e^2+DeltaTime_etilde^2;
            kappa(p_tar)=sign(LorenzInvariant_dS(p_tar))*kappa(p_tar);            
        end 
        LorenzInvariant_dS(p_tar-denominator_e_tar)=LorenzInvariant_dS(p_tar);
        kappa(p_tar-denominator_e_tar)=kappa(p_tar);
    end % if 
end
%% calculate kappa for non-boundary edges
disp('calculating kappa for non-boundary edges')
for e_tar=1:MeshNum.E  
    if att.bound_e(e_tar)==true
         continue;
    end

    %row_sC_e_tar=find(sC(:,e_tar));

   % disp(e_target)
   % disp('is not boundary')
    % check exception-attribute
    if exception==1 %non-boundary exception (non-orthogonal)
    else %non-exceptional, non-boundary
        denominator_e_tar=denominator.e(e_tar);
        p_tar=first_pIdx.e(e_tar);
        
        endpoints_e = endpoint(e_tar,sG);

        dn1=denominator.n(endpoints_e(1));
        dn2=denominator.n(endpoints_e(2));
        stn_PastEdge_ep1=first_stn_n(endpoints_e(1))-dn1/denominator_e_tar;
        stn_PastEdge_ep2=first_stn_n(endpoints_e(2))-dn2/denominator_e_tar;
        for i=1:denominator_e_tar
            p_tar=p_tar+1;
%            timing_p_target=(1.0d0/denominator_e_tar)*(i+0.5d0);
            stn_PastEdge_ep1=stn_PastEdge_ep1+dn1/denominator_e_tar;
            stn_PastEdge_ep2=stn_PastEdge_ep2+dn2/denominator_e_tar;
                        
            BulgeArea1=CalcBulgeArea(stn_PastEdge_ep1,endpoints_e(1),e_tar,edgevec,sG,denominator,stnInfo,cdt);
            BulgeArea2=CalcBulgeArea(stn_PastEdge_ep2,endpoints_e(2),e_tar,edgevec,sG,denominator,stnInfo,cdt);

            Length_e_tar=sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec);
            Area_p_tar =Length_e_tar*cdt/denominator_e_tar + (BulgeArea1+BulgeArea2);
            Length_tilde_e=sqrt(edgevec.dual(e_tar).vec.'*edgevec.dual(e_tar).vec);
            kappa(p_tar) = Length_tilde_e/Area_p_tar;
            
            %LorenzInvariant_dS(p_tar)=-length_tilde_e^2;            
            
            % this two lines must be alternatively commented out 
            % kappa(p_target) = sign(LorenzInvariant_dS(p_tar))*kappa(p_target);
            kappa(p_tar) = -kappa(p_tar);
            
        end %loop for i
        kappa(p_tar-denominator_e_tar)=kappa(p_tar);
        LorenzInvariant_dS(p_tar-denominator_e_tar)=LorenzInvariant_dS(p_tar);
    end % if for exception
end % loop for e_target

%% calculate kappa for non-boundary faces
disp('calculating kappa for faces')
for f_tar=1:MeshNum.F
   % disp(f_target)
    %sC_f_tar=sC(f_tar,:);
    %col_sC_f_tar=find(sC_f_tar);
 
    if exception==1
    else
        denominator_f_tar=denominator.f(f_tar);
        p_tar=first_pIdx.f(f_tar)-1;

        for i=0:denominator_f_tar-1
            p_tar=p_tar+1;
            Area_p_tar=CalArea_p_f(f_tar,i,stnInfo,sC,sG,denominator,edgevec,first_stn_n);
            kappa(p_tar)=(cdt/denominator_f_tar)/Area_p_tar;
            %stdistance(p_target)=(cdt/denominator_f_target)^2;
            %kappa(p_target) = sign(stdistance(p_target))*kappa(p_target);
        end %loop
        %stdistance(first_p_for_f(f_target)+denominator_f_target)=stdistance(first_p_for_f(f_target));
        kappa(first_pIdx.f(f_tar)+denominator_f_tar)=kappa(first_pIdx.f(f_tar));
        b_area(f_tar)=(cdt/denominator_f_tar)/kappa(first_pIdx.f(f_tar));
    end % if for exceptions
end % loop for f_target

%% STOP

disp('Constitutive:ENDED')


end

function [endpoints_e] = endpoint(e_tgt,sG)
endpoints_e=zeros(2,1);
% find endpoints n1,n2 of e_target from sG
sG_e_target=sG(e_tgt,:);
col_sG_e_target=find(sG_e_target);
for jj=1:size(col_sG_e_target,2)
    endpoints_e(jj)=col_sG_e_target(jj);
end
end

function [e_connec2ntar]=e_connec2ntar_shareincfwith_e_tar(n_tar,e_tar,sC,sG,edgevec)
e_connec2ntar=struct('eIdx',[],'ShareIncf',[],'DirVec',[]);
sG_n_tar=sG(:,n_tar);
row_sG_n_tar=find(sG_n_tar);
e_con=1;
for ee=1:size(row_sG_n_tar,1) % loop for e
    e=row_sG_n_tar(ee);
    if e==e_tar  % find e that connects to endpoints_e(1)
        continue
    end
    sC_e=sC(:,e);
    row_sC_e=find(sC_e);
    for ff=1:size(row_sC_e,1) % loop for f
        f=row_sC_e(ff);
        if sC(f,e_tar)~=0 % if sC(f,e)~=0 && sC(f,e_target)~=0
            e_connec2ntar(e_con).eIdx=e;
            e_connec2ntar(e_con).ShareIncf=f;
            e_connec2ntar(e_con).DirVec=DirectionVec(e,edgevec);
            e_con=e_con+1;
            break;
        end %if
    end %for
    if (e_con>2) % if both of e was found then exit the whole loop1
        break; %loop1
    end %if
end %for loop1
end

function DirectionVecRetu=DirectionVec(e_tar,edgevec)
%compute direction vectors
DirectionVecRetu=edgevec.prim(e_tar).vec;
DirectionVecRetu=(1.0/sqrt( DirectionVecRetu.'*DirectionVecRetu ))*DirectionVecRetu;
end

function [Parallel]=ParaCheck_e_cons(e_connec_to_ep1,e_connec_to_ep2)
global EPSILON

Parallel=struct('is',false,'e_con_ep1_Idx',[],'e_con_ep2_Idx',[]);
for e_con_ep2=1:2
    for e_con_ep1=1:2
        if e_connec_to_ep1(e_con_ep1).ShareIncf==e_connec_to_ep2(e_con_ep2).ShareIncf
            ErrVec=e_connec_to_ep1(e_con_ep1).DirVec-e_connec_to_ep2(e_con_ep2).DirVec;
            if abs(ErrVec)<EPSILON
                Parallel.is=true;
                Parallel.e_con_ep1_Idx=e_con_ep1;
                Parallel.e_con_ep2_Idx=e_con_ep2;
                break
            end
        end
    end
    if  Parallel.is==true
        break
    end
end
end


function BulgeArea=CalcBulgeArea(stn_PastEdge_ep1,n_tar,e_tar,edgevec,sG,denominator,stnInfo,cdt)
% calculate CutOffArea for ep1 side
DirVec_e_tar=(1.0/sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec))*edgevec.prim(e_tar).vec;
DirVec_e_tar = sG(e_tar,n_tar)*DirVec_e_tar;
BulgeArea=0.0;
dn1   =denominator.n(n_tar);
de_tar=denominator.e(e_tar);
for j=1:dn1/de_tar
    LowerEdgeLength = stnInfo(stn_PastEdge_ep1+j-1).NodeDeltaVec.'*DirVec_e_tar;
    UpperEdgeLength = stnInfo(stn_PastEdge_ep1+j  ).NodeDeltaVec.'*DirVec_e_tar;
    BulgeArea=BulgeArea+(cdt/dn1)*(LowerEdgeLength+UpperEdgeLength)/2.0;
end
end

function Area_p_tar=CalArea_p_f(f_tar,timesection,stnInfo,sC,sG,denominator,edgevec,first_stn_n)
Area_p_tar=0.0;
denominator_f_tar=denominator.f(f_tar);

sC_f_tar=sC(f_tar,:);
% find edge1
LastEdge=find(sC_f_tar,1);
Sum_EdgeVec_Last=edgevec.prim(LastEdge).vec;
sG_InitEdge=sG(LastEdge,:);
n_start_NewEdge=find(sG_InitEdge==-1);
%n_origin=n_start_NewEdge;
n_tar_NewEdge=find(sG_InitEdge==1);
stn_start_NewEdge=first_stn_n(n_start_NewEdge)+timesection*denominator.n(n_start_NewEdge)/denominator_f_tar;
stn_tar_NewEdge=first_stn_n(n_tar_NewEdge)+timesection*denominator.n(n_tar_NewEdge)/denominator_f_tar;
% add delta to edge1_vector
Sum_EdgeVec_Last=Sum_EdgeVec_Last...
    +sG(LastEdge,n_start_NewEdge)*stnInfo(stn_start_NewEdge).NodeDeltaVec...
    +sG(LastEdge,n_tar_NewEdge)*stnInfo(stn_tar_NewEdge).NodeDeltaVec;
% end add delta to edge1_vector

n_start_NewEdge=n_tar_NewEdge;
%sG_node1=sG(:,n_start_NewEdge);
%row_sG_node1=find(sG_node1);
stn_start_NewEdge=stn_tar_NewEdge;

while true
    sC_f_tar(LastEdge)=0;
    col_sC_f_tar=find(sC_f_tar);
    if size(col_sC_f_tar,2)==0
        break
    end
    for ee=1:size(col_sC_f_tar,2)
        e=col_sC_f_tar(ee);
        if sG(e,n_start_NewEdge)~=0
            NewEdge=e;
            NewEdge_vector=edgevec.prim(NewEdge).vec;
            break
        end
    end
    
    
    
%     for ee=1:size(row_sG_node1,1)
%         e=row_sG_node1(ee);
%         if e~=LastEdge && sC(f_tar,e)~=0 % && sG(e,node1)~=0
%             NewEdge=e;
%             NewEdge_vector=edgevec.prim(NewEdge).vec;
%             break;
%         end
%     end
    
    sG_NewEdge=sG(NewEdge,:);
    sG_NewEdge(n_start_NewEdge)=0;
    n_tar_NewEdge=find(sG_NewEdge);

    
%     if n_tar_NewEdge==n_origin
%         break
%     end
    
    stn_tar_NewEdge=first_stn_n(n_tar_NewEdge)+timesection*denominator.n(n_tar_NewEdge)/denominator_f_tar;
    
    % add delta to edge2_vector
    NewEdge_vector=NewEdge_vector...
        +sG(NewEdge,n_start_NewEdge)*stnInfo(stn_start_NewEdge).NodeDeltaVec...
        +sG(NewEdge,n_tar_NewEdge)*stnInfo(stn_tar_NewEdge).NodeDeltaVec;
    % end add delta to edge2_vector
    
    % correct the direction of edge2_vector
    NewEdge_vector=sG(NewEdge,n_tar_NewEdge)*NewEdge_vector;
    
    Sum_EdgeVec_New=Sum_EdgeVec_Last+NewEdge_vector;
    
    % adding (1/2)*(cross product)_z
    Area_p_tar=Area_p_tar+0.5*Sum_EdgeVec_Last.'*[0 1;-1 0]*Sum_EdgeVec_New;
    
    n_start_NewEdge=n_tar_NewEdge;
    %sG_node1=sG(:,n_start_NewEdge);
    %row_sG_node1=find(sG_node1);
    LastEdge=NewEdge;
    %sG_edge1=sG(edge1,:);
    %col_sG_edge1=find(sG_edge1);
    stn_start_NewEdge=stn_tar_NewEdge;
    Sum_EdgeVec_Last=Sum_EdgeVec_New;
    
end % while for the boundary edges of f_target
Area_p_tar=abs(Area_p_tar);
end
