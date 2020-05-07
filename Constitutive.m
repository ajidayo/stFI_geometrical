function [kappa,b_area,att,MeshNum] ...
    =Constitutive(cdt,sC,sG,UpdateNum,edgevec,first_pIdx,att,MeshNum)
%disp('Constitutive: CALLED')

global EPSILON
global DIM

exception=0;% dummy
%% Parameters
%cdt=0.3;
%% Allocation
%endpoints_tilde_e=zeros(2,1);

att.e.nonorthogonal=logical(sparse(MeshNum.E,1));
att.e.corner=logical(sparse(MeshNum.E,1));
for e_tar=1:MeshNum.E
   if DirectionVecPrim(e_tar,edgevec).'*DirectionVecDual(e_tar,edgevec)>EPSILON
       disp(['nonorthogonal edge found:e=',num2str(e_tar),...
           'cos of the edge pair =',...
           num2str(DirectionVecPrim(e_tar,edgevec).'*DirectionVecDual(e_tar,edgevec))])
       pause(1)
       att.e.nonorthogonal(e_tar)=true;
   end
end

DirVec_NodeDelta=zeros(DIM,1); %(xy)

first_stn_n=zeros(MeshNum.N,1);
kappa=zeros(MeshNum.P,1);
LorenzInvariant_dS=sparse(MeshNum.P,1);

b_area=zeros(MeshNum.F,1);
%% generate st nodes

stn=1;
for n=1:MeshNum.N
    first_stn_n(n)=stn;
    stn=stn+UpdateNum.n(n)+1;
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
%disp('calculating kappa for boundary edges')
for e_tar=1:MeshNum.E
    if att.e.boundaryedge(e_tar)==false || att.e.nonorthogonal(e_tar)==true
         continue;
    end

    row_sC_e_tar=find(sC(:,e_tar));
  
    % check exception-attribute
    if exception==1 %exception (corner)
    else %non-exceptional dt-boundary edge
        UpdNum_e_tar=UpdateNum.e(e_tar);
        p_tar=first_pIdx.e(e_tar);        

        endpoints_e_tar = endpoint(e_tar,sG);
        endpoints_e_tar_tilde=find(row_sC_e_tar);

        DirVec_e_tar=DirectionVecPrim(e_tar,edgevec);
        %DirVec_etilde_tar=DirectionVecDual(e_tar,edgevec);
                
        e_connec_to_ep1=e_connec2ntar_shareincfwith_e_tar(endpoints_e_tar(1),e_tar,sC,sG,edgevec);
        e_connec_to_ep2=e_connec2ntar_shareincfwith_e_tar(endpoints_e_tar(2),e_tar,sC,sG,edgevec);

        % correct directions of vec_e1c, vec_e2c: the direction is the same
        % as the direction of edgevec.dual.(the direction of g, or tilde e)
        for e_con=1:2
            e_connec_to_ep1(e_con).DirVec=e_connec_to_ep1(e_con).DirVec...
                *sign(DirVec_e_tar.'*e_connec_to_ep1(e_con).DirVec);
            e_connec_to_ep2(e_con).DirVec=e_connec_to_ep2(e_con).DirVec...
                *sign(DirVec_e_tar.'*e_connec_to_ep2(e_con).DirVec);
        end
        
        Parallel=ParaCheck_e_cons(e_connec_to_ep1,e_connec_to_ep2);
        if Parallel.is==true
            DirVec_NodeDelta=e_connec_to_ep1(Parallel.e_con_ep1_Idx).DirVec;
        else
            % maybe corner of the subgrid
            disp(['error; Cannot define DirVec_NodeDelta for e_tar = ', num2str(e_tar)])
            att.e.corner(e_tar)=true;
            continue
        end
        costheta=DirVec_NodeDelta.'*DirVec_e_tar;
        sintheta=sqrt(1-costheta^2);
        
        for timesection=1:UpdNum_e_tar
            p_tar=p_tar+1;
            timing_p_target=(timesection-0.5)/UpdNum_e_tar;
            Length_e_tar=sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec);
            Area_p_tar=Length_e_tar*cdt/UpdNum_e_tar;
            Length_tilde_e_tar=sqrt(edgevec.dual(e_tar).vec.'*edgevec.dual(e_tar).vec);
            kappa(p_tar) = Length_tilde_e_tar/Area_p_tar;
            
            stn_FutureEdge_ep1=first_stn_n(endpoints_e_tar(1))+timesection;
            stn_FutureEdge_ep2=first_stn_n(endpoints_e_tar(2))+timesection;
            
            UdNf1=UpdateNum.f(endpoints_e_tar_tilde(1));
            % calculate difference_time
            for j=1:UdNf1
                % can be written without if
                if (1.0/UdNf1)*(j-1)< timing_p_target && timing_p_target < (1.0/UdNf1)*j
                    timing_ep1tilde=(1.0/UdNf1)*(j-0.5);
                    break;
                end
            end
            UdNf2=UpdateNum.f(endpoints_e_tar_tilde(1));
            for j=1:UdNf2
                if (1.0d0/UdNf2)*(j-1)< timing_p_target && timing_p_target < (1.0d0/UdNf2)*j
                    timing_ep2tilde=(1.0/UdNf2)*(j-0.5);
                    break;
                end
            end
            
            DeltaTime_etilde=cdt*(...
                sC(endpoints_e_tar_tilde(1),e_tar)*timing_ep1tilde+sC(endpoints_e_tar_tilde(2),e_tar)*timing_ep2tilde);
            DeltaArea=DeltaTime_etilde/kappa(p_tar); % N.B. delta_area has sign
            delta_position=DeltaArea/(Length_e_tar*sintheta);          
            stnInfo(stn_FutureEdge_ep1).NodeDeltaVec=...
                stnInfo(stn_FutureEdge_ep1).NodeDeltaVec+delta_position*DirVec_NodeDelta;
            stnInfo(stn_FutureEdge_ep2).NodeDeltaVec=...
                stnInfo(stn_FutureEdge_ep2).NodeDeltaVec+delta_position*DirVec_NodeDelta;
            
            LorenzInvariant_dS(p_tar)=-Length_tilde_e_tar^2+DeltaTime_etilde^2;
            kappa(p_tar)=sign(LorenzInvariant_dS(p_tar))*kappa(p_tar);            
        end 
        LorenzInvariant_dS(p_tar-UpdNum_e_tar)=LorenzInvariant_dS(p_tar);
        kappa(p_tar-UpdNum_e_tar)=kappa(p_tar);
    end % if 
end
%% calculate kappa for corners

for e_tar=1:MeshNum.E
    if att.e.boundaryedge(e_tar)==false || att.e.corner(e_tar)==false
         continue;
    end
    UpdNum_e_tar=UpdateNum.e(e_tar);
    p_tar=first_pIdx.e(e_tar);
        
    endpoints_e_tar = endpoint(e_tar,sG);
    % find stn from which we fetch NodeDeltaVec
    e_adjto_e_tar = find(sG(e_tar,:)*sG.');
    for ee=1:size(e_adjto_e_tar,2) % all adjacent edges
        e = e_adjto_e_tar(ee);
        if att.e.corner(e)==false
            endpoints_e = endpoint(e,sG);
            for ep=1:2
                if sG(e_tar,endpoints_e(ep))~=0
                    nfetch=endpoints_e(ep);
                end
            end
        end
    end  
    for timesection=1:UpdNum_e_tar
        p_tar=p_tar+1;
        stn_FutureEdge_ep1=first_stn_n(endpoints_e_tar(1))+timesection;
        stn_FutureEdge_ep2=first_stn_n(endpoints_e_tar(2))+timesection;
        stn_FutureEdge_nfetch=first_stn_n(nfetch)+timesection;
        stnInfo(stn_FutureEdge_ep1).NodeDeltaVec=...
            stnInfo(stn_FutureEdge_nfetch).NodeDeltaVec;
        stnInfo(stn_FutureEdge_ep2).NodeDeltaVec=...
            stnInfo(stn_FutureEdge_nfetch).NodeDeltaVec;
    end
end

for e_tar=1:MeshNum.E
    if att.e.boundaryedge(e_tar)==false || att.e.corner(e_tar)==false
         continue;
    end
    UpdNum_e_tar=UpdateNum.e(e_tar);
    p_tar=first_pIdx.e(e_tar);
    
    endpoints_e = endpoint(e_tar,sG);
    endpoints_etilde=find(row_sC_e_tar);
        
    for timesection=1:UpdNum_e_tar
            p_tar=p_tar+1;
            timing_p_target=(timesection-0.5)/UpdNum_e_tar;

            for ep =1:2
                if sG(e_tar,endpoints_e(ep))==-1
                    nstart=endpoints_e(ep);
                else 
                    ntarget=endpoints_e(ep);
                end
            end
            stn_FutureEdge_nstart=first_stn_n(nstart)+timesection;
            stn_FutureEdge_ntarget=first_stn_n(ntarget)+timesection;
            stn_PastEdge_nstart=first_stn_n(nstart)+timesection-1;
            stn_PastEdge_ntarget=first_stn_n(ntarget)+timesection-1;
            
            Futureedgevec_prim_e_tar = edgevec.prim(e_tar);
            Pastedgevec_prim_e_tar = edgevec.prim(e_tar);
            Futureedgevec_prim_e_tar=Futureedgevec_prim_e_tar...
                +stnInfo(stn_FutureEdge_nstart).NodeDeltaVec...
                +stnInfo(stn_FutureEdge_ntarget).NodeDeltaVec;
            Patsedgevec_prim_e_tar=Pastedgevec_prim_e_tar...
                +stnInfo(  stn_PastEdge_nstart).NodeDeltaVec...
                +stnInfo(  stn_PastEdge_ntarget).NodeDeltaVec;

            OrthogLength_FtrEdge ...
                =Futureedgevec_prim_e_tar.'*[0 1;-1 0]*DirectionVecDual(e_tar,edgevec);
            OrthogLength_PstEdge ...
                =  Patsedgevec_prim_e_tar.'*[0 1;-1 0]*DirectionVecDual(e_tar,edgevec);
            stP = 0.5*(OrthogLength_FtrEdge+OrthogLength_PstEdge)*cdt/UpdNum_e_tar; 
          
            ssP = ...
                abs(Pastedgevec_prim_e_tar.'*[0 1;-1 0]*...
                ( stnInfo(stn_FutureEdge_nstart).NodeDeltaVec...
                - stnInfo(  stn_PastEdge_nstart).NodeDeltaVec)...
                ) + abs(-Futureedgevec_prim_e_tar.'*[0 1;-1 0]*...
                ( stnInfo(  stn_PastEdge_ntarget).NodeDeltaVec )...
                - stnInfo(stn_FutureEdge_ntarget).NodeDeltaVec...
                );
            P = [stP;ssP];

            UdNf1=UpdateNum.f(endpoints_etilde(1));
            % calculate difference_time
            for j=1:UdNf1
                % can be written without if
                if (1.0/UdNf1)*(j-1)< timing_p_target && timing_p_target < (1.0/UdNf1)*j
                    timing_ep1tilde=(1.0/UdNf1)*(j-0.5);
                    break;
                end
            end
            UdNf2=UpdateNum.f(endpoints_etilde(1));
            for j=1:UdNf2
                if (1.0d0/UdNf2)*(j-1)< timing_p_target && timing_p_target < (1.0d0/UdNf2)*j
                    timing_ep2tilde=(1.0/UdNf2)*(j-0.5);
                    break;
                end
            end
            DeltaTime_etilde=cdt*(...
                sC(endpoints_etilde(1),e_tar)*timing_ep1tilde+sC(endpoints_etilde(2),e_tar)*timing_ep2tilde);
            edgevec_dual_e_tar = edgevec.dual(e_tar);
            Length_tilde_e = edgevec_dual_e_tar.'*edgevec_dual_e_tar;
            Direction_tildeP = [DeltaTime_etilde;Length_tilde_e];
            Direction_tildeP = (1.0/sqrt(Direction_tildeP.'*Direction_tildeP))...
             *Direction_tildeP;
         
            Area_p_tar=P.'*Direction_tildeP;
            kappa(p_tar) = sqrt(Direction_tildeP.'*Direction_tildeP)/Area_p_tar;

            LorenzInvariant_dS(p_tar)=-Length_tilde_e^2+DeltaTime_etilde^2;
            kappa(p_tar)=sign(LorenzInvariant_dS(p_tar))*kappa(p_tar);
    end
    LorenzInvariant_dS(p_tar-UpdNum_e_tar)=LorenzInvariant_dS(p_tar);
    kappa(p_tar-UpdNum_e_tar)=kappa(p_tar);
end

%% calculate kappa for non-boundary edges
%disp('calculating kappa for non-boundary edges')
for e_tar=1:MeshNum.E  
    if att.e.boundaryedge(e_tar)==true
         continue;
    end

    if exception==1 
    else %non-exceptional, non-boundary
        UpdNum_e_tar=UpdateNum.e(e_tar);
        p_tar=first_pIdx.e(e_tar);
        
        endpoints_e = endpoint(e_tar,sG);

        UdNn1=UpdateNum.n(endpoints_e(1));
        UdNdn2=UpdateNum.n(endpoints_e(2));
        stn_PastEdge_ep1=first_stn_n(endpoints_e(1))-UdNn1/UpdNum_e_tar;
        stn_PastEdge_ep2=first_stn_n(endpoints_e(2))-UdNdn2/UpdNum_e_tar;
        for timesection=1:UpdNum_e_tar
            p_tar=p_tar+1;
            stn_PastEdge_ep1=stn_PastEdge_ep1+UdNn1/UpdNum_e_tar;
            stn_PastEdge_ep2=stn_PastEdge_ep2+UdNdn2/UpdNum_e_tar;
                        
            BulgeArea1=CalcBulgeArea(stn_PastEdge_ep1,endpoints_e(1),e_tar,edgevec,sG,UpdateNum,stnInfo,cdt);
            BulgeArea2=CalcBulgeArea(stn_PastEdge_ep2,endpoints_e(2),e_tar,edgevec,sG,UpdateNum,stnInfo,cdt);

            Length_e_tar=sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec);
            Area_p_tar =Length_e_tar*cdt/UpdNum_e_tar + (BulgeArea1+BulgeArea2);
            Length_tilde_e=sqrt(edgevec.dual(e_tar).vec.'*edgevec.dual(e_tar).vec);
            
            
            sin_p2d_e=DirectionVecPrim(e_tar,edgevec).'*[0 1;-1 0]*DirectionVecDual(e_tar,edgevec);
         
            Area_p_tar=sin_p2d_e*Area_p_tar;
            kappa(p_tar) = Length_tilde_e/Area_p_tar;
            %LorenzInvariant_dS(p_tar)=-length_tilde_e^2;            
            
            % this two lines must be alternatively commented out 
            % kappa(p_target) = sign(LorenzInvariant_dS(p_tar))*kappa(p_target);
            kappa(p_tar) = -kappa(p_tar);
            
        end %loop for i
        kappa(p_tar-UpdNum_e_tar)=kappa(p_tar);
        %LorenzInvariant_dS(p_tar-denominator_e_tar)=LorenzInvariant_dS(p_tar);
    end % if for exception
end % loop for e_target

%% calculate kappa for non-boundary faces
%disp('calculating kappa for faces')
for f_tar=1:MeshNum.F
   
    if exception==1
    else
        UdN_f_tar=UpdateNum.f(f_tar);
        p_tar=first_pIdx.f(f_tar)-1;

        for timesection=0:UdN_f_tar-1
            p_tar=p_tar+1;
            Area_p_tar=CalArea_p_f(f_tar,timesection,stnInfo,sC,sG,UpdateNum,edgevec,first_stn_n);
            kappa(p_tar)=(cdt/UdN_f_tar)/Area_p_tar;
            %LorenzInvariant_dS(p_tar)=(cdt/denominator_f_target)^2;
            %kappa(p_target) = sign(LorenzInvariant_dS(p_target))*kappa(p_target);
        end %loop
        %LorenzInvariant_dS(first_pIdx.f(f_tar)+denominator_f_tar)=LorenzInvariant_dS(first_pIdx.f(f_target));
        kappa(first_pIdx.f(f_tar)+UdN_f_tar)=kappa(first_pIdx.f(f_tar));
        b_area(f_tar)=(cdt/UdN_f_tar)/kappa(first_pIdx.f(f_tar));
    end % if 
end %  f_tar

%% STOP

%disp('Constitutive:ENDED')


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
            e_connec2ntar(e_con).DirVec=DirectionVecPrim(e,edgevec);
            e_con=e_con+1;
            break;
        end %if
    end %for
    if (e_con>2) % if both of e was found then exit the whole loop1
        break; %loop1
    end %if
end %for loop1
end

function DirectionVecRetu=DirectionVecPrim(e_tar,edgevec)
%compute direction vectors
DirectionVecRetu=edgevec.prim(e_tar).vec;
DirectionVecRetu=(1.0/sqrt( DirectionVecRetu.'*DirectionVecRetu ))*DirectionVecRetu;
end


function DirectionVecRetu=DirectionVecDual(e_tar,edgevec)
%compute direction vectors
DirectionVecRetu=edgevec.dual(e_tar).vec;
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


function BulgeArea=CalcBulgeArea(stn_PastEdge_ep1,n_tar,e_tar,edgevec,sG,UpdaneNum,stnInfo,cdt)
% calculate CutOffArea for ep1 side
DirVec_e_tar=(1.0/sqrt(edgevec.prim(e_tar).vec.'*edgevec.prim(e_tar).vec))*edgevec.prim(e_tar).vec;
DirVec_e_tar = sG(e_tar,n_tar)*DirVec_e_tar;
BulgeArea=0.0;
UdNn1   =UpdaneNum.n(n_tar);
UdNe_tar=UpdaneNum.e(e_tar);
for j=1:UdNn1/UdNe_tar
    LowerEdgeLength = stnInfo(stn_PastEdge_ep1+j-1).NodeDeltaVec.'*DirVec_e_tar;
    UpperEdgeLength = stnInfo(stn_PastEdge_ep1+j  ).NodeDeltaVec.'*DirVec_e_tar;
    BulgeArea=BulgeArea+(cdt/UdNn1)*(LowerEdgeLength+UpperEdgeLength)/2.0;
end
end

function Area_p_tar=CalArea_p_f(f_tar,timesection,stnInfo,sC,sG,UpdateNum,edgevec,first_stn_n)
Area_p_tar=0.0;
UdNf_tar=UpdateNum.f(f_tar);

sC_f_tar=sC(f_tar,:);
% find edge1
LastEdge=find(sC_f_tar,1);
Sum_EdgeVec_Last=edgevec.prim(LastEdge).vec;
sG_InitEdge=sG(LastEdge,:);
n_start_NewEdge=find(sG_InitEdge==-1);
%n_origin=n_start_NewEdge;
n_tar_NewEdge=find(sG_InitEdge==1);
stn_start_NewEdge=first_stn_n(n_start_NewEdge)+timesection*UpdateNum.n(n_start_NewEdge)/UdNf_tar;
stn_tar_NewEdge=first_stn_n(n_tar_NewEdge)+timesection*UpdateNum.n(n_tar_NewEdge)/UdNf_tar;
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
    
    stn_tar_NewEdge=first_stn_n(n_tar_NewEdge)+timesection*UpdateNum.n(n_tar_NewEdge)/UdNf_tar;
    
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
%disp('Constitutive: ENDED')
end
