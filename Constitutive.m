function [kappa,b_area,MeshNum] ...
    =Constitutive(cdt,sC,sG,denominator,edgevec,first_p,tilde_node_position,MeshNum)

global EPSILON
global dim

denominator_face=denominator.f;
denominator_edge=denominator.e;
denominator_node=denominator.n;

edge_vector=edgevec.prim;
tilde_edge_vector=edgevec.dual;

first_p_for_f=first_p.f;
first_p_for_e=first_p.e;

exception=0;% dummy
%% Parameters
%cdt=0.3;
%% Allocation
endpoints_e=zeros(2,1);
endpoints_tilde_e=zeros(2,1);
error=zeros(2,1);
e1c=zeros(2,1);
e2c=zeros(2,1);
vec_e1c=zeros(2,dim);vec_e2c=zeros(2,dim); % (n1 or n2,xy)
vec_e_target=zeros(dim,1); %(xy)
direction_vector=zeros(dim,1); %(xy)
edge1_vector=zeros(dim,1);
edge2_vector=zeros(dim,1);
edge1plus2_vector=zeros(dim,1);

first_stn_for_n=zeros(MeshNum.N,1);
%z=zeros(MeshNum.P,1);
kappa=zeros(MeshNum.P,1);
stdistance=zeros(MeshNum.P,1);

b_area=zeros(MeshNum.F,1);
%% generate st nodes

stn=1;
for n=1:MeshNum.N
    first_stn_for_n(n)=stn;
    stn=stn+denominator_node(n)+1;
end
stNNum=stn-1;
denominator.stN=stNNum;

stnode_delta_vector=sparse(stNNum,dim);

%% calculate kappa for boundary edges
disp('calculating kappa for boundary edges')
for e_target=1:MeshNum.E
    sC_e_target=sC(:,e_target);
    row_sC_e_target=find(sC_e_target);
    % check boundary-attribute
    check_boundary=0.0;
    for ff=1:size(row_sC_e_target,1)
        f=row_sC_e_target(ff);
        check_boundary=check_boundary+sC(f,e_target)*denominator_face(f);
    end
    if abs(check_boundary)<EPSILON % not a boundary
         continue;
    end
%    disp('edge:')
 %   disp(e_target)
  %  disp('is boundary; next to')
%    for ff=1:size(row_sC_e_target,1)
 %       f=row_sC_e_target(ff)
  %  end
  
    % check exception-attribute
    if exception==1 %exception (corner)
    else %boundary_attribute(e_target) .eqv. .TRUE.%non-exceptional boundary edge
        denominator_e_target=denominator_edge(e_target);
        p_target=first_p_for_e(e_target);
        
        % find endpoints n1,n2 of e from sG
        sG_e_target=sG(e_target,:);
        col_sG_e_target=find(sG_e_target);
        for jj=1:size(col_sG_e_target,2)
            endpoints_e(jj)=col_sG_e_target(jj);
        end
        % find e1c(1), e1c(2)
        
        %loop1
        sG_endpoint1=sG(:,endpoints_e(1));
        row_sG_endpoint1=find(sG_endpoint1);
        j=1;
        for ee=1:size(row_sG_endpoint1,1) % loop for e
            e=row_sG_endpoint1(ee);
            if e==e_target  % find e that connects to endpoints_e(1)
                continue
            end
            sC_e=sC(:,e);
            row_sC_e=find(sC_e);
            for ff=1:size(row_sC_e,1) % loop for f
                f=row_sC_e(ff);
                if sC(f,e_target)~=0 % if sC(f,e)~=0 && sC(f,e_target)~=0
                    e1c(j)=e;
                    j=j+1;
                    break;
                end %if
            end %for
            if (j>2) % if both of e was found then exit the whole loop1
                break; %loop1
            end %if
        end %for loop1
        % end find e1c(1), e1c(2)
        
        % find e2c(1), e2c(2)
        %loop1
        sG_endpoint2=sG(:,endpoints_e(2));
        row_sG_endpoint2=find(sG_endpoint2);
        j=1;
        for ee=1:size(row_sG_endpoint2,1) % loop for e
            e=row_sG_endpoint2(ee);
            if e==e_target  % find e that connects to endpoints_e(1)
                continue
            end
            sC_e=sC(:,e);
            row_sC_e=find(sC_e);
            for ff=1:size(row_sC_e,1) % loop for f
                f=row_sC_e(ff);
                if sC(f,e_target)~=0 % if sC(f,e)~=0 && sC(f,e_target)~=0
                    if sC(f,e1c(1))~=0
                        e2c(1)=e;
                    elseif sC(f,e1c(2))~=0
                        e2c(2)=e;
                    end
                    j=j+1;
                    break;
                end %if
            end %for
            if j>2 % if both of e was found then exit the whole loop1
                break; %loop1
            end %if
        end %for loop1
        % end  find e2c(1), e2c(2)

        % compute direction vectors
        for xy=1:dim
            for j=1:2
                vec_e1c(j,xy)=edge_vector(e1c(j),xy);
                vec_e2c(j,xy)=edge_vector(e2c(j),xy);
            end
        end
        for xy=1:dim
            vec_e_target(xy)=edge_vector(e_target,xy);
        end
        for xy=1:dim
            for j=1:2
                vec_e1c(j,xy)=vec_e1c(j,xy)/sqrt( vec_e1c(j,1)^2 + vec_e1c(j,2)^2 );
                vec_e2c(j,xy)=vec_e2c(j,xy)/sqrt( vec_e2c(j,2)^2 + vec_e2c(j,1)^2 );
            end
        end
        for xy=1:dim
            vec_e_target(xy)=vec_e_target(xy)/sqrt( vec_e_target(1)^2 + vec_e_target(2)^2 );
        end
        % correct directions of vec_e1c, vec_e2c: the direction is the same
        % as the direction of tilde_edge_vector.(the direction of g, or tilde e)
        for xy=1:dim
            for j=1:2
                vec_e1c(j,xy)=vec_e1c(j,xy)...
                    *sign(vec_e1c(j,1)*tilde_edge_vector(e_target,1)+vec_e1c(j,2)*tilde_edge_vector(e_target,2));
                vec_e2c(j,xy)=vec_e2c(j,xy)...
                    *sign(vec_e2c(j,1)*tilde_edge_vector(e_target,1)+vec_e2c(j,2)*tilde_edge_vector(e_target,2));
            end
        end
        % end compute direction vectors

        % check that the two edges are parallel: if parallel then the direction of the two edges is the shift direction for boundary nodes
        error(1)=0.0;
        for xy=1:dim
            error(1)=error(1)+(vec_e1c(1,xy)-vec_e2c(1,xy))^2;
        end
        error(2)=0.0;
        for xy=1:dim
            error(2)=error(2)+(vec_e1c(2,xy)-vec_e2c(2,xy))^2;
        end
        if error(1)<EPSILON
            for xy =1:dim
                direction_vector(xy)=vec_e1c(1,xy);
            end
        end
        if error(2)<EPSILON
            for xy =1:dim
                direction_vector(xy)=vec_e1c(2,xy);
            end
        end
        costheta=0.0;
        for xy=1:dim
            costheta=costheta+direction_vector(xy)*vec_e_target(xy);
        end
        sintheta=sqrt(1-costheta^2);
        % end check parallel
        
        % find endpoints of tilde_e from sC
        j=1;
        for ff=1:size(row_sC_e_target,1)
            f=row_sC_e_target(ff);
            endpoints_tilde_e(j)=f;
            j=j+1;
            if j>2
                break;
            end
        end
        % end find endpoints of tilde_e from sC
        
        for i=0:denominator_e_target-1
            p_target=p_target+1;
            timing_p_target=(i+0.5)/denominator_e_target;
            length_e_target=sqrt(edge_vector(e_target,1)^2+edge_vector(e_target,2)^2);
            area_p_target=length_e_target*cdt/denominator_e_target;
            length_tilde_e=sqrt(tilde_edge_vector(e_target,1)^2+tilde_edge_vector(e_target,2)^2);
            kappa(p_target) = length_tilde_e/area_p_target;
            
            % set stn_1 and stn_2; the endpoints of the %%%%%%%%%%%% UPPER %%%%%%%%%%% edge of p_target
            stn_target_1=first_stn_for_n(endpoints_e(1))+i+1;
            stn_target_2=first_stn_for_n(endpoints_e(2))+i+1;
            % end set stn_1 and stn_2
        
            % calculate difference_time
            df1=denominator_face(endpoints_tilde_e(1));
            for j=0:df1-1
                % can be written without if
                if (1.0/df1)*j< timing_p_target && timing_p_target < (1.0/df1)*(j+1.0)
                    timing_1=(1.0/df1)*(j+0.5);
                    break;
                end
            end
            df2=denominator_face(endpoints_tilde_e(2));
            for j=0:df2-1
                if (1.0d0/df2)*j< timing_p_target && timing_p_target < (1.0d0/df2)*(j+1.0d0)
                    timing_2=(1.0/df2)*(j+0.5);
                    break;
                end
            end
            difference_time=...
                sC(endpoints_tilde_e(1),e_target)*timing_1+sC(endpoints_tilde_e(2),e_target)*timing_2;
            % -sC(endpoints_tilde_e(1),e_target)*timing_1-sC(endpoints_tilde_e(2),e_target)*timing_2;
            difference_time=cdt*difference_time;
            % end calculate difference_time
       
            delta_area=difference_time/kappa(p_target); % N.B. delta_area has sign
            delta_position=delta_area/(length_e_target*sintheta);          
            for xy=1:dim
                stnode_delta_vector(stn_target_1,xy)=...
                    stnode_delta_vector(stn_target_1-1,xy)+delta_position*direction_vector(xy);
            end
            for xy=1:dim
                stnode_delta_vector(stn_target_2,xy)=...
                    stnode_delta_vector(stn_target_2-1,xy)+delta_position*direction_vector(xy);
            end
            stdistance(p_target)=-length_tilde_e^2+difference_time^2;
            kappa(p_target)=sign(stdistance(p_target))*kappa(p_target);            
        end % loop for i
        stdistance(p_target-denominator_e_target)=stdistance(p_target);
        kappa(p_target-denominator_e_target)=kappa(p_target);
    end % if for exceptional or not
end
%% calculate kappa for non-boundary edges
disp('calculating kappa for non-boundary edges')
for e_target=1:MeshNum.E  
    sC_e_target=sC(:,e_target);
    row_sC_e_target=find(sC_e_target);
    % check boundary-attribute
    check_boundary=0.0;
    for ff=1:size(row_sC_e_target,1)
        f=row_sC_e_target(ff);
        check_boundary=check_boundary+sC(f,e_target)*denominator_face(f);
    end
    if abs(check_boundary)>EPSILON % if boundary
        continue;
    end
    
   % disp(e_target)
   % disp('is not boundary')
    % check exception-attribute
    if exception==1 %non-boundary exception (non-orthogonal)
    else %non-exceptional, non-boundary
        denominator_e_target=denominator_edge(e_target);
        p_target=first_p_for_e(e_target);
            
        
        % find endpoints n1,n2 of e from sG
        sG_e_target=sG(e_target,:);
        col_sG_e_target=find(sG_e_target);
        for jj=1:size(col_sG_e_target,2)
            endpoints_e(jj)=col_sG_e_target(jj);
        end
        % end find endpoints
        
        dn1=denominator_node(endpoints_e(1));
        dn2=denominator_node(endpoints_e(2));
       
        stn_target_1=first_stn_for_n(endpoints_e(1))-dn1/denominator_e_target;
        stn_target_2=first_stn_for_n(endpoints_e(2))-dn2/denominator_e_target;
                
        for i=0:denominator_e_target-1
            p_target=p_target+1;
            timing_p_target=(1.0d0/denominator_e_target)*(i+0.5d0);
            stn_target_1=stn_target_1+dn1/denominator_e_target;
            stn_target_2=stn_target_2+dn2/denominator_e_target;
        
            % calculate area1 for n1 side
            norm=0.0;
            for xy=1:dim
                direction_vector(xy)=edge_vector(e_target,xy);
                norm=norm+direction_vector(xy)^2;
            end
            norm=sqrt(norm);
            for xy=1:dim
                direction_vector(xy) = direction_vector(xy)/norm;
                direction_vector(xy) = sG(e_target,endpoints_e(1))*direction_vector(xy);
            end
            area1=0.0;
            
            for j=1:dn1/denominator_e_target
                length_lower_edge=0.0;
                length_upper_edge=0.0;
                for xy=1:dim
                    length_lower_edge = length_lower_edge...
                        + stnode_delta_vector(stn_target_1+j-1,xy)*direction_vector(xy);
                    length_upper_edge = length_upper_edge...
                        + stnode_delta_vector(stn_target_1+j  ,xy)*direction_vector(xy);
                end
                area1=area1+(cdt/dn1)*(length_lower_edge+length_upper_edge)/2.0;
            end
            % end calculate area1
            
            % calculate area2 for n2 side
            norm=0.0;
            for xy=1:dim
                direction_vector(xy)=edge_vector(e_target,xy);
                norm=norm+direction_vector(xy)^2;
            end
            norm=sqrt(norm);
            for xy=1:dim
                direction_vector(xy)=direction_vector(xy)/norm;
                direction_vector(xy) = sG(e_target,endpoints_e(2))*direction_vector(xy);
            end
            area2=0.0d0;
            for j=1:dn2/denominator_e_target
                length_lower_edge=0.0;
                length_upper_edge=0.0;
                for xy=1:dim
                    length_lower_edge = length_lower_edge...
                        + stnode_delta_vector(stn_target_2+j-1,xy)*direction_vector(xy);
                    length_upper_edge = length_upper_edge...
                        + stnode_delta_vector(stn_target_2+j  ,xy)*direction_vector(xy);
                end
                area2=area2+(cdt/dn2)*(length_lower_edge+length_upper_edge)/2.0;
            end
            % end calculate area2
            
            length_e_target=sqrt(edge_vector(e_target,1)^2+edge_vector(e_target,2)^2);
            area_p_target =length_e_target*cdt/denominator_e_target + (area1+area2);
            length_tilde_e=sqrt(tilde_edge_vector(e_target,1)^2+tilde_edge_vector(e_target,2)^2);
            kappa(p_target) = length_tilde_e/area_p_target;
            %stdistance(p_target)=-length_tilde_e^2;
            
            % this two lines must be alternatively commented out 
            % kappa(p_target) = sign(stdistance(p_target))*kappa(p_target);
            kappa(p_target) = -kappa(p_target);
            
        end %loop for i
        kappa(p_target-denominator_e_target)=kappa(p_target);
        stdistance(p_target-denominator_e_target)=stdistance(p_target);
    end % if for exception
end % loop for e_target

%% calculate kappa for non-boundary faces
disp('calculating kappa for faces')
for f_target=1:MeshNum.F
   % disp(f_target)
    sC_f_target=sC(f_target,:);
    col_sC_f_target=find(sC_f_target);
 
    if exception==1
    else
        denominator_f_target=denominator_face(f_target);
        p_target=first_p_for_f(f_target)-1;

        for i=0:denominator_f_target-1
            p_target=p_target+1;
            area_p_target=0.0;
            
            % find edge1
            edge1=col_sC_f_target(1);
            for xy = 1:dim
                edge1_vector(xy)=edge_vector(edge1,xy);
            end
            % end find edge1
            
            % find nstart==node1, node2
            sG_edge1=sG(edge1,:);
            col_sG_edge1=find(sG_edge1);
            for nn=1:size(col_sG_edge1,2)
                n=col_sG_edge1(nn);
                if sG(edge1,n)==-1
                    nstart=n;
                    node1=n;
                elseif sG(edge1,n)==1
                    node2=n;
                end
            end
            % end find nstart==node1, node2

            stn_node1=first_stn_for_n(node1)+i*denominator_node(node1)/denominator_f_target;
            stn_node2=first_stn_for_n(node2)+i*denominator_node(node2)/denominator_f_target;
            
            % add delta to edge1_vector
            for xy=1:dim
                edge1_vector(xy)=edge1_vector(xy)...
                    +sG(edge1,node1)*stnode_delta_vector(stn_node1,xy)...
                    +sG(edge1,node2)*stnode_delta_vector(stn_node2,xy);
            end
            % end add delta to edge1_vector
            
            node1=node2;
            sG_node1=sG(:,node1);
            row_sG_node1=find(sG_node1);
            
            stn_node1=stn_node2;
            
            while true
                % find edge2
                for ee=1:size(row_sG_node1,1)
                    e=row_sG_node1(ee);
                    if e~=edge1 && sC(f_target,e)~=0 % && sG(e,node1)~=0 
                        edge2=e;
                        for xy = 1:dim
                            edge2_vector(xy)=edge_vector(edge2,xy);
                        end
                        break;
                    end
                end
                % end find edge2

                sG_edge2=sG(edge2,:);
                col_sG_edge2=find(sG_edge2);
                
                % find node2
                for nn=1:size(col_sG_edge2,2)
                    n=col_sG_edge2(nn);
                    if n~=node1 %% && sG(edge2,n)~=0
                        node2 = n;
                        break
                    end
                end
                % end find node2
                
                if node2==nstart
                    break
                end
                
                stn_node2=first_stn_for_n(node2)+i*denominator_node(node2)/denominator_f_target;
                
                % add delta to edge2_vector
                for xy=1:dim
                    edge2_vector(xy)=edge2_vector(xy)...
                        +sG(edge2,node1)*stnode_delta_vector(stn_node1,xy)...
                        +sG(edge2,node2)*stnode_delta_vector(stn_node2,xy);
                end
                % end add delta to edge2_vector
                
                % correct the direction of edge2_vector
                for xy=1:dim
                    edge2_vector(xy)=sG(edge2,node2)*edge2_vector(xy);
                end
                % end correct the direction of edge2_vector
                
                for xy=1:dim
                    edge1plus2_vector(xy)=edge1_vector(xy)+edge2_vector(xy);
                end
                
                area_cross_product...
                    =edge1_vector(1)*edge1plus2_vector(2)...
                    -edge1_vector(2)*edge1plus2_vector(1);
                area_p_target=area_p_target+0.5d0*area_cross_product;
                
                node1=node2;
                sG_node1=sG(:,node1);
                row_sG_node1=find(sG_node1);
                
                edge1=edge2;
                sG_edge1=sG(edge1,:);
                col_sG_edge1=find(sG_edge1);
                
                
                
                stn_node1=stn_node2;
                for xy=1:dim
                    edge1_vector(xy)=edge1plus2_vector(xy);
                end
            end % while for the boundary edges of f_target
            
            area_p_target=abs(area_p_target);
            kappa(p_target)=(cdt/denominator_f_target)/area_p_target;
            %stdistance(p_target)=(cdt/denominator_f_target)^2;
            %kappa(p_target) = sign(stdistance(p_target))*kappa(p_target);
        end %loop for i
        %stdistance(first_p_for_f(f_target)+denominator_f_target)=stdistance(first_p_for_f(f_target));
        kappa(first_p_for_f(f_target)+denominator_f_target)=kappa(first_p_for_f(f_target));
        b_area(f_target)=(cdt/denominator_f_target)/kappa(first_p_for_f(f_target));
    end % if for exceptions
end % loop for f_target

  
%% END; FREEING MEMORIES

%clear denominator_node;

%% STOP

disp('Constitutive:END')

end