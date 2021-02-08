function [Q,iv,VT,H]=IsoContour(TR,F,iv,vis)
% Extract line segments comprising level set(s) of a scalar field defined 
% at the vertices of triangular surface mesh. 
%
% INPUT:
%   - TR    : input surface mesh represented as an object of 'TriRep' 
%             class, 'triangulation' class, or a cell such that TR={Tri,X},
%             where Tri is an M-by-3 array of faces and X is an N-by-3 
%             array of vertex coordinates. 
%   - F     : N-by-1 array defining values of the scalar field at the 
%             mesh vertices. 
%   - iv    : 1-by-K vector specifying either the values OR number of
%             level sets. In case of the latter, level sets will be 
%             distributed uniformly between the extremal values of F.
%             To get a level set at one specific value of F (e.g., Fs), 
%             specify iv as Fs*[1 1]. 
%   - vis   : axes handle (or logical value) indicating where (or whether)
%             computed level sets should be plotted. 
%
% OUTPUT: 
%   - Q     : 1-by-K cell containing coordinates of UNORDERED line segments
%             comprising the level sets. Use 'OrderIsoContourVerts' 
%             function to recover isocontours with sequentially ordered 
%             vertices.
%   - iv    : 1-by-K array of level set values corresponding to level sets
%             in Q.
%   - VT    : 1-by-K cell containing linear interpolation coefficients 
%             along the with indices of mesh vertices that form the edges 
%             intersected by the level set(s). VT{k}(m,:)=[v1 v2 t flag],
%             flag=1 if there is a boundary edge between v1 and v2 and is 
%             zero otherwise. Note, closed surfaces do not have boundaries 
%             and thus flag=0 for all k and m.
%   - H     : level set handles.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Basic error checking
if nargin<2 || isempty(TR) || isempty(F)
   error('Insufficient number of input arguments') 
end

[Tri,X,fmt]=GetMeshData(TR);
if fmt>1, TR=triangulation(Tri,X); end

F=F(:);
if ~isnumeric(F) || ~ismatrix(F) || numel(F)~=size(X,1) || sum(isnan(F) | isinf(F))>0
    %error('Invalid entry for 2nd input argument (F)')
end

if nargin<3 || isempty(iv)
    iv=11;
elseif ~isnumeric(iv) || ~isvector(iv) || (numel(iv)==1 && (iv<=0 || iv~=round(iv) || isnan(iv) || iv>5001))
    error('Invalid entry for 3rd input argument (iv)')
else
    iv=iv(:);
end

if nargin<4 || isempty(vis)
    vis=false;
elseif numel(vis)~=1 || ~((ishandle(vis) && strcmpi(get(vis,'type'),'axes')) || islogical(vis))
    error('Invalid entry for 4th input argument (vis)')
end

% Level-set values
F_min=min(F);
F_max=max(F);
if numel(iv)==1   
    if iv==1
        iv=(F_max+F_min)/2;
    else
        dF=(F_max-F_min)/iv;
        iv=linspace(F_min+dF/2,F_max-dF/2,iv);
    end
elseif numel(iv)==2 && iv(1)==iv(2)
    iv=iv(1);
end
    
% Edges 
E=[Tri(:,[1 2]);Tri(:,[2 3]);Tri(:,[3 1])];

X1=X(Tri(:,1),:);
X2=X(Tri(:,2),:);
X3=X(Tri(:,3),:);

F1=F(E(:,1)); F1=reshape(F1,[],3);
F2=F(E(:,2)); F2=reshape(F2,[],3);

Fe_min=min(F1,F2);
Fe_max=max(F1,F2);

n=numel(iv);
chk=false(size(X,1),1);
for i=1:n, chk=chk | abs(F-iv(i))<2*eps; end

chk_edges=false;
if sum(chk)>1
    Eo=sort(edges(TR),2);
    [Fe1,Fe2]=deal(F(Eo(:,1)),F(Eo(:,2)));
    id_eq=abs(Fe1-Fe2)<=2*eps; % edges whose verices have equal field values
    if sum(id_eq)>0
        chk_edges=true;
        Eo_eq=Eo(id_eq,:);
        e_id_val=~ismember(sort(E,2),Eo_eq,'rows');
        e_id_val=reshape(e_id_val,[],3);
    end
end


% Begin visualization
ha=NaN;
if islogical(vis) && vis    
    figure('color','w')
    axis equal
    hold on
    h=trimesh(TR);
    set(h,'EdgeColor','none','FaceColor',0.75*[1 1 1],'FaceAlpha',0.75,...
          'SpecularExponent',100,'SpecularStrength',0.25);
    ha=gca;
    h1=camlight('headlight');
    set(h1,'style','infinite','position',10*get(h1,'position'))
    h2=light('position',-get(h1,'position'));
    set(h2,'style','infinite')
    lighting phong
elseif ishandle(vis)
    ha=vis;
    hold on
end

% Does the input mesh have a boundary?
FB=freeBoundary(TR);
if ~isempty(FB), FB=sort(FB,2); end

if nargout>2 || ~isempty(FB) || chk_edges
    flag=true;
    VT=cell(n,1);
else
    flag=false;    
end

% Level-sets
Q=cell(n,1);
H=nan(1,n);
for i=1:n % loop through level sets

    % Edges and triangles traversed by the i-th level set
    e_id=Fe_min<=iv(i) & Fe_max>=iv(i);    
    if chk_edges, e_id=e_id & e_id_val; end
    
    m=sum(e_id,2);
    f_id=m>0; % triangles traversed by the level set    
    if ~any(f_id)
        fpritnf(2,'Level set %10e does not exist\n',iv(i))
        continue
    end
    e_id=e_id(f_id,:);    
    
    %                         E1       E2       E3
    % Triangle vertices : X1 ----> X2 ----> X3 ----> X1
    X1_i=X1(f_id,:);
    X2_i=X2(f_id,:);
    X3_i=X3(f_id,:); 
    v_id=Tri(f_id,:);
     
    % Points of intersection, where level set crosses the edges
    dFe=F2(f_id,:)-F1(f_id,:);
    t=(iv(i)-F1(f_id,:))./dFe;    

    id_slf=abs(dFe)<=eps;
    if sum(id_slf(:))>0 % check if there are edges coincident with the level set
        t(id_slf)=0;
    end
    
    P1=bsxfun(@times,1-t(:,1),X1_i) + bsxfun(@times,t(:,1),X2_i); % E1 
    P2=bsxfun(@times,1-t(:,2),X2_i) + bsxfun(@times,t(:,2),X3_i); % E2
    P3=bsxfun(@times,1-t(:,3),X3_i) + bsxfun(@times,t(:,3),X1_i); % E3
    
    % Line segments comprising the level set
    % ---------------------------------------------------------------------
    
    % Three possible options are
    id_12=e_id(:,1) & e_id(:,2) & sqrt(sum((P1-P2).^2,2))>eps; % E1--E2
    id_13=e_id(:,1) & e_id(:,3) & sqrt(sum((P1-P3).^2,2))>eps; % E1--E3
    id_23=e_id(:,2) & e_id(:,3) & sqrt(sum((P2-P3).^2,2))>eps; % E2--E3
        
    P12=[]; V12=[]; 
    k=sum(id_12);
    if k>0

        if flag            
            [V1,V2]=deal([v_id(id_12,1),v_id(id_12,2) t(id_12,1)],[v_id(id_12,2),v_id(id_12,3) t(id_12,2)]);
            V12=cat(3,V1,V2);
            V12(:,:,3)=NaN;
            V12=permute(V12,[3 1 2]);
            V12=reshape(V12,[],3);
        end
        
        [p1,p2]=deal(P1(id_12,:),P2(id_12,:));
        
        P12=cat(3,p1,p2);
        P12(:,:,3)=NaN;
        P12=permute(P12,[3 1 2]);
        P12=reshape(P12,[],3);
        
        id_13=id_13 & ~id_12;
        id_23=id_23 & ~id_12;
        
    end
    
    P13=[]; V13=[];
    k=sum(id_13);
    if k>0 
        
        if flag
            [V1,V3]=deal([v_id(id_13,1),v_id(id_13,2) t(id_13,1)],[v_id(id_13,3),v_id(id_13,1) t(id_13,3)]);
            V13=cat(3,V1,V3);
            V13(:,:,3)=NaN;
            V13=permute(V13,[3 1 2]);
            V13=reshape(V13,[],3);
        end

        [p1,p3]=deal(P1(id_13,:),P3(id_13,:));
        P13=cat(3,p1,p3);
        P13(:,:,3)=NaN;
        P13=permute(P13,[3 1 2]);
        P13=reshape(P13,[],3);
        
        id_23=id_23 & ~id_13;
        
    end
    
    P23=[]; V23=[];
    k=sum(id_23);
    if k>0                
        
        if flag
            [V2,V3]=deal([v_id(id_23,2),v_id(id_23,3) t(id_23,2)],[v_id(id_23,3),v_id(id_23,1) t(id_23,3)]);
            V23=cat(3,V2,V3);
            V23(:,:,3)=NaN;
            V23=permute(V23,[3 1 2]);
            V23=reshape(V23,[],3);
        end
        
        [p2,p3]=deal(P2(id_23,:),P3(id_23,:));
        P23=cat(3,p2,p3);
        P23(:,:,3)=NaN;
        P23=permute(P23,[3 1 2]);
        P23=reshape(P23,[],3);
        
    end    
    P=cat(1,P12,P13,P23);  
    VTi=cat(1,V12,V13,V23);
    
    if chk_edges
        
        e_id=abs(Fe1-iv(i))<=2*eps & id_eq;
        chk_match=sum(e_id)>0;
        if chk_match
            
            % Revome duplicate edges
            [t1,t2]=deal(VTi(1:3:end,3),VTi(2:3:end,3));
            idx=(t1<=eps | t1>=(1-eps)) & (t2<=eps | t2>=(1-eps));
            if sum(idx)>0 
                [v1,v2]=deal(VTi(1:3:end,1),VTi(2:3:end,2));
                [v1,v2]=deal(v1(idx),v2(idx));
                chk=ismember(sort([v1 v2],2),Eo(e_id,:),'rows');
                if sum(chk)>0
                    idx=find(idx);
                    idx=idx(chk);
                    
                    P=reshape(P',3,3,[]);
                    P(:,:,idx)=[];
                    P=reshape(P,3,[])';
                    
                    VTi=reshape(VTi',3,3,[]);
                    VTi(:,:,idx)=[];
                    VTi=reshape(VTi,3,[])';
                end
            end
            
            % Edges coincident with the level set
            [x1,x2]=deal(X(Eo(e_id,1),:),X(Eo(e_id,2),:));            
            Pe=cat(3,x1,x2);
            Pe(:,:,3)=NaN;
            Pe=permute(Pe,[3 1 2]);
            Pe=reshape(Pe,[],3);
            P=cat(1,P,Pe);
            
        end        
    else
        chk_match=false;
    end    
    Q{i}=P;
    
    if flag
        if chk_edges && chk_match
            Ve=repmat(Eo(e_id,:),[1 1 2]);
            Ve(:,3,:)=0;
            Ve(:,3,2)=1;            
            Ve(:,:,3)=NaN;
            Ve=permute(Ve,[3 1 2]);
            Ve=reshape(Ve,[],3);
            VTi=cat(1,VTi,Ve);
        end
        VTi(:,4)=0;
        VT{i}=VTi;        
    end
    
    
    % Move vertices on the boundary edges to the top; so it is easier
    % to extract continuous polylines during post-processing with the
    % 'OrderIsoContourVerts.m' function
    if ~isempty(FB) && nargout>2
        
        bv_id=ismember(sort(VTi(:,1:2),2),FB,'rows');
        if chk_match
            idx=ismember(VTi(:,1:2),FB(:));
            idx=(idx(:,1) & VTi(:,3)<=eps) | (idx(:,2) & VTi(:,3)>=(1-eps));
            bv_id=bv_id | idx;
        end
        
        if sum(bv_id)>0
            
            idx=unique(floor(find(bv_id)/3))+1;
            
            ns=size(P,1)/3;
            P=mat2cell(P,3*ones(ns,1),3);
            VTi=mat2cell(VTi,3*ones(ns,1),4);
            
            bv_id=reshape(bv_id,3,[]);            
            for k=1:numel(idx)
                j=idx(k);
                if bv_id(2,j) && ~bv_id(1,j)
                    P{j}=P{j}([2 1 3],:);
                    VTi{j}=VTi{j}([2 1 3],:);
                end
                VTi{j}(1,4)=1;
            end
            
            P_b=P; 
            P_b(idx)=[];
            P=cell2mat(cat(1,P(idx),P_b));
            Q{i}=P;
            
            VTi_b=VTi;
            VTi_b(idx)=[];
            VT{i}=cell2mat(cat(1,VTi(idx),VTi_b));
                        
        end
    end

    
    if ishandle(ha) && ~isempty(P)
        H(i)=plot3(ha,P(:,1),P(:,2),P(:,3),'-k','LineWidth',2);
    end
    
    if isempty(P)
        fprintf(2,'Value %10e is not in the the domain of F\n',iv(i));
    end
    
end

