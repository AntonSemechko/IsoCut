function [TRc,LF] = InsertIsoCut(TR,C,FV)
% Locally modify connectivity of a triangular surface mesh so that it 
% contains edges coincident with an iso-contour computed with the 
% 'IsoContour' function or another function that immitates its output.
%
% INPUT
%   - TR    : input surface mesh represented as an object of 'TriRep' 
%             class, 'triangulation' class, or a cell such that TR={Tri,X},
%             where Tri is an M-by-3 array of faces and X is an N-by-3 
%             array of vertex coordinates.  
%   - C     : N-by-3 ORDERED list of coordinates of the cut you wish to
%             make in the input mesh. 
%   - FV    : N-by-4 list of barycentric coordinates corresponding to C.
%             Note, both C and FV are the outputs of the 'OrderIsoContourVerts'
%             function; which pre-processes the iso-contour(s) genenerated 
%             by the 'IsoContour' function to ensure sequential vertex order.  
%
% OUTPUT
%   - TRc   : input surface mesh whose connectivity has been locally 
%             modified to insert the cut specified by C and FV.     
%   - LF    : map between the faces of TR and TRc.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Basic error checking
if nargin<3 || isempty(TR) || isempty(C) || isempty(FV)
    error('Insufficent number of input arguments')
end

if ~isnumeric(C) || ~ismatrix(C) || size(C,2)~=3 || size(C,1)<2
    error('Invalid entry for 2nd inpupt argument (C)')
end

if ~isnumeric(C) || ~ismatrix(FV) || size(FV,2)~=4 || size(FV,1)~=size(C,1)
    error('Invalid entry for 3rd inpupt argument (FV)')
end
    

% Mesh data
[Tri,X,fmt] = GetMeshData(TR);
clear TR


% Identify and remove repeating vertices in the cut
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tol = 1E-12;
Nc = size(C,1);

v_id = zeros(Nc,1);
t = FV(:,3);

idx = t<tol; 
v_id(idx) = FV(idx,1);

idx = t>(1-tol);
v_id(idx) = FV(idx,2);

idx = false(Nc,1);
for i = 1:(Nc-3)    
    
    if idx(i) || v_id(i)==0, continue; end   
    
    if v_id(i)==v_id(i+1) % look ahead by one
        idx(i+1) = true;
    elseif v_id(i)==v_id(i+2) % look ahead by two
        idx(i+2) = true;
    end    
end

if any(idx)
    warning('Detected %d repeating vertices. Need to improve logic for extracting sequential point samples.',nnz(idx))
    C(idx,:) = [];
    FV(idx,:) = [];
    Nc = size(C,1);
end


% Deal with open cuts. If they terminate at the boundary or existing 
% vertices, great, if not, have to add extra points before and after the 
% cut to ensure correct triangulation.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chk_contr_closed = norm(C(1,:)-C(end,:))<tol;
if ~chk_contr_closed    

    
    id = [1 2; Nc Nc-1]; % end point indices into FV
    for j = 1:2
        
        fv1 = FV(id(j,1),:); % [v1 v2 t flag]
        
        chk_vtx = fv1(3)<=tol | fv1(3)>=(1-tol);
        if chk_vtx || fv1(4)==1, continue; end % boundary or existing vertex?

        % Opposing vertex
        fv2 = FV(id(j,2),:);         
        
        tri = Tri(sum(ismember(Tri,fv1(1:2)),2)==2,:); % faces attached to fv1(1:2)
        tri(sum(ismember(tri,fv2(1:2)),2)==2,:) = [];  % remove face that will be cut
        
        v = tri(~ismember(tri,fv1(1:2)));            % vertex that opposes fv1

        % Insert opposing vertex into C and FV
        tri = circshift(tri,1-find(tri==v)); 
        if j==1 % insert before
            C = cat(1, X(v,:), C);
            FV = cat(1, [tri(1:2) 0 0], FV);            
            id = id(2,:) + 1;
        else % insert after
            C = cat(1, C, X(v,:));
            FV = cat(1, FV, [tri(1:2) 0 0]);
        end
        Nc = Nc + 1;
                
    end   
       
    chk_contr_closed = norm(C(1,:) - C(end,:))<=eps;
end


% Assign indices to contour vertices. This automatically accounts for 
% unprobable, but possible, cases where contour passes though existing 
% vertices.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if chk_contr_closed
    Nc = Nc - 1;  % # of unique vertices
    Ne = Nc;      % max # of new edges that will be inserted 
else    
    Ne = Nc - 1;  % # of edges that will be inserted
end
v_id = zeros(Nc + 1,1);

Xc = zeros(0,3);  % unique list of contour vertices that has NULL intersection with X 

k = 0;
Nx = size(X,1);
for i = 1:Nc
    t = FV(i,3);
    if t<=tol
        v_id(i) = FV(i,1);
    elseif t>=(1-tol)
        v_id(i) = FV(i,2);
    else
        k = k+1;
        v_id(i) = Nx + k;
        Xc(k,:) = C(i,:);
    end
end
v_id(Nc + 1) = v_id(1); % indices assigned to points in C

% Update vertex list
X = cat(1,X,Xc);


% Modify connectivity
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Tri_new = zeros(0,3);  f_tri_new = zeros(0,1);
Quad_new = zeros(0,4); f_quad_new = zeros(0,1);

[cnt_tri,cnt_quad] = deal(0);
for n = 1:Ne

    % Points containing the line segment 
    fv1 = FV(n,:);    
    fv2 = FV(n+1,:);   
    
    % Case1: Cut passes through an existing edge so connectivity doesn't 
    % have to be modified
    if (fv1(3)<=tol || fv1(3)>=(1-tol)) && (fv2(3)<=tol || fv2(3)>=(1-tol))
        continue        
    end

    % Case 2: Cut passes through an existing vertex --> Triangle gets split
    % into two smaller triangles
    chk = true;
    if fv1(3)<=tol
        v = fv1(1);        % 1st point of the cut passes thought this vertex
        fv = fv2;          % 2nd point of the cut passeses though here 
        v_new = v_id(n+1); % index of the 2nd point
    elseif fv1(3)>=(1-tol)
        v = fv1(2);
        fv = fv2;
        v_new = v_id(n+1);
    elseif fv2(3)<=tol
        v = fv2(1);
        fv = fv1;
        v_new = v_id(n);
    elseif fv2(3)>=(1-tol)
        v = fv2(2);
        fv = fv1;
        v_new = v_id(n);
    else
        chk = false; % must be Case 3 
    end
    
    if chk
        % Modify connectivy of the cut triangle. Modification is reflected
        % in both Tri and Tri_new. Tri_new contains the 2nd triangle
        % generated by the cut.        
        cnt_tri = cnt_tri + 1;   
        [Tri,Tri_new(cnt_tri,:),f_tri_new(cnt_tri,1)] = split_face_1(v,fv,Tri,v_new);
        continue
    end
    
    % Case 3: Cut passes through two distinct edges. Thus triangle gets 
    % split into one smaller triangle and one quadrilateral. This is by far 
    % the most common case. The new, smaller triangle gets saved into Tri
    % and quadrilateral into Quad_new.
    cnt_quad = cnt_quad + 1;
    [Tri,Quad_new(cnt_quad,:),f_quad_new(cnt_quad,1)] = split_face_2(fv1,fv2,Tri,[v_id(n) v_id(n+1)]);
    
end

% Split quads into triangles
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if ~isempty(Quad_new)
    TriQuad = Quad2Tri(Quad_new,X);
    Tri_new = cat(1,Tri_new,TriQuad);
    
    f_quad_new = [f_quad_new f_quad_new]';
    f_tri_new = cat(1,f_tri_new(:),f_quad_new(:));    
end


% Final mesh
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Nf = size(Tri,1);
LF = [(1:Nf)'; f_tri_new(:)];

switch fmt
    case 1
        TRc = triangulation(cat(1,Tri,Tri_new),X);
    case 2
        TRc = TriRep(cat(1,Tri,Tri_new),X); %#ok<*DTRIREP>
    case 3
        TRc = {cat(1,Tri,Tri_new) X};
    case 4
        TRc = struct('faces',cat(1,Tri,Tri_new),'vertices',X);
end



function [Tri,quad_new,f_mod] = split_face_2(fv1,fv2,Tri,v_new)
% Split triangle into one smaller triangle and one quadrilateral using a 
% line segment that cuts through two of its edges.
%
%   - fv1, fv2  : edges containg end-points of the line segment being
%                 insterted into the mesh 
%   - Tri       : mesh connectivity
%   - v_new     : 1-by-2 vector of vertex indices assigned to fv1 and fv2 


fv1 = fv1(1:2);
fv2 = fv2(1:2);
fv = [fv1 fv2];

% Identify triangle
idx = ismember(Tri,fv);
idx = sum(idx,2)==3;
%idx = find(idx);
f_mod = find(idx);

% Put vertex common to fv1 and fv2 at the top
tri = Tri(f_mod,:);
if nnz(fv==tri(2))==2
    tri = circshift(tri,[0 -1]);
elseif nnz(fv==tri(3))==2
    tri = circshift(tri,[0 -2]);
end

if nnz(ismember(tri([1 3]),fv1))==2
    v_new = v_new([2 1]);
end

% Modify connectivity
Tri(f_mod,:) = [tri(1) v_new];
quad_new = [v_new(1) tri(2:3) v_new(2)];


function [Tri,tri_new,f_mod] = split_face_1(v,fv,Tri,v_new)
% Split triangle into two smaller triangles using a line segment that 
% starts at vertex v and terminates at the opposing edge.
%
%   - v     : existing vertex coinsident with the cut
%   - fv    : opposing edge
%   - Tri   : mesh connectivity
%   - v_new : index of the new vertex on the opposing edge
  

% Identify triangle containing the cut
idx = ismember(Tri,[v fv(1:2)]);
idx = sum(idx,2)==3;
%idx=find(idx);
f_mod = find(idx);

% Put v at the top
tri = Tri(f_mod,:);
id = find(tri==v);
tri = circshift(tri,[0 1-id]);

% Modify connectivity
Tri(f_mod,:) = [tri([1 2]) v_new];
tri_new = [tri([3 1]) v_new];

