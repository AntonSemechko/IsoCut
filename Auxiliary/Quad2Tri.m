function Tri = Quad2Tri(F,V)
% Convert a quadrilateral surface mesh into a triangular surface mesh so
% that the minimum aspect ratio of the resulting triangles is maximized. 
%
% INPUT:
%   - F     : M-by-4 array of quad face-vertex connectivities
%   - V     : N-by-3 array of mesh vertex coordinates
%
% OUTPUT:
%   - Tri   : (2*M)-by-3 array of triangle face-vertex connectivities
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


flag = false;
if nargin==1
    if isa(F,'double')
        error('Insufficient number of input arguments')
    else
        flag = true;
        try
            [F,V,fmt] = GetMeshData(F);
        catch
            error('Unrecognized input mesh format')
        end
    end
end

% Quad split pattern A
F1a = F(:,[1 2 3]);
F2a = F(:,[3 4 1]);

% Quad split pattern B
F1b = F(:,[1 2 4]);
F2b = F(:,[2 3 4]);

% For every quad, select split pattern than maximizes minimum aspect ratio
if size(V,2)==2, V(:,3) = 0; end

AR1a = TriangleAspectRatios({F1a V});
AR2a = TriangleAspectRatios({F2a V});
AR1 = min(AR1a,AR2a);

AR1b = TriangleAspectRatios({F1b V});
AR2b = TriangleAspectRatios({F2b V});
AR2 = min(AR1b,AR2b);

Tri = cat(1,F1a',F2a');
idx = (AR1<AR2)';

if any(idx)
    Tri_b = cat(1,F1b',F2b');
    Tri(:,idx) = Tri_b(:,idx);
end
Tri = reshape(Tri,3,[])';

if flag
    Tri = MeshOut(Tri,V,fmt);
end


function TR = MeshOut(Tri,X,fmt)

switch fmt
    case 1
        TR = triangulation(Tri,X);
    case 2
        TR = TriRep(Tri,X); %#ok<*DTRIREP>
    case 3
        TR = {Tri X};
    case 4
        TR = struct('faces',Tri,'vertices',X);
end

