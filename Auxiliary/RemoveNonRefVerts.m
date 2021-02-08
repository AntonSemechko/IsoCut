function [TR,idx_v]=RemoveNonRefVerts(TR)
% Remove non-referenced mesh vertices.
%
% INPUT:
%   - TR    : input surface mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%
% OUTPUT:
%   - TR    : output mesh; same format as the input mesh.
%   - idx_v : list of referenced vertices.
%   
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,V,fmt]=GetMeshData(TR);
d=size(Tri,2);

[idx_v,~,idx]=unique(Tri(:));
V=V(idx_v,:);

Tri=(1:size(V,1))';
Tri=Tri(idx);
Tri=reshape(Tri,[],d);

switch fmt
    case 1
        TR=triangulation(Tri,V);
    case 2
        TR=TriRep(Tri,V); %#ok<*DTRIREP>
    case 3
        TR={Tri V};
    case 4
        TR=struct('faces',Tri,'vertices',V);
end

