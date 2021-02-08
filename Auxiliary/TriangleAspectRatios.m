function [AR,TA]=TriangleAspectRatios(TR,ARdef)
% Compute aspect ratios (AR) of surface mesh triangles.
%
% INPUT :
%   - TR    : surface mesh represented as an object of 'TriRep' class,
%             'triangulation' class, or a cell such that TR={Tri,V}, where
%             Tri is an M-by-3 array of faces and V is an N-by-3 array of 
%             vertex coordinates.
%   - ARdef : integer (1 or 2) used to specify preferred definition of the
%             triangle AR:
%             1. AR = length of the longest edge divided by the length of 
%                     the shortest altitude 
%             2. AR = 2*Rin/Rout, where Rin and Rout are the triangle 
%                     inradius and circumradius, respectively. <--- {DEFAULT}
%
% OUTPUT:
%   - AR    : M-by-1 array of triangle ARs, where M is the total number of
%             mesh triangles. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Face and vertex lists
[Tri,V]=GetMeshData(TR);
if size(Tri,2)~=3
    error('This function is intended ONLY for triangular surface meshes')
end

% Triangle aspect ratio definition
if nargin<2 || isempty(ARdef), ARdef=2; end

% Edge lengths
V1=V(Tri(:,1),:);
V2=V(Tri(:,2),:);
V3=V(Tri(:,3),:);

D12=V2-V1;
D23=V3-V2;
D31=V1-V3;

E1=sqrt(sum(D12.^2,2));
E2=sqrt(sum(D23.^2,2));
E3=sqrt(sum(D31.^2,2));

% Triangle areas
TA=cross(D12,-D31,2);
TA=sqrt(sum(TA.^2,2))/2;

% Compute triangle aspect ratios according to the specified definition
switch ARdef
        
    case 1 % 1st (default) definition
        
        % Which of the three edges is the longest?
        Emax=max([E1 E2 E3],[],2);
        
        % Aspect ratios
        AR=(Emax.^2)./(2*TA);
        
    case 2 % 2nd definition
            
        %Rin=2*TA./(E1+E2+E3);
        %Rout=(E1.*E2.*E3)./(4*TA);
        %AR=2*Rin./Rout;        
        
        AR=16*(TA.^2)./((E1+E2+E3).*(E1.*E2.*E3));        
        AR(isnan(AR))=0;
        
    otherwise
    
    error('Invalid format for 2nd input argument (ARdef)')

end

