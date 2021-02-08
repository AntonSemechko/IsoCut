function IsoCut_demo1(n,k)
% Example of how to extract coordinates of planar iso-contours 
% corresponding to the intersection of a surface mesh with multiple, 
% uniformly-spaced, cutting planes.
%
% INPUT:
%   - n     : (optional) integer in the range 1 to 3, inclusive, specifying
%             prinicpal axis along which to cut the mesh. n=1 is the
%             default setting.
%   - k     : (optional) number of cutting planes. k=20 is the default 
%             setting. For the purposes of the demo 1=<k<=100.
%             
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


%% STEPS:
% 1. Load sample mesh
% 2. Define direction along which to slice the mesh. Cutting planes will be
%    perpendicular to this direction.
% 3. Define scalar field at the mesh vertices. This is done by evaluating 
%    equation of the plane. 
% 4. Extract line segments compirising the intersection contours using the
%    'IsoContour' function. Contours are level sets of the scalar field 
%     from step 3.  
% 5. Stitch togther the line segments of each countour into ordered vertex
%    sequences using the 'OrderIsoContourVerts' function. This step is
%    optional. 


%% NOTES:
% 1. The maximum number of cutting planes is 5000, but for the pupose of 
%    the demo, the upper limit has been set to 100.  
% 2. Planes do not have to be spaced uniformly . See documentation for the 
%    'IsoContour' function for additional information.
% 4. A given level set will produce at least one disconnected iso-contour,
%    but there may be several due to non-convex geometry of the surface.  
% 3. To get contour(s) corresponding to the intersection of mesh with a 
%    single cutting plane, obtain sclar field from the equation of the plane
%       F = X*V(:) + P(:)' , where
%           X : N-by-3 array of mesh vertex cooridnates
%           V : plane normal
%           P : point on the plane
%    The cut occurs along the zero-level set of F. 


%% Check the inputs
% -------------------------------------------------------------------------
if nargin<1 || isempty(n)
    n=1;
elseif ~isnumeric(n) || ~isscalar(n) || ~ismember(n,1:3)
    n=1;
    fprintf(2,"Invalid entry for 'n'. Using default setting n=1.")
end

if nargin<2 || isempty(k)
    k=20;
elseif ~isnumeric(k) || ~isscalar(k) || k<1 || k>100 || k~=round(k)
    k=20;
    fprintf(2,"Invalid entry for 'k'. Using default setting k=20.")
end


%% Step 1: Load mesh
TR=load('sample_meshes.mat','tooth');
TR=TR.tooth;


%% Step 2: Cut-plane normal. Use principal axis of covariance matrix for demo
X=TR.vertices;
Xo=mean(X);
dX=bsxfun(@minus,X,Xo);
C=dX'*dX;
[R,~]=svd(C); % principal axes are stacked along columns of R


%% Step 3: Define scalar field. Mesh mesh will be cut along level sets of
% this field.
F=dX*R(:,n); % project vertices on the n-th pricipal axis

F=F-min(F);
[hm,ha,hb,hl]=VisualizeScalarFieldOnTriMesh(TR,F); % visualize mesh and scalar field
colormap(ha,'parula')


%% Step 4: Cut the mesh using k uniformly spaced planes
[Q,iv,VT]=IsoContour(TR,F,k,ha); % information about the iso-contours is stored in the Q and VT variables


%% Step 5: Make sure vertices making up the iso-contours are ordered sequentially.
% This step is optional, and only necessary if you are using the iso-contours for
% something other than visualization.
 [C,VT_out]=OrderIsoContourVerts(Q,VT); % coordinates of the contours are in the cell C. To associate contours with specific iso-value/level-set, call OrderIsoContourVerts(Q{i},VT{i}) separately for each iso-value iv(i), 1<=i<=k 

 
