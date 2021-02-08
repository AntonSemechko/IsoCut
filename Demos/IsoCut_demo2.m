function IsoCut_demo2(k,a)
% Example of how to extract coordinates of non-planar iso-contours.
%
% INPUT:
%   - k     : (optional) number of level sets. k=1 is the default 
%             setting. For the purposes of the demo 1=<k<=20.
%   - a     : (optional) Lame curve parameter. a=0.9 is the default
%             setting.
%             
% AUTHOR: Anton Semechko (a.semechko@gmail.com)


%% Check the inputs
% -------------------------------------------------------------------------

if nargin<1 || isempty(k)
    k=1;
elseif ~isnumeric(k) || ~isscalar(k) || k<1 || k>20 || k~=round(k)
    k=1;
    fprintf(2,"Invalid entry for 'k'. Using default setting k=1.")
end

if nargin<2 || isempty(a)
    a=0.9;
elseif ~isnumeric(a) || ~isscalar(a) 
    a=0.9;
    fprintf(2,"Invalid entry for 'a'. Using default setting a=0.9.")
end


%% Step 1: Load mesh
TR=load('sample_meshes.mat','sphere');
TR=TR.sphere;


%% Step 2: Define scalar field at the mesh vertices

% Stereo projection through y=1
X=TR.vertices;
Y=X(:,[1 3]);
D=1+X(:,2);

P=bsxfun(@rdivide,Y,D+eps); 

% Super ellipse
F=abs(P(:,1)).^a + abs(P(:,2)).^a - 1;

% Constrain vales for vusualization purposes
f=2.5;
F(F>f)=f;
F(F<-f)=-f;
%F=F-f;

% Deal with point at infinity
F(D<eps)=f;



%% Step 3: Visualize field on mesh
[hm,ha,hb,hl]=VisualizeScalarFieldOnTriMesh(TR,F); % visualize mesh and scalar field
colormap(ha,'parula')
view([-20 0])
try
    hl(1).Position=[0 -100 0];
    hl(2).Position=[0 100 0];
catch
end
set(hm,'FaceAlpha',0.9)


%% Step 4: Visualize level sets
if k==1 % cut along zero level set
    ls=f-0.4;
    [Q,iv,VT]=IsoContour(TR,F,ls*[1 1],ha); % information about the cut is stored in the Q and VT variables
else         
    iv=linspace(fix(10*min(F))/10 + 0.1, f-0.4, 20);
    if k<20
        iv(1:(20-k))=[];
    end
    [Q,iv,VT]=IsoContour(TR,F,iv,ha);
end


%% Step 5: Make sure vertices making up the iso-contour are ordered sequentially.
% This step is optional, and only necessary if you are using the iso-contours for
% something other than visualization.
 [C,VT]=OrderIsoContourVerts(Q,VT); 

 
