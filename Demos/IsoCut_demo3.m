function IsoCut_demo3(iv)
% Modify mesh connectivity so mesh contains edges conicident with a cut.
%
% INPUT:
%   - k     : (optional) level set value. Must be a real number between 
%             0.01 and 0.99. iv=0.60 is the default setting.
%             
% AUTHOR: Anton Semechko (a.semechko@gmail.com)


%% Check the inputs
% -------------------------------------------------------------------------

if nargin<1 || isempty(iv)
    iv = 0.6;
elseif ~isnumeric(iv) || ~isscalar(iv) || iv<0.01 || iv>0.99 
    iv = 0.6;
    fprintf(2,"Invalid entry for 'iv'. Using default setting iv=0.6.")
end


%% Step 1: Load mesh
TR = load('sample_isocut_meshes.mat','bunny');
TR = TR.bunny;


%% Step 2: Define cutting plane normal
P = [ 0.034488  0.036484  0.006912; ...
      0.019104  0.055041  0.047894; ...
     -0.008596  0.123650  0.033786];

N = cross(P(1,:) - P(2,:), P(3,:) - P(2,:));
N = N/norm(N);


%% Step 3: Scalar field corresponding to the equation of the plane

X = TR.vertices;

F = bsxfun(@minus,X,P(1,:))*N(:);

Fmin = min(F);
Fmax = max(F);

F = F-Fmin;
F = F/(Fmax-Fmin);


hf = figure('color','w');
maximize_fig(hf);


ha1 = subtightplot(1,3,1);
h = patch(TR); 
set(h,'FaceColor',0.5*[1 1 1],'EdgeColor','k','FaceAlpha',0.95,'EdgeAlpha',1)
axis equal
hold on

[hm1,~,hb1] = VisualizeScalarFieldOnTriMesh(TR,F,ha1); % visualize mesh and scalar field

set(hm1,'FaceAlpha',0.9) 

vis_cutting_plane(X,N,iv,ha1)


%% Step 4: Insert cut into the mesh
[TRc,Fc,C,Hc] = IsoCut(TR,F,iv,ha1); 
set(Hc,'color','r','LineWidth',1)

% Visulize modified mesh
ha2 = subtightplot(1,3,2);
h = patch(TRc); 
set(h,'FaceColor',0.5*[1 1 1],'EdgeColor','k','FaceAlpha',0.95,'EdgeAlpha',1)
axis equal
hold on

[hm2,~,hb2] = VisualizeScalarFieldOnTriMesh(TRc,Fc,ha2); % visualize mesh and scalar field
set(hm2,'FaceAlpha',0.9) 



%% Step 5: Retain portion of the cut mesh with F<=iv


[Tri,X] = GetMeshData(TRc);

chk_v = Fc<=iv;                 % check is F(i)<=iv for each vertex
chk_f = sum(chk_v(Tri),2)==3;  % faces whos vertices have F<=iv

Tri2 = Tri(chk_f,:);           % list of face corresponding to chk_f
TR2 = RemoveNonRefVerts({Tri2 X});

TR2 = struct('faces',TR2{1},'vertices',TR2{2});
F2 = Fc(chk_v);

% Visulize retained portion of the mesh  
ha3 = subtightplot(1,3,3);
h = patch(TR2); 
set(h,'FaceColor',0.5*[1 1 1],'EdgeColor','k','FaceAlpha',0.95,'EdgeAlpha',1)
axis equal
hold on

[hm3,~,hb3] = VisualizeScalarFieldOnTriMesh(TR2,F2,ha3); % visualize mesh and scalar field
set(hm3,'FaceAlpha',0.9) 

for i = 1:numel(C)
    plot3(ha2,C{i}(:,1),C{i}(:,2),C{i}(:,3),'-r','LineWidth',1)
    plot3(ha3,C{i}(:,1),C{i}(:,2),C{i}(:,3),'-r','LineWidth',2)
end

colormap(hf,'parula')
delete([hb1 hb2 hb3]) % delete colorbars
drawnow
pause(0.1)
set([ha1 ha2 ha3],'CLim',[0 1]) % enfore consistent color maps


% Sync the axes
hLink12 = MatchAxesView(ha1,ha2);
hLink13 = MatchAxesView(ha1,ha3);
hf.UserData = [hLink12 hLink13];



function vis_cutting_plane(X,N,iv,ha)

N = N(:);


Xo = mean(X,1);
X = bsxfun(@minus,X,Xo);       % shift centroid to origin
x = X*N;                       % distance to the plane
dX = X - bsxfun(@times,x,N');  % subspace orthogonal to N

% Principal vectors ortognal to cutting plane
[V,~] = svd(dX'*dX);

% Vertices of quadrilateral used to represent the cutting plane
dY = dX*V;
dYmin = min(dY);
dYmax = max(dY);

P1 = dYmin;
P2 = [dYmax(1) dYmin(2) 0];
P3 = [dYmax(1) dYmax(2) 0];
P4 = [dYmin(1) dYmax(2) 0];
P = 1.15*[P1;P2;P3;P4];

% Transform vertices back into the orignal frame 
x_min = min(x);
x_max = max(x);
P = bsxfun(@plus,P*V',Xo);
P = bsxfun(@plus,P,x_min*N' + iv*(x_max-x_min)*N'); 

% Visualize
h = patch(ha,'faces',[1 2 3 4],'vertices',P);
set(h,'FaceColor','r','FaceAlpha',0.5,'EdgeColor','k','LineWidth',1)

