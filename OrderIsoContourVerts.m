function [C,VT_out] = OrderIsoContourVerts(Q,VT,tol)
% Get continuous polyline contours from unordered level set line segments
% generated by the 'IsoContour' function. 
%
% INPUT:
%   - Q      : K-by-1 cell containing unordered line segments comprising K 
%			   level-sets.
%   - VT     : K-by-1 cell generated by the 'IsoContour' function; contains
%              info used to keep track of the edges intersected by
%              the level sets. Set VT=[] if this info is unavailable or 
%              you wish to omit it (NOT RECOMMENDED).
%   - tol    : (optional) smallest allowable (Euclidean) distance between
%              two vertices; tol=1E-11 is the default setting.
%
% OUTPUT:
%   - C      : N-by-1 cell containing ordered lists of vertices comprising
%              the contours in Q, where N>=K.  
%   - VT_out : correspoding version of VT after post-processing.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


% Basic error checking
if nargin<1 || isempty(Q) 
    error('Insufficient number of input arguments')
elseif ~iscell(Q)
    error('Invalid format for 1st input argument (Q)')
end

if nargin<2 || isempty(VT)
    VT = [];
    chk_vt = false;
elseif ~isequal(numel(Q),numel(VT)) || ~iscell(VT)
    error('Invalid entry for 2nd input argument (VT)')
else
    chk_vt = true;
end

if nargin<3 || isempty(tol)
    tol = 1E-11;
elseif ~isnumeric(tol) || numel(tol)~=1 || tol<eps
    error('Invalid entry for 3rd input argument (tol)')
end
tol = max(tol,1E-11);


% Loop through level-sets
[C,VT_out] = deal(cell(1));
N = 0;
for k = 1:numel(Q)
    
    if isempty(Q{k}), continue; end
    
    % Coordinates of contour(s) with k-th isovalue
    Q1 = Q{k}(1:3:end,:);
    Q2 = Q{k}(2:3:end,:);
     
    if chk_vt
        VT1 = VT{k}(1:3:end,:);
        VT2 = VT{k}(2:3:end,:);
    end 
    
    % Order vertices to generate either closed or open contour(s). Open 
    % contours can be identified only if the user provides the VT cell. The
    % first and last end points of open contours must be situated on mesh 
    % boundaries, and it is assumed that these are the only boundary points 
    % contained in open contours.  
    
    n = 0; 
    while ~isempty(Q1) 
        
       n = n + 1;
              
       Nq = size(Q1,1);
       if Nq==1
           if norm(Q1-Q2)>=tol % just one line segment remains by itself
               N = N + 1;
               C{N} = cat(1,Q1,Q2);
               if chk_vt
                   VT_out{N} = cat(1,VT1,VT2);
               end
           end
           break
       else
           N = N + 1;
       end
       
       X = [Q1(1,:); Q2(1,:)];
       Q1(1,:) = [];
       Q2(1,:) = [];
       Nq = Nq - 1;
       
       if chk_vt
           Y = [VT1(1,:); VT2(1,:)];
           VT1(1,:) = [];
           VT2(1,:) = [];
       end    
       
       x_end = X(end,:);
       chk_opn = false;
       for i = 1:Nq

           % Find line segment attached to the last contour vertex
           D1 = bsxfun(@minus,Q1,x_end);
           [D1,idx1] = min(sum(D1.^2,2));
           
           D2 = bsxfun(@minus,Q2,x_end);
           [D2,idx2] = min(sum(D2.^2,2));
           
           if D1<D2
                x = Q2(idx1,:);
                if chk_vt, y = VT2(idx1,:); end
                idx = idx1;
           else
               x = Q1(idx2,:);
               if chk_vt, y = VT1(idx2,:); end
               idx = idx2;
           end
           
           Q1(idx,:) = [];
           Q2(idx,:) = [];
           if chk_vt
               VT1(idx,:) = [];
               VT2(idx,:) = [];
               if y(4)==1 % x is on the boundary, break loop                   
                   X = cat(1,X,x);
                   Y = cat(1,Y,y);
                   chk_opn = true;
                   break 
               end
           end           
           
           if norm(X(1,:)-x)<=1E-12 % contour closed in on itself, can break loop
                break           
           elseif norm(x-x_end)>=tol
               X = cat(1,X,x);
               if chk_vt, Y = cat(1,Y,y); end               
               x_end = x;
           end
           
       end % end for
       
       if chk_opn
           C{N,1} = X;
           VT_out{N,1} = Y;
       else
           C{N,1} = cat(1,X,X(1,:));
           if chk_vt, VT_out{N,1} = cat(1,Y,Y(1,:)); end
       end
       
    end % end while
    
end
