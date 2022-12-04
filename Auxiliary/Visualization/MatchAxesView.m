function hLink = MatchAxesView(ha_ref,ha_src,zoom_off,consenLim)
% Match view of source axes to a reference axes.
%
% INPUT
%   - ha_ref    : reference axes handle OR an axes-view-props structure
%                 containing desired camara properties 
%   - ha_src    : handle of the source axes
%   - zoom_off  : (optional) set zoom_off = true to disable synching of the 
%                 CameraViewAngle property. zoom_off = false is the defaut
%                 setting.
%   - consenLim : (optional) specify consenLim = true to set XLim, YLim,
%                 and ZLim properties so they contain limits of both the 
%                 source and reference axes. For example,
%                 XLim = [min(ha_ref.XLim(1),ha_src.XLim(1)) max(ha_ref.XLim(2),ha_src.XLim(2))]
%                 consenLim = false is the default setting. Note this
%                 option is used only if ha_ref is a handle.
% 
% OUTPUT
%   - hLink     : if ha_ref is a handle, hLink is an object that links
%                 camera properties of the source and reference axes. Note
%                 this object MUST remain in memory in order for the
%                 axes properties to remain linked. Deleting hLink will
%                 break the link.
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<3 || isempty(zoom_off), zoom_off = false; end
if nargin<4 || isempty(consenLim), consenLim = false; end

if nargin<2 || isempty(ha_src) || ~ishandle(ha_src) || ~strcmpi(get(ha_src,'type'),'axes')
    return
end

prop_names = fieldnames(GetAxesViewProps);
ha_ref = ha_ref(1);

chk_ref_handle = false;
hLink = [];
if ishandle(ha_ref) && strcmpi(get(ha_ref,'type'),'axes')
    chk_ref_handle = true;
    avp = GetAxesViewProps(ha_ref);
elseif isstruct(ha_ref) 
    if isequal(fieldnames(ha_ref),prop_names)
        avp = ha_ref;
    else
        return
    end
else
    return
end


chk_link = false(size(prop_names));
for i = 1:numel(prop_names)
    
    if strcmpi(prop_names{i},'CameraViewAngle') && zoom_off
        continue
    end
    
    chk_link(i) = true;
    
    if chk_ref_handle && consenLim && contains(prop_names{i},'Lim')
        limVal_ref = get(ha_ref,prop_names{i});
        limVal_src = get(ha_src,prop_names{i});
        limVal = [min(limVal_ref(1),limVal_src(1)) max(limVal_ref(2),limVal_src(2))];
        set([ha_ref ha_src],prop_names{i},limVal)
    else
        set(ha_src,prop_names{i},getfield(avp,prop_names{i})) %#ok<*GFLD>
    end
    
end
drawnow

if chk_ref_handle
    hLink = linkprop([ha_ref ha_src],prop_names(chk_link));
end

%if nargout<1, clear hLink; end

