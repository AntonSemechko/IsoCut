function MatchAxesView(ha_ref,ha_src,zoom_off)

if nargin<3 || isempty(zoom_off), zoom_off=false; end

if nargin<2 || isempty(ha_src) || ~ishandle(ha_src) || ~strcmpi(get(ha_src,'type'),'axes')
    return
end

prop_names=fieldnames(GetAxesViewProps);
ha_ref=ha_ref(1);
if ishandle(ha_ref) && strcmpi(get(ha_ref,'type'),'axes')
    avp=GetAxesViewProps(ha_ref);
elseif isstruct(ha_ref) 
    if isequal(fieldnames(ha_ref),prop_names)
        avp=ha_ref;
    else
        return
    end
else
    return
end

for i=1:numel(prop_names)
    if strcmpi(prop_names{i},'CameraViewAngle') && zoom_off
        continue
    end
    set(ha_src,prop_names{i},getfield(avp,prop_names{i})) %#ok<*GFLD>
end
drawnow
