function avp = GetAxesViewProps(ha)
% Get axes properties used by the 'MatchAxesView' function. 
%
% AUTHOR: Anton Semechko (a.semechko@gmail.com)
%


if nargin<1 || isempty(ha) || ~ishandle(ha) || ~strcmpi(get(ha,'type'),'axes')
    avp = struct('CameraPosition',[],'CameraTarget',[],'CameraUpVector',[],...
               'CameraViewAngle',[],'XLim',[],'YLim',[],'ZLim',[]);
    return   
end
    
avp.CameraPosition = get(ha,'CameraPosition');
avp.CameraTarget = get(ha,'CameraTarget'); 
avp.CameraUpVector = get(ha,'CameraUpVector'); 
avp.CameraViewAngle = get(ha,'CameraViewAngle');

avp.XLim = get(ha,'XLim');
avp.YLim = get(ha,'YLim');
avp.ZLim = get(ha,'ZLim');
