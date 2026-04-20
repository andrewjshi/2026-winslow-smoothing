function fancycamera_dolly(varargin)
%FANCYCAMERA Interactively orbit/zoom/pan with the mouse buttons.

hfig=gcbf;
if isempty(hfig)
  hfig=gcf;
end
if nargin==0
  down=get(hfig, 'windowbuttondownfcn');
  if ~strcmp(down,'fancycamera_dolly(''down''  )')
    cameramenu
  end
  set(hfig, 'windowbuttondownfcn',   'fancycamera_dolly(''down''  )')
  set(hfig, 'windowbuttonmotionfcn', 'fancycamera_dolly(''motion'')')
  set(hfig, 'windowbuttonupfcn',     'fancycamera_dolly(''up''    )')
else
  h = findobj(get(hfig,'children'), 'type', 'uimenu', 'tag', 'cm598');
  if isempty(h)
    fancycamera;
    return;
  end
  Udata = get(h(1), 'userdata');
  switch get(hfig,'selectiontype')
   case 'normal'
    Udata.mode='Orbit';
   case 'extend'
    Udata.mode='Dolly In/Out';
   case 'alt'
    Udata.mode='Dolly Horiz/Vert';
  end
  set(h(1), 'userdata', Udata);
  cameramenu(varargin{:});
end
