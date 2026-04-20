function dgimgen(fn,frame,sz,res,format)

if nargin<3, sz=[7,7]; end
if nargin<4, res=100; end
if nargin<5, format='jpg'; end

colormap(jet(65536));
delete(findobj(gcf,'tag','Colorbar'));
axis off
fn=sprintf('%s%05d.%s',fn,frame,format);
set(gcf,'paperpos',[0,0,sz]);
%axis equal,axis tight

switch format
 case 'png'
  print('-dpng',sprintf('-r%d',res),fn);
 case 'jpg'
  print('-djpeg90',sprintf('-r%d',res),fn);
 otherwise
  error('Unknown format');
end

