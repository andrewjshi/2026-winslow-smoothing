function d=dnacabox(p,a,hs,theta,varargin)

d1=dnaca(protate(p,-theta),a,hs,theta,varargin{:});
d2=dcircle(p,-1,0,5);
d3=dellipse(pshift(p,-4,0),[5,2.5]);
d=ddiff(dintersect(d3,d2),d1);
