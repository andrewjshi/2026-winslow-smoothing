function d=dnacabox2(p,a,hs,theta,pv,varargin)

d1=dnaca(protate(p,-theta),a,hs,theta,pv);
d2=dcircle(p,3,0,5);
d=ddiff(d2,d1);
