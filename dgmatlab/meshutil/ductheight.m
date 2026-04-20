function h=ductheight(x)

h=0*x;
ix=x>.75 & x<=1.25;
h(ix)=.125*(cos(pi+(x(ix)-.75)/.5*pi)+1);
ix=x>=1.25 & x<1.75;
h(ix)=.25;
ix=x>=1.75 & x<2.25;
h(ix)=.125*(cos((x(ix)-1.75)/.5*pi)+1);
