porder = 4;
msh = mshcircle(1);
msh = nodealloc(msh,porder);
curvep1 = msh.p1;

msh = mshcurved(msh,[]);
msh = nodealloc(msh,porder);
linearp1 = msh.p1;

doplot = true;
msh1 = elliptic_smoothing(msh, curvep1, doplot);
