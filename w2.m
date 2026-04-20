porder = 4;
n = 10;
msh = mshsquare(n+1,n+1);
msh = nodealloc(msh, porder);

[~,edg] = getbndnodes(msh, dginit(msh), 4);
x = msh.p1(:,1,:);
y = msh.p1(:,2,:);

newx = x;
newy = y;
newx(edg) = x(edg) + 0.5*sin(pi*y(edg));
curvep1 = cat(2, newx, newy);
msh_tangled = msh;
msh_tangled.p1 = curvep1;
dgmeshplot_curved(msh_tangled, 4, 1, 0); hold on; set(gcf,'Position',[114 1 560 420])

doplot = true;
msh1 = elliptic_smoothing(msh, curvep1, doplot);
