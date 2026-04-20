function msh = WBmesh( nlay, nh1, nh2, nv1, nv2, nv3, nl1, nl2, nl3)
%WBMESH  Generate tetrahedral Wedge-Box mesh.
%
%   Symtax: msh = WBmesh( nlay, nh1, nh2, nh3, nv1, nv2, nl1, nl2, nl3)
%
%   Default: msh = WBmesh(5,5,3,3,5,3,3,7,3);

%Number of subdivisions

if nargin < 1,
%   nlay = 9;  % Number of viscous layers + 1
%   nh1 = 5;   % Number of horizontal points y-direction in inner region
%   nh2 = 3;   % Number of horizontal points y-direction in inner region
%   nv1 = 3;   % Number of vertical points z-direction in lower region
%   nv2 = 5;   % Number of vertical points z-direction in middle region
%   nv3 = 3;   % Number of vertical points z-direction in upper region
%   nl1 = 3;   % Number of horizontal points x-direction pre-body region
%   nl2 = 7;   % Number of horizontal points x-direction body region
%   nl3 = 3;   % Number of horizontal points x-direction aft-body region
    nlay = 5;  % Number of viscous layers + 1
    nh1 = 5;   % Number of horizontal points y-direction in inner region
    nh2 = 5;   % Number of horizontal points y-direction in inner region
    nv1 = 5;   % Number of vertical points z-direction in lower region
    nv2 = 5;   % Number of vertical points z-direction in middle region
    nv3 = 5;   % Number of vertical points z-direction in upper region
    nl1 = 5;   % Number of horizontal points x-direction pre-body region
    nl2 = 8;  % Number of horizontal points x-direction body region
    nl3 = 5;   % Number of horizontal points x-direction aft-body region
end

%Dimensions

bl = 10;
l1 = 100;
l2 = 200;
h1 = 60;
h2 = 80;
w  = 100;

pwedge = [    0,    0,    0;
              9,    0, 22.5;
             63,    0, 22.5;
             63,    0,    0;
             18,   18,    0;
             27,   18, 22.5;
             63,   18, 22.5;
             63,   18,    0];
   
% Boundary layer region

abl = bl*sin(pi/8);
cbl = bl/tan(pi/8);
dbl = 9*(22.5+bl)/22.5;

pwedbl = [       -cbl,     -bl,       0;
              dbl-cbl,     -bl, 22.5+bl;
                63+bl,     -bl, 22.5+bl;
                63+bl,     -bl,       0;
               18-abl,   18+bl,       0;
           18+dbl-abl,   18+bl, 22.5+bl;
                63+bl,   18+bl, 22.5+bl;
                63+bl,   18+bl,       0];
            
%Outer Box


pc1 =    [   -l1,           -h1,      0;
             -l1,           -h1,  0.5*w;
             -l1,           -h1,      w;
             -l1,  (-2*h1+h2)/3,      0;
             -l1,  (-2*h1+h2)/3,  0.5*w;
             -l1,  (-2*h1+h2)/3,      w;
             -l1,  (-h1+2*h2)/3,      0;
             -l1,  (-h1+2*h2)/3,  0.5*w;
             -l1,  (-h1+2*h2)/3,      w;
             -l1,            h2,      0;
             -l1,            h2,  0.5*w;
             -l1,            h2,      w];
         
pc2 = pc1; pc2(:,1) = l2;

pm1 = [ (-2*l1+l2)/3,             -h1,      0;
        (-2*l1+l2)/3,             -h1,  0.5*w;
        (-2*l1+l2)/3,             -h1,      w;
        (-2*l1+l2)/3,    (-2*h1+h2)/3,      w;
        (-2*l1+l2)/3,    (-h1+2*h2)/3,      w;
        (-2*l1+l2)/3,              h2,      0;
        (-2*l1+l2)/3,              h2,  0.5*w;
        (-2*l1+l2)/3,              h2,      w];
pm2 = pm1; pm2(:,1) = (-l1+2*l2)/3;

fwedge = [  1, 4, 3, 2;
            5, 6, 7, 8;
            1, 2, 6, 5;
            2, 3, 7, 6;
            3, 4, 8, 7];
        
points = [ pwedge; pwedbl; pc1; pm1; pm2; pc2];

hex(1:22,:) = [  9, 10,  1,  2, 13, 14,  5,  6;
                 2, 10,  3, 11,  6, 14,  7, 15;
                 4,  3, 12, 11,  8,  7, 16, 15;
                 5,  6,  8,  7, 13, 14, 16, 15;
                 9, 10, 12, 11,  1,  2,  4,  3;
                20, 21,  9, 10, 23, 24, 13, 14;
                17, 18, 29, 30, 20, 21,  9, 10;
                18, 19, 30, 31, 21, 22, 10, 32;
                21, 22, 10, 32, 24, 25, 14, 33;
                24, 25, 14, 33, 27, 28, 35, 36;
                23, 24, 13, 14, 26, 27, 34, 35;
                29, 30, 37, 38,  9, 10, 12, 11;
                30, 31, 38, 39, 10, 32, 11, 40;
                10, 32, 11, 40, 14, 33, 15, 41;
                14, 33, 15, 41, 35, 36, 43, 44;
                13, 14, 16, 15, 34, 35, 42, 43;
                37, 38, 45, 46, 12, 11, 48, 49;
                38, 39, 46, 47, 11, 40, 49, 50;
                11, 40, 49, 50, 15, 41, 52, 53;
                15, 41, 52, 53, 43, 44, 55, 56;
                16, 15, 51, 52, 42, 43, 54, 55;
                12, 11, 48, 49, 16, 15, 51, 52];
            
X1  = block( [nh1, nlay, nv2], [-2,-5,-2,-5, 0.2,0.2,0.2,0.2, -2,-2,-5,-5], points(hex(1,:),:));
X2  = block( [nlay, nl2, nv2], [5,5,5,5, -5,-2,-5,-2, -5,-2,-5,-2], points(hex(2,:),:));
X3  = block( [nh1, nlay, nv2], [-5,-2,-5,-2, 5,5,5,5, -5,-5,-2,-2], points(hex(3,:),:));
X4  = block( [nh1, nl2, nlay], [-5,-5,-2,-2, -5,-5,-2,-2, 5,5,5,5], points(hex(4,:),:)); 
X5  = block( [nh1, nl2, nlay], [-2,-2,-5,-5, -2,-2,-5,-5, 0.2,0.2,0.2,0.2], points(hex(5,:),:));
X6  = block( [nh1, nl1, nv2], [1,-2,1,-2, 0.2,0.2,0.2,0.2, 1,1,-2,-2], points(hex(6,:),:));
X7  = block( [nh1, nl1, nv1], [1,1,1,-2, 0.2,0.2,0.2,0.2, 1,1,0.15,0.15], points(hex(7,:),:));
X8  = block( [nh2, nl1, nv1], [1,1,1,8, 0.2,0.2,0.2,0.2, 1,1,0.15,1], points(hex(8,:),:));
X9  = block( [nh2, nl1, nv2], [1,8,1,8, 0.2,0.2,0.2,0.2, 1,1,-2,1], points(hex(9,:),:));
X10 = block( [nh2, nl1, nv3], [1,8,1,1, 0.2,0.2,0.2,0.2, 1,1,6,1], points(hex(10,:),:));
X11 = block( [nh1, nl1, nv3], [1,-2,1,1, 0.2,0.2,0.2,0.2, 1,1,6,6], points(hex(11,:),:));
X12 = block( [nh1, nl2, nv1], [1,1,-2,-2, 1,1,-2,-2, 0.15,0.15,0.15,0.15], points(hex(12,:),:));
X13 = block( [nh2, nl2, nv1], [1,1,8,8, 1,1,-2,1, 0.15,1,0.15,1], points(hex(13,:),:));
X14 = block( [nh2, nl2, nv2], [8,8,8,8, -2,1,-2,1, -2,1,-2,1], points(hex(14,:),:));
X15 = block( [nh2, nl2, nv3], [8,8,1,1, -2,1,1,1, 6,1,6,1], points(hex(15,:),:));
X16 = block( [nh1, nl2, nv3], [-2,-2,1,1, -2,-2,1,1, 6,6,6,6], points(hex(16,:),:));
X17 = block( [nh1, nl3, nv1], [1,1,-2,1, 2,2,5,5, 0.15,0.15,1,1], points(hex(17,:),:));
X18 = block( [nh2, nl3, nv1], [1,1,8,1, 2,1.0,5,2, 0.15,1,1,1], points(hex(18,:),:));
X19 = block( [nh2, nl3, nv2], [8,1,8,1, 5,2,5,2, -2,1,1,1], points(hex(19,:),:));
X20 = block( [nh2, nl3, nv3], [8,1,1,1, 5,2,2,1, 6,1,1,1], points(hex(20,:),:));
X21 = block( [nh1, nl3, nv3], [-2,1,1,1, 5,5,2,2, 6,6,1,1], points(hex(21,:),:));
X22 = block( [nh1, nl3, nv2], [-2,1,-2,1, 5,5,5,5, -2,-2,1,1], points(hex(22,:),:));

%X1  = block( [nh1, nlay, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(1,:),:));
%X2  = block( [nlay, nl2, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(2,:),:));
%X3  = block( [nh1, nlay, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(3,:),:));
%X4  = block( [nh1, nl2, nlay], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(4,:),:)); 
%X5  = block( [nh1, nl2, nlay], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(5,:),:));
%X6  = block( [nh1, nl1, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(6,:),:));
%X7  = block( [nh1, nl1, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(7,:),:));
%X8  = block( [nh2, nl1, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(8,:),:));
%X9  = block( [nh2, nl1, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(9,:),:));
%X10 = block( [nh2, nl1, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(10,:),:));
%X11 = block( [nh1, nl1, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(11,:),:));
%X12 = block( [nh1, nl2, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(12,:),:));
%X13 = block( [nh2, nl2, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(13,:),:));
%X14 = block( [nh2, nl2, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(14,:),:));
%X15 = block( [nh2, nl2, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(15,:),:));
%X16 = block( [nh1, nl2, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(16,:),:));
%X17 = block( [nh1, nl3, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(17,:),:));
%X18 = block( [nh2, nl3, nv1], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(18,:),:));
%X19 = block( [nh2, nl3, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(19,:),:));
%X20 = block( [nh2, nl3, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(20,:),:));
%X21 = block( [nh1, nl3, nv3], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(21,:),:));
%X22 = block( [nh1, nl3, nv2], [ 1,1,1,1,1,1,1,1,1,1,1,1], points(hex(22,:),:));

% Triangulation size

n = [size(X1 ); size(X2 ); size(X3 ); size(X4 ); size(X5 ); size(X6 ); size(X7 ); size(X8 ); ...
     size(X9 ); size(X10); size(X11); size(X12); size(X13); size(X14); size(X15); size(X16); ...
     size(X17); size(X18); size(X19); size(X20); size(X21); size(X22)];
 
n = (n(:,2:4)-1)'; nhex = sum(prod(n));

% XT = zeros(3,4,5*nhex);

istart=1;
ind=istart;                       XT = reshape(block2tets(X1,ind),12,[]);
ind=ind*(-1)^(nh1-1);             XT = [XT,reshape(block2tets(X2,ind),12,[])];
ind=ind*(-1)^(nl2-1);             XT = [XT,reshape(block2tets(X3,ind),12,[])];
ind=istart*(-1)^(nlay-1+nv2-1);   XT = [XT,reshape(block2tets(X4,ind),12,[])];
ind=istart;                       XT = [XT,reshape(block2tets(X5,ind),12,[])];

i1=istart*(-1)^(nl1-1);           XT = [XT,reshape(block2tets(X7,i1),12,[])];
i2=i1*(-1)^(nh1-1);               XT = [XT,reshape(block2tets(X8,i2),12,[])];
i3=i1*(-1)^(nv1-1);               XT = [XT,reshape(block2tets(X6,i3),12,[])];
i4=i3*(-1)^(nh1-1);               XT = [XT,reshape(block2tets(X9,i4),12,[])];
i5=i3*(-1)^(nv2-1);               XT = [XT,reshape(block2tets(X11,i5),12,[])];
i6=i5*(-1)^(nh1-1);               XT = [XT,reshape(block2tets(X10,i6),12,[])];

i1=i1*(-1)^(nl1-1);               XT = [XT,reshape(block2tets(X12,i1),12,[])];
i2=i2*(-1)^(nl1-1);               XT = [XT,reshape(block2tets(X13,i2),12,[])];
i3=i3*(-1)^(nl1-1);
i4=i4*(-1)^(nl1-1);               XT = [XT,reshape(block2tets(X14,i4),12,[])];
i5=i5*(-1)^(nl1-1);               XT = [XT,reshape(block2tets(X16,i5),12,[])];
i6=i6*(-1)^(nl1-1);               XT = [XT,reshape(block2tets(X15,i6),12,[])];

i1=i1*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X17,i1),12,[])];
i2=i2*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X18,i2),12,[])];
i3=i3*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X22,i3),12,[])];
i4=i4*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X19,i4),12,[])];
i5=i5*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X21,i5),12,[])];
i6=i6*(-1)^(nl2-1);               XT = [XT,reshape(block2tets(X20,i6),12,[])];

XT = reshape(XT,3,4,[]);

p0=reshape(XT,3,[])';
t0=reshape(1:prod(size(XT))/3,4,size(XT,3))';
[p,t]=fixmesh(p0,t0,1.0e-6);
bndexpr={'all(p(:,3)<1.0e-5)', ...
         'all((p(:,1) > -1) & (p(:,1) < 64) & (p(:,2) > -1) & (p(:,2) < 19) & (p(:,3) < 23))', ...
         'any((p(:,1) < -1) | (p(:,1) > 64) | (p(:,2) < -1) | (p(:,2) > 19) | (p(:,3) > 23))'};
msh=ml2msh(p,t,bndexpr);
msh=mshcurved(msh,[]);

% Check boundary faces

alay = nlay-1;
ah1 = nh1-1;
ah2 = nh2-1;
av1 = nv1-1;
av2 = nv2-1;
av3 = nv3-1;
al1 = nl1-1;
al2 = nl2-1;
al3 = nl3-1;

symmetry = 2*((al1+al2+al3)*(av1+av2+av3) - al2*av2) + 4*alay*(av2+al2);
wall = 4*ah1*(av2+al2) + 2*av2*al2;
farfield = 4*(ah1+ah2)*(av1+av2+av3) + 4*(ah1+ah2)*(al1+al2+al3) + 2*(al1+al2+al3)*(av1+av2+av3);

if (sum(msh.t2t(:) == -1) ~= symmetry),  error('Parity error 1');  end
if (sum(msh.t2t(:) == -2) ~= wall),  error('Parity error 2');  end
if (sum(msh.t2t(:) == -3) ~= farfield),  error('Parity error 3');  end
