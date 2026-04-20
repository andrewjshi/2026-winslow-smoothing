function qual=quadqual(p,q,type,signed)
%QUADQUAL Quadrilateral quality.
%   QUAL=QUADQUAL(P,Q,TYPE)
%
%   TYPE == 1: Lo (default)

if nargin<3
  type=1;
end
if nargin<4
  signed=0;
end

switch type
 case 1
  alpha1=simpqual(p,q(:,[1,2,4]),2,1);
  alpha2=simpqual(p,q(:,[1,4,3]),2,1);
  alpha3=simpqual(p,q(:,[1,2,3]),2,1);
  alpha4=simpqual(p,q(:,[2,4,3]),2,1);
  alphas=[alpha1,alpha2,alpha3,alpha4];
  inverted=any(alphas<0,2);
  [foo,ix]=sort(abs(alphas),2);
  nq=size(q,1);
  alphas=alphas(repmat((1:nq)',[1,4])+nq*(ix-1));
  qual=alphas(:,1).*alphas(:,2)./alphas(:,3)./alphas(:,4);
  qual(inverted) = 0;
  case 2
    vab = p(q(:,2),:) - p(q(:,1),:);
    vbc = p(q(:,4),:) - p(q(:,2),:);
    vcd = p(q(:,3),:) - p(q(:,4),:);
    vda = p(q(:,1),:) - p(q(:,3),:);
    ab = sqrt(sum(vab.^2,2));
    bc = sqrt(sum(vbc.^2,2));
    cd = sqrt(sum(vcd.^2,2));
    da = sqrt(sum(vda.^2,2));
    sin1 = (vda(:,1).*vab(:,2) - vda(:,2).*vab(:,1)) ./ da ./ ab;
    sin2 = (vab(:,1).*vbc(:,2) - vab(:,2).*vbc(:,1)) ./ ab ./ bc;
    sin3 = (vbc(:,1).*vcd(:,2) - vbc(:,2).*vcd(:,1)) ./ bc ./ cd;
    sin4 = (vcd(:,1).*vda(:,2) - vcd(:,2).*vda(:,1)) ./ cd ./ da;
    ka = (da.^2 + ab.^2) ./ da ./ ab ./ sin1;
    kb = (ab.^2 + bc.^2) ./ ab ./ bc ./ sin2;
    kc = (bc.^2 + cd.^2) ./ bc ./ cd ./ sin3;
    kd = (cd.^2 + da.^2) ./ cd ./ da ./ sin4;
    qual = 4 ./ sqrt(sum([ka,kb,kc,kd].^2,2));
    inverted = any([ka,kb,kc,kd]<0,2);
    qual(inverted) = 0;
 otherwise
  error('Incorrect type.');
end

function alpha = sin_of_angle_between_edges(v1,v2)

v1abs = sqrt(sum(v1.^2,2));
v2abs = sqrt(sum(v2.^2,2));
v1v2 = cross(v1,v2,2);
alpha = v1v2./v1abs./v2abs;
