function [p,t] = simprefine(p,t,crit,varargin)

if nargin<3,
  warning('No refinement function -- nothing done.');
end

dim = size(p,2);
switch dim
 case 2
  [p,t] = refine2D(p,t,crit,varargin{:});
 case 3
  [p,t] = refine3D(p,t,crit,varargin{:});
 otherwise
  error('Dimension not implemented');
end

function [p,t] = refine2D(p,t,crit,varargin)

se = [1,2,3; 2,3,1];
refe3 = [1,1,1,0,0,0; 1,0,0,1,1,0; 0,1,0,0,1,1; 0,0,1,1,0,1];
refe2 = [1,0,0,0,0,1; 0,1,0,1,0,0; 0,0,1,0,1,0];

sid = [t(:,se(:,1)); t(:,se(:,2)); t(:,se(:,3))];

sid = [min(sid,[],2), max(sid,[],2)];
sid = unique(sid,'rows');

psid = 0.5*(p(sid(:,1),:)+p(sid(:,2),:));

fprintf('Total number of edges     %10d\n', size(sid,1));

if isnumeric(crit) | islogical(crit)
  if size(crit,2)==2 % edges
    [foo,ix]=ismember(sort(crit,2),sid,'rows');
    ref=logical(zeros(size(sid(:,1))));
    ref(ix)=true;
  elseif size(crit,2)==1 & size(crit,1)==size(p,1) % nodes
    refpoint=false(size(p,1),1);
    refpoint(crit)=true;
    ref = max( refpoint(sid(:,1)), refpoint(sid(:,2)));
  elseif size(crit,2)==1 & size(crit,1)==size(t,1) % triangle
    error('not implemented');
%    refpoint=false(size(p,1),1);
%    refpoint(crit)=true;
%    ref = max( refpoint(sid(:,1)), refpoint(sid(:,2)));
  end
else
  % Old way
%  refpoint = feval(crit,p,varargin{:});
%  ref = max( refpoint(sid(:,1)), refpoint(sid(:,2)));
  % New way - function gives sizes
  hsid=feval(crit,psid,varargin{:});
  L=sqrt(sum((p(sid(:,2),:)-p(sid(:,1),:)).^2,2));
  ref=hsid<L;
end

fprintf('Number of edges flagged   %10d\n', sum(ref));

a = sparse( double([sid(:,1); sid(:,2)]), double([sid(:,2); sid(:,1)]), [1:size(sid,1),1:size(sid,1)]);

more = true;
while more,
   more = false;
   newe = 0;
   for i=1:size(t,1)
       side = reshape(t(i,se(:,:)),2,3); 
       for j=1:3
           refe(j) = ref(a(side(1,j),side(2,j)));
       end
       sume = sum(refe);
       switch sume
              case 3,
                  newe = newe + 3;
              case 2,
                  more = true;
                  for j=1:3
                      ref(a(side(1,j),side(2,j))) = 1;
                  end
              case 1,
                  newe = newe + 1;
       end
   end
end

fprintf('Total edges refined       %10d\n', sum(ref));

refp = zeros(size(ref));
np = size(p,1);
p = [p; psid(ref,:)];
for i=1:size(ref,1)
    if ref(i),
       np = np+1;
       refp(i) = np;
    end
end

ne = size(t,1);
fprintf('Intial number of elements %10d\n', ne);
fprintf('Number of new elements    %10d\n', newe);

t = [t;zeros(newe,3)];
next = ne+1;

for i=1:ne
    side = reshape(t(i,se(:,:)),2,3);
    for j=1:3
        refe(j) = ref(a(side(1,j),side(2,j)));
    end
    sume = sum(refe);
    switch sume
           case 0,
           case 1, 
               n = t(i,:);
               if refe(1),
                  new = refp(a(side(1,1),side(2,1)));
                  t(i,2) = new;
                  t(next,:) = [new,n(2),n(3)]; next = next+1;
               elseif refe(2),
                  new = refp(a(side(1,2),side(2,2)));
                  t(i,3) = new;
                  t(next,:) = [n(1),new,n(3)]; next = next+1;
               else
                  new = refp(a(side(1,3),side(2,3)));
                  t(i,3) = new;
                  t(next,:) = [new,n(2),n(3)]; next = next+1;
               end
           case 3,
               n = t(i,:);
               n3(1) = refp(a(side(1,1),side(2,1)));
               n3(2) = refp(a(side(1,2),side(2,2)));
               n3(3) = refp(a(side(1,3),side(2,3)));
               t(i,:) = [n(1),n3(1),n3(3)];
               t(next:next+2,:) = [n3(1),n(2),n3(2); n3(2),n(3),n3(3); n3(1),n3(2),n3(3)];
               next = next+3;
           otherwise,
               error('Something went really wrong!'); 
    end
end

fprintf('Total number of elements  %10d\n', size(t,1));




function [p,t] = refine3D(p,t,crit,varargin)

se = [1,2,3,1,2,3; 2,3,1,4,4,4];
refe3 = [1,1,1,0,0,0; 1,0,0,1,1,0; 0,1,0,0,1,1; 0,0,1,1,0,1];
refe2 = [1,0,0,0,0,1; 0,1,0,1,0,0; 0,0,1,0,1,0];

sid = [t(:,se(:,1)); t(:,se(:,2)); t(:,se(:,3)); t(:,se(:,4)); t(:,se(:,5)); t(:,se(:,6))];

sid = [min(sid,[],2), max(sid,[],2)];
sid = unique(sid,'rows');

psid = 0.5*(p(sid(:,1),:)+p(sid(:,2),:));

fprintf('Total number of edges     %10d\n', size(sid,1));

if isnumeric(crit)
  refpoint=false(size(t,1),1);
  refpoint(crit)=true;
  refpoint=repmat(refpoint,3,1);
else
  refpoint = feval(crit,p,varargin{:});
end
ref = max( refpoint(sid(:,1)), refpoint(sid(:,2)));

fprintf('Number of edges flagged   %10d\n', sum(ref));

a = sparse( double([sid(:,1); sid(:,2)]), double([sid(:,2); sid(:,1)]), [1:size(sid,1),1:size(sid,1)]);

more = true;
while more,
   more = false;
   newe = 0;
   for i=1:size(t,1)
       side = reshape(t(i,se(:,:)),2,6); 
       for j=1:6
           refe(j) = ref(a(side(1,j),side(2,j)));
       end
       sume = sum(refe);
       switch sume
              case 6,
                  newe = newe + 7;
              case 5,
                  more = true;
                  for j=1:6
                      ref(a(side(1,j),side(2,j))) = 1;
                  end
              case 4,
                  more = true;
                  for j=1:6
                      ref(a(side(1,j),side(2,j))) = 1;
                  end
              case 3,
                  if sum(refe == refe3(1,:)) ~= 6 & sum(refe == refe3(2,:)) ~= 6 &  ...
                     sum(refe == refe3(3,:)) ~= 6 & sum(refe == refe3(4,:)) ~= 6, 
                     more = true;
                     for j=1:6
                         ref(a(side(1,j),side(2,j))) = 1;
                     end
                  else
                     newe = newe + 3;
                  end
              case 2,
                  if sum(refe == refe2(1,:)) == 6 | sum(refe == refe2(2,:)) == 6 | sum(refe == refe2(3,:)) == 6,
                     more = true;
                     for j=1:6
                         ref(a(side(1,j),side(2,j))) = 1;
                     end
                  else 
                     for k=1:4
                         if sum(refe == refe3(k,:)) == 5
                            more = true;
                            for j=1:6
                                if refe3(k,j) == 1, ref(a(side(1,j),side(2,j))) = 1; end
                            end
                         end
                     end
                  end
              case 1,
                  newe = newe + 1;
       end
   end
end

fprintf('Total edges refined       %10d\n', sum(ref));

refp = zeros(size(ref));
np = size(p,1);
p = [p; psid(ref,:)];
for i=1:size(ref,1)
    if ref(i),
       np = np+1;
       refp(i) = np;
    end
end

ne = size(t,1);
fprintf('Intial number of elements %10d\n', ne);
fprintf('Number of new elements    %10d\n', newe);

t = [t;zeros(newe,4)];
next = ne+1;

for i=1:ne
    side = reshape(t(i,se(:,:)),2,6);
    for j=1:6
        refe(j) = ref(a(side(1,j),side(2,j)));
    end
    sume = sum(refe);
    switch sume
           case 0,
           case 1, 
               n = t(i,:);
               if refe(1),
                  new = refp(a(side(1,1),side(2,1)));
                  t(i,2) = new;
                  t(next,:) = [new,n(2),n(3),n(4)]; next = next+1;
               elseif refe(2),
                  new = refp(a(side(1,2),side(2,2)));
                  t(i,3) = new;
                  t(next,:) = [n(1),new,n(3),n(4)]; next = next+1;
               elseif refe(3),
                  new = refp(a(side(1,3),side(2,3)));
                  t(i,3) = new;
                  t(next,:) = [new,n(2),n(3),n(4)]; next = next+1;
               elseif refe(4),
                  new = refp(a(side(1,4),side(2,4)));
                  t(i,4) = new;
                  t(next,:) = [new,n(2),n(3),n(4)]; next = next+1;
               elseif refe(5),
                  new = refp(a(side(1,5),side(2,5)));
                  t(i,4) = new;
                  t(next,:) = [n(1),new,n(3),n(4)]; next = next+1;
               else
                  new = refp(a(side(1,6),side(2,6)));
                  t(i,4) = new;
                  t(next,:) = [n(1),n(2),new,n(4)]; next = next+1;
               end
           case 3,
               n = t(i,:);
               if sum(refe == refe3(1,:)) == 6,
                   n3(1) = refp(a(side(1,1),side(2,1)));
                   n3(2) = refp(a(side(1,2),side(2,2)));
                   n3(3) = refp(a(side(1,3),side(2,3)));
                   t(i,:) = [n(1),n3(1),n3(3),n(4)];
                   t(next:next+2,:) = [n3(1),n(2),n3(2),n(4); n3(2),n(3),n3(3),n(4); n3(1),n3(2),n3(3),n(4)];
                   next = next+3;
               elseif sum(refe == refe3(2,:)) == 6,
                   n3(1) = refp(a(side(1,1),side(2,1)));
                   n3(2) = refp(a(side(1,4),side(2,4)));
                   n3(3) = refp(a(side(1,5),side(2,5)));
                   t(i,:) = [n(1),n3(1),n(3),n3(2)];
                   t(next:next+2,:) = [n3(1),n(2),n(3),n3(3); n3(2),n3(3),n(3),n(4); n3(2),n3(1),n(3),n3(3)];
                   next = next+3;
               elseif sum(refe == refe3(3,:)) == 6,
                   n3(1) = refp(a(side(1,2),side(2,2)));
                   n3(2) = refp(a(side(1,5),side(2,5)));
                   n3(3) = refp(a(side(1,6),side(2,6)));
                   t(i,:) = [n(1),n(2),n3(1),n3(2)];
                   t(next:next+2,:) = [n(1),n3(1),n(3),n3(3); n(1),n3(2),n3(3),n(4); n(1),n3(2),n3(1),n3(3)];
                   next = next+3;
               else
                   n3(1) = refp(a(side(1,3),side(2,3)));
                   n3(2) = refp(a(side(1,4),side(2,4)));
                   n3(3) = refp(a(side(1,6),side(2,6)));
                   t(i,:) = [n(1),n(2),n3(1),n3(2)];
                   t(next:next+2,:) = [n3(1),n(2),n(3),n3(3); n3(2),n(2),n3(3),n(4); n3(2),n(2),n3(1),n3(3)];
                   next = next+3;
               end
           case 6,
               n = t(i,:);
               for j=1:6
                   nn(j) = refp(a(side(1,j),side(2,j)));
               end
               t(i,:) = [n(1),nn(1),nn(3),nn(4)];
               t(next:next+6,:) = [nn(1),n(2),nn(2),nn(5); nn(2),n(3),nn(3),nn(6); nn(4),nn(5),nn(6),n(4); ...
                                   nn(1),nn(2),nn(3),nn(6); nn(1),nn(2),nn(6),nn(5); nn(1),nn(3),nn(4),nn(6); ... 
                                   nn(1),nn(6),nn(4),nn(5)];
               next = next+7;
           otherwise,
               error('Something went really wrong!'); 
    end
end

fprintf('Total number of elements  %10d\n', size(t,1));


function ref = dummy(p)
ref = logical(zeros(size(p(:,1))));

