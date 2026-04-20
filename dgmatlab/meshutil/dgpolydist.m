function d=dgpolydist(p1,poly)

d=permute(reshape(dpoly(reshape(permute(p1,[1,3,2]),[],2),poly),[size(p1,1),size(p1,3)]),[1,3,2]);
