function out = ac_inverse(A, detA, B)
% Inverse
% out = ac_inverse(A, detA, B)
%
% A and B are 3-by-3

out = {'BBB11 = (AAA22*AAA33 - AAA23*AAA32) / detAAA';
       'BBB12 = (AAA13*AAA32 - AAA12*AAA33) / detAAA';
       'BBB13 = (AAA12*AAA23 - AAA13*AAA22) / detAAA';
       'BBB21 = (AAA23*AAA31 - AAA21*AAA33) / detAAA';
       'BBB22 = (AAA11*AAA33 - AAA13*AAA31) / detAAA';
       'BBB23 = (AAA13*AAA21 - AAA11*AAA23) / detAAA';
       'BBB31 = (AAA21*AAA32 - AAA22*AAA31) / detAAA';
       'BBB32 = (AAA12*AAA31 - AAA11*AAA32) / detAAA';
       'BBB33 = (AAA11*AAA22 - AAA12*AAA21) / detAAA'};

out = strrep(out, 'detAAA', detA);
out = strrep(out, 'AAA', A);
out = strrep(out, 'BBB', B);
