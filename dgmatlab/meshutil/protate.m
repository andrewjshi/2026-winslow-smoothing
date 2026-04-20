function p=protate(p,phi)

A=[cos(phi),-sin(phi);sin(phi),cos(phi)];
p=p*A;
