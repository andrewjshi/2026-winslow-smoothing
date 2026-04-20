%%%%%%%%%%%%%
%% Fluxes

in={'U',{'u','v'}; ...
    'Q',{'ux','vx', ...
         'uy','vy'}};
pars={'pars',{'lambda','mu','rho'}};
def={'e11=ux';
     'e12=(vx+uy)/2';
     'e21=(uy+vx)/2';
     'e22=vy';
     '';
     's11=lambda*(e11+e22)+2*mu*e11';
     's12=2*mu*e12';
     's21=2*mu*e21';
     's22=lambda*(e11+e22)+2*mu*e22'};

out={'FF',{'-s11/rho','-s21/rho', ...
           '-s12/rho','-s22/rho'}};

clear opts
opts.dojac=true;
opts.jacperm={{[2,2,2],[1,3,2]},{[2,2,2,2],[1,3,2,4]}};
[flux2d,flux2d_deriv,min,mout]=autoc_new(in,pars,def,out,'flux2d',opts);

%%%%%%%%%%%%%%%%%%
%% Output file

createc({flux2d,flux2d_deriv},'linelast_autogen.cpp');
