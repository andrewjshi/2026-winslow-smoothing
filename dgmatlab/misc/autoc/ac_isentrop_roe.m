function out = ac_isentrop_roe(UR,UL,N,Fr,GAMMA,k_pr,VG)
% ROE Fluxes
% out = ac_roe(UR,UL,N,F,GAMMA,k_pr,VG)
%
% Computes inviscid fluxes in the manner of Roe.
% UR = right state
% UL = left state
% N = outward normal
% F = resulting fluxes
% GAMMA = 1.4
% k_pr = 
% VG = grid velocity

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end
if (isnumeric(k_pr))
  k_pr = num2str(k_pr);
end

out = {
'ndotxdot=[nx]*[xdot]+[ny]*[ydot]+[nz]*[zdot]';
'rr1  = 1.0/[rr]';
'ur   = [rur]*rr1';
'vr   = [rvr]*rr1';
'wr   = [rwr]*rr1';
'u2r  = ur*ur+vr*vr+wr*wr';
'pr   = [k_pr]*[rr]^([GAMMA])';
'Er   = pr/([GAMMA]-1)+[rr]/2*u2r';
'rEr  = [rr]*Er';
'hr   = Er+pr*rr1';
'unr  = ur*[nx]+vr*[ny]+wr*[nz]';
'';
'runr = [rur]*[nx]+[rvr]*[ny]+[rwr]*[nz]';
'rhr  = rEr+pr';
'';
'rl1  = 1.0/[rl]';
'ul   = [rul]*rl1';
'vl   = [rvl]*rl1';
'wl   = [rwl]*rl1';
'u2l  = ul*ul+vl*vl+wl*wl';
'pl   = [k_pr]*[rl]^([GAMMA])';
'El   = pl/([GAMMA]-1)+[rl]/2*u2l';
'rEl  = [rl]*El';
'hl   = El+pl*rl1';
'unl  = ul*[nx]+vl*[ny]+wl*[nz]';
'';
'runl = [rul]*[nx]+[rvl]*[ny]+[rwl]*[nz]';
'rhl  = rEl+pl';
'';
'di    = sqrt([rr]*rl1)';
'd1    = 1.0/(di+1.0)';
'ui    = (di*ur+ul)*d1';
'vi    = (di*vr+vl)*d1';
'wi    = (di*wr+wl)*d1';
'hi    = (di*hr+hl)*d1';
'ci2   = ([GAMMA]-1.0)*(hi-0.5*(ui*ui+vi*vi+wi*wi))';
'ci    = sqrt(ci2)';
'af    = 0.5*(ui*ui+vi*vi+wi*wi)';
'uni   = ui*[nx]+vi*[ny]+wi*[nz]';
'';
'dr    = [rr]-[rl]';
'dru   = [rur]-[rul]';
'drv   = [rvr]-[rvl]';
'drw   = [rwr]-[rwl]';
'drE   = rEr-rEl';
'';
'rlam1 = abs(uni+ci-ndotxdot)';
'rlam2 = abs(uni-ci-ndotxdot)';
'rlam3 = abs(uni-ndotxdot)';
'';
's1    = 0.5*(rlam1+rlam2)';
's2    = 0.5*(rlam1-rlam2)';
'al1x  = ([GAMMA]-1.0)*(af*dr-ui*dru-vi*drv-wi*drw+drE)';
'al2x  = -uni*dr+dru*[nx]+drv*[ny]+drw*[nz]';
'cc1   = ((s1-rlam3)*al1x/ci2)+(s2*al2x/ci)';
'cc2   = (s2*al1x/ci)+(s1-rlam3)*al2x';
'';
'unr2=unr-ndotxdot';
'unl2=unl-ndotxdot';
'[F1]    = 0.5*([rr]*unr2+[rl]*unl2)                    - 0.5*(rlam3*dr+cc1)';
'[F2]    = 0.5*([rur]*unr2+[rul]*unl2)+0.5*[nx]*(pr+pl) - 0.5*(rlam3*dru+cc1*ui+cc2*[nx])';
'[F3]    = 0.5*([rvr]*unr2+[rvl]*unl2)+0.5*[ny]*(pr+pl) - 0.5*(rlam3*drv+cc1*vi+cc2*[ny])';
'[F4]    = 0.5*([rwr]*unr2+[rwl]*unl2)+0.5*[nz]*(pr+pl) - 0.5*(rlam3*drw+cc1*wi+cc2*[nz])';
};

replacements = {
% UR
'[rr]',  [UR '1'];
'[rur]', [UR '2'];
'[rvr]', [UR '3'];
'[rwr]', [UR '4'];
% UL
'[rl]',  [UL '1'];
'[rul]', [UL '2'];
'[rvl]', [UL '3'];
'[rwl]', [UL '4'];
% N
'[nx]',  [N  '1'];
'[ny]',  [N  '2'];
'[nz]',  [N  '3'];
% Fr
'[F1]',  [Fr '1'];
'[F2]',  [Fr '2'];
'[F3]',  [Fr '3'];
'[F4]',  [Fr '4'];
% GAMMA
'[GAMMA]', GAMMA;
% k_pr
'[k_pr]', k_pr;
% VG
'[xdot]', [VG '1'];
'[ydot]', [VG '2'];
'[zdot]', [VG '3'];
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end
