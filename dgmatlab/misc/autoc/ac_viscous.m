function out = ac_viscous(U, Ux, Fv, GAMMA, RE, Pr)
% VISCOUS Fluxes
% out = ac_viscous(U, Ux, Fv, GAMMA, RE, Pr)
%
% WATCH OUT FOR SIGNS
%  -- this might return the negative of what you would expect.

if (isnumeric(GAMMA))
  GAMMA = num2str(GAMMA);
end
if (isnumeric(RE))
  RE = num2str(RE);
end
if (isnumeric(Pr))
  Pr = num2str(Pr);
end

out = { ...
'u=[ru]/[r]';
'v=[rv]/[r]';
'w=[rw]/[r]';
'E=[rE]/[r]';
'ux=([rux]-[rx]*u)/[r]';
'uy=([ruy]-[ry]*u)/[r]';
'uz=([ruz]-[rz]*u)/[r]';
'vx=([rvx]-[rx]*v)/[r]';
'vy=([rvy]-[ry]*v)/[r]';
'vz=([rvz]-[rz]*v)/[r]';
'wx=([rwx]-[rx]*w)/[r]';
'wy=([rwy]-[ry]*w)/[r]';
'wz=([rwz]-[rz]*w)/[r]';
'Ex=([rEx]-[rx]*E)/[r]';
'Ey=([rEy]-[ry]*E)/[r]';
'Ez=([rEz]-[rz]*E)/[r]';
'';
'pxx=2/3*(2*ux-vy-wz)';
'pxy=vx+uy';
'pxz=wx+uz';
'pyz=wy+vz';
'pyy=2/3*(2*vy-ux-wz)';
'pzz=2/3*(2*wz-ux-vy)';
'ex=Ex-u*ux-v*vx-w*wx';
'ey=Ey-u*uy-v*vy-w*wy';
'ez=Ez-u*uz-v*vz-w*wz';
'[FF]11=-0';
'[FF]21=-pxx/[RE]';
'[FF]31=-pxy/[RE]';
'[FF]41=-pxz/[RE]';
'[FF]51=-(pxx*u+pxy*v+pxz*w+[GAMMA]/[Pr]*ex)/[RE]';
'[FF]12=-0';
'[FF]22=-pxy/[RE]';
'[FF]32=-pyy/[RE]';
'[FF]42=-pyz/[RE]';
'[FF]52=-(pxy*u+pyy*v+pyz*w+[GAMMA]/[Pr]*ey)/[RE]';
'[FF]13=-0';
'[FF]23=-pxz/[RE]';
'[FF]33=-pyz/[RE]';
'[FF]43=-pzz/[RE]';
'[FF]53=-(pxz*u+pyz*v+pzz*w+[GAMMA]/[Pr]*ez)/[RE]';
};

replacements = {
% U
'[r]',  [U '1'];
'[ru]', [U '2'];
'[rv]', [U '3'];
'[rw]', [U '4'];
'[rE]', [U '5'];
% Ux
% _x
'[rx]',  [Ux '11'];
'[rux]', [Ux '21'];
'[rvx]', [Ux '31'];
'[rwx]', [Ux '41'];
'[rEx]', [Ux '51'];
% _y
'[ry]',  [Ux '12'];
'[ruy]', [Ux '22'];
'[rvy]', [Ux '32'];
'[rwy]', [Ux '42'];
'[rEy]', [Ux '52'];
% _z
'[rz]',  [Ux '13'];
'[ruz]', [Ux '23'];
'[rvz]', [Ux '33'];
'[rwz]', [Ux '43'];
'[rEz]', [Ux '53'];
% Fv
'[FF]',  Fv;
% GAMMA
'[GAMMA]', GAMMA;
% RE
'[RE]', RE;
% Pr
'[Pr]', Pr;
};

for r=replacements'
  out = strrep(out,r{1},r{2});
end

