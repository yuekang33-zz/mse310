%  hsdLPsolversub.m
%
% find the optimal solution
%
cvx = ((gamma*mu)*ee)./x - s;
cvz = ((gamma*mu)*eez)./z - w;
woz = w./z;
r1  = cvx -rd;
r1(bindx) = r1(bindx) - cvz - rb.*woz;
r2  = c;
r22 = c;
r2(bindx)  = r2(bindx)  - u.*woz;
r22(bindx) = r22(bindx) + u.*woz;
d   = s./x;
d(bindx) = d(bindx)+woz;
d   = sqrt(d);
AD = sparse(i,j,v./d(j), m, n);
%
% Find a scaling parameter
%
%alpha = (max([abs(max(vd));abs(min(vd))]))/1000;
%
r1  = [  r1./d  ; rp];
r2  = [ -r2./d  ;  b];
r22 = [-r22./d  ;  b];
%
ss=[speye(n) AD'; AD sparse(m,m)]\[r1 r2];
clear r1 r2 AD
ss(n+1:n+m,:)=-ss(n+1:n+m,:);
%
% get dtau
%
dtau = (gamma*mu/tau - kappa + rg + u'*cvz + u'*(woz.*rb) - r22'*ss(:,1));
dtau = dtau/(kappa/tau+u'*(woz.*u)+ r22'*ss(:,2));
clear r22
ss   = ss(:,1)+dtau*ss(:,2);
%
% get dx 
%
dx  = ss(1:n)./d;
clear d
%
dy  = ss(n+1:n+m);
clear ss
%
% get ds 
%
ds  = cvx - (s.*dx)./x;
clear cvx
%
% get dz and dw if bounds exist
%
if nb > 0.5,
  dz  = dtau*u - dx(bindx) - rb;
  dw  = cvz - woz.*dz;
  clear cvz woz
end;
%
% get dkappa
%
dkappa= gamma*mu/tau - kappa - kappa*dtau/tau;
%
if nb > 0.5;
  ratp  = beta/abs(min([dx./x;dz./z;dtau/tau;dkappa/kappa]));
  ratd  = beta/abs(min([ds./s;dw./w;dtau/tau;dkappa/kappa]));
else
  ratp  = beta/abs(min([dx./x;dtau/tau;dkappa/kappa]));
  ratd  = beta/abs(min([ds./s;dtau/tau;dkappa/kappa]));
end;
x    = x    + ratp*dx;
y    = y    + ratd*dy;
s    = s    + ratd*ds;
clear dx dy ds
if nb > 0.5,
  z    = z    + ratp*dz;
  w    = w    + ratd*dw;
  clear dz dw
end;
taup = tau  + ratp*dtau;
taud = tau  + ratd*dtau;
tau  = min([taup; taud]);
if taup <= taud,
  kappa = kappa + ratp*dkappa;
  y    = (tau/taud)*y;
  s    = (tau/taud)*s;
  if nb > 0.5, w = (tau/taud)*w; end;
else
  kappa = kappa + ratd*dkappa;
  x     = (tau/taup)*x;
  if nb > 0.5, z = (tau/taud)*z; end;
end;
%return





