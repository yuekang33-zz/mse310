%  LP Solver: hsdLPsolver.m
%
%  See user guide at the end of this file.
%
 if exist('toler') ~= 1
   toler = 1.e-8;
 end;
 if exist('beta') ~= 1
   beta = .995;
 end;
 if exist('bindx') ~= 1
   bindx = [];
 end;
 if exist('findx') ~= 1
   findx = [];
 end;
 [nb,x]  = size(bindx);
 [nf,x]  = size(findx);
 norbc = 1+max([c;b])-min([c;b]);
 c = c/norbc;
 b = b/norbc;
 if nb > .5
     u=u/norbc;
 end
 if nf > .5,
   A = [A  -A(:,findx)];
   c = [c ; - c(findx)];
 end;
 [m,n] = size(A);
 [i,j,v] = find(A);
 ee = ones(n,1);
%
% Initialization
%
% b    = b/(norm(b)+1);
% c    = c/(norm(c)+1);
 x    = ((norm(b)+1)/sqrt(n))*ones(n,1);
 s    = ((norm(c)+1)/sqrt(n))*ones(n,1);
 if nb > 0.5,
   z    = ((norm(b)+1)/sqrt(n))*ones(nb,1);
   w    = ((norm(c)+1)/sqrt(n))*ones(nb,1);
   eez  = ones(nb,1);
 else
   z   = 1;
   w   = 0;
   eez = 0;
   rb  = 0;
   u   = 0;
 end;
 y    = zeros(m,1);
 tau0  = 1;
 kappa0= (norm(c)+1)*(norm(b)+1)/n;
% tau0  = (norm(b)+1)/sqrt(n);
% kappa0= (norm(c)+1)/sqrt(n);
 tau   = tau0;
 kappa = kappa0;
%
% Compute residuals
%
 mu0  = (x'*s+z'*w+tau*kappa)/(n+nb+1);;
 mu   = mu0;
 rp  = tau*b - A*x;
%
 if nb > 0.5,  
   rb  = x(bindx) + z - tau*u; 
 end;
%
 rd  = tau*c - A'*y -s;
 rd(bindx) = rd(bindx) + w;
 obp = c'*x;
 obd = b'*y-u'*w;
 rg  = obp - obd + kappa;
 zh   = (obp-obd)/tau;
%
% Iteration
%
 gamma=1/sqrt(n+nb+1);
 go   = 1;
 iter=1;
 while go >= toler,
%  [x,y,s,tau,kappa]=sphslfsubb(A,b,c,x,y,s,tau,kappa,mu,rp,rd,rg,beta,gamma);
   hsdLPsolversub
%
%  Compute new residuals
%
%   tau/kappa
   mu  = (x'*s+z'*w+tau*kappa)/(n+nb+1),
   rp  = tau*b - A*x;
   if nb > 0.5,  
     rb  = x(bindx) + z - tau*u; 
   end;
   rd  = tau*c - A'*y -s;
   rd(bindx) = rd(bindx) + w;
   obp = c'*x;
   obd = b'*y-u'*w;
   rg  = obp - obd + kappa;
%   pause
   go = max([norm([rp;rb])/(tau+norm([x;z])); norm(rd)/(tau+norm([s;w]))]);
   go = max([go; abs(obp-obd)/(tau+abs(obd))]);
%
%  Check infeasibility
% 
   if (tau*kappa0/(tau0*kappa) < toler) & (mu/mu0 < toler/n),
     if (obp < 0) & (obd > 0),
      disp('Both primal and dual are infeasible.'); go=-1;
     end;
     if obp < 0, 
      disp('Dual is infeasible.')  ; go=-1; 
     end;
     if obd > 0,
      disp('Primal is infeasible.'); go=-1;
     end;
   end;
   zh   = [zh (obp-obd)/tau];
%   
%  Adjust centering weight
%
   gamma = 1/sqrt(n+nb+1);
   if nb > 0.5, 
     if min([x.*s;z.*w;tau*kappa])/mu >= (1-beta),
       gamma = 0;
     end;
   else
     if min([x.*s;tau*kappa])/mu >= (1-beta),
       gamma = 0;
     end;
   end;
%
   iter = iter + 1;
 end;
%
% Make an output
%
% b    = b/(1-norm(b));
% c    = c/(1-norm(c));
 s(bindx) = s(bindx) - w;
 if go >= 0,
   x = norbc*x/tau;
   z = norbc*z/tau;
   s = norbc*s/tau;
   y = norbc*y/tau;
   w = norbc*w/tau;
 end
 if nf > .5,
   n = n - nf;
   x(findx) = x(findx) - x(n+1:n+nf);
   x = x(1:n);
   s = s(1:n);
   A = A(:,1:n);
   c = c(1:n);
 end;
 c = c*norbc;
 b = b*norbc;
 if nb > .5
     u=u*norbc;
 end;
 clear norbc ee eez rp rb rd rg i j v
%return
%
%  This program solves linear program:
%
%      minimize    c'*x
%      subject to  A*x=b, 
%                  x_f free, x_p >= 0, u >= x_b >= 0.
%
%  Input 
%      A: sparse constraint matrix
%      b: constraint right-hand column vector
%      c: objective column vector
%
%  Optional Input
%      u: upp bound vector for a set of nonnegative variable x_b, whose
%         indices are represented by bindx.
%      bindx: the index set for those upp-bounded nonnegative variables
%             Default value: []
%      findx: the index set for free variables x_f
%             Default value: []
%      toler: stopping tolerance, the objective value close to the      
%             optimal one in the range of the tolerance. 
%             Default value: 1.0e-6
%      beta : step size, 0 < beta < 1. 
%             Default value: .995
%
%
%  Output
%     x  : optimal solution
%     y  : optimal dual solution (shadow price) for equality constraints
%     w  : optimal dual solution for uppbound constraints
%     s  : optimal dual slack solution
%     zh : (infeasible) primal-dual gap history vs iteration. 
%     
%  Subroutine called : hsdLPsolversub
%
%
%  How to use it? Just type "sphslfbf" and hit the "return".
%
%
%  This program is the MATLAB implementation of the homogeneous and
%  selfdual interior-point algorithm for LP. The result is a
%  strict complementarity solution, i.e., a solution in the relative
%  interior of the LP optimal face.
%
%  Technical References
%
%  Y. Ye, M. J. Todd and S. Mizuno, "An $O(\sqrt{n}L)$-iteration
%  homogeneous and self-dual linear programming algorithm," (1992),
%  to appear in Math. of OR.
%
%  
%  X. Xu, P. Hung and Y. Ye, "A simplified homogeneous and self-dual
%  linear programming algorithm and its implementation," manuscript,
%  Department of Management Sciences, The University of Iowa
%  (Iowa City, IA 52242, 1993).
%
%  Changes from basic version: 
%     
%   1) Adaptively adjust gamma by checking min(x.*s)/mu. If it is
%      bounded from below by constant (1-beta)/10, then gamma=0.
%
%   2) Set eta = 1 to speed up feasibility improvement.
%
%   3) Sparse version.
%
%   4) Handle up bounds x_i <= u_i, for i in bindx.
%
%   5) Handle free variables x_i for i in findx.
%
%   6) cleaned on 11/16/93
