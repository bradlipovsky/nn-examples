function [u,u0] = ssaflowline(p,H)
% SSAFLOWLINE  Computes the velocity from the SSA in the flow-line case
% form:
%   [u,u0] = ssaflowline(p,J,H,b,uleft,initchoice)
% all 6 input arguments are required:
%   p = struct with parameters A,C,g,L,m,n,rho,rhow,secpera
%   J = number of grid subintervals
%   H = thickness (m), a length J+1 column vector
%   b = bed elevation (m), same size as H
%   uleft = velocity boundary value at x=0 end of domain
%   initchoice = 1,2,3 chooses initial guess scheme;
%                use 1 for ice shelves and 2 for ice streams
% outputs:
%   u = the velocity solution (m s-1), a length J+1 column vector
%   u0 = the velocity initial guess (m s-1), same size as u
% does "Picard" iteration to solve nonlinear SSA problem
% calls SSAINIT to get initial guess for solving the SSA
% calls FLOWLINE to solve inner linear PDE boundary value problem
% example:  TESTSHELF

%if nargin ~= 6, error('exactly 6 input arguments required'), end

dx = p.L / p.J;
x = (0:dx:p.L)';
xstag = (dx/2:dx:p.L-dx/2)';

% create parts of PDE problem to solve; see flowline.m
% coefficient for dragging:
alpha = p.C * ones(size(x));
h = H + p.b;
hx = regslope(dx,h);
% driving stress is right-hand side:
beta = p.rho * p.g * H .* hx;
% value in calving front stress boundary condition:
gamma = ( 0.25 * p.A^(1/p.n) * (1 - p.rho/p.rhow) *...
          p.rho * p.g * H(end) )^p.n;

u0 = ssainit(p,x,beta,gamma,p.initchoice);
u = u0;

% "outer" iteration for solution of nonlinear equation;
%   "Picard" iteration
Hstag = stagav(H);
tol = 1.0e-10;  % m s-1; = 3 mm a-1 velocity error
maxdiff = Inf;
while maxdiff > tol
  % find coefficient W(x) on staggered grid using "old" u
  uxstag = stagslope(dx,u);
  W = 2 * p.A^(-1/p.n) * Hstag .* (abs(uxstag)).^((1/p.n)-1);

  % solve problem for this W
  unew = flowline(p.L,p.J,gamma,W,alpha,beta,p.uleft);
  maxdiff = max(abs(unew-u));
  u = unew;
  fprintf('.')
end
fprintf('\n')

  function fav = stagav(f)
  % average regular grid values onto staggered grid
  fav = 0.5 * (f(1:end-1) + f(2:end));

  %function slope = regslope(dx,f)
  % compute regular grid values of slope  f_x = f'
  % if length(h)=J+1, returns vector with length J+1
  %J = length(f) - 1;
  %slope = [(f(2)-f(1))/dx; (f(3:J+1)-f(1:J-1))/(2*dx); (f(J+1)-f(J))/dx];

  function slope = stagslope(dx,f)
  % compute staggered grid values of slope  f_x = f'
  % if length(h)=J+1, returns vector with length J
  slope = (f(2:end) - f(1:end-1)) / dx; 
