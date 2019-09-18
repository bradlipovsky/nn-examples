function [u,time,cond] = flowline(L,J,gamma,W,alpha,beta,V0)
% FLOWLINE Compute finite difference solution to elliptic PDE problem
%   (W(x) u_x)_x - alpha(x) u = beta(x)     on  0 < x < L
%   u(0) = V0,   u_x(L) = gamma
% requires:
%   W     a vector of length J (staggered values W_{j+1/2}, j=1:J)
%   alpha a vector of length J+1 (reg values alpha_j, j=1:J+1)
%   beta  like alpha
% returns:
%   u     a solution vector (reg values alpha_j, j=1:J+1)
% optionally returns computation time in seconds
% optionally returns 2-norm condition number of A, BUT the
%   evaluation of cond(A) is very slow compared to solve time
% examples: see TESTFLOWLINE, CONVANALYSIS, SSAFLOWLINE

if nargout>=2, tic; end % start timer

dx = L / J;

rhs = dx^2 * beta(:); % a (J+1) length column vector
rhs(1) = V0;
rhs(J+1) = rhs(J+1) - 2 * gamma * dx * W(J);

A = sparse(J+1,J+1);  % allocates no space yet
A(1,1) = 1.0;
for j=2:J
  A(j,j-1:j+1) = [ W(j-1), -(W(j-1) + W(j) + alpha(j) * dx^2), W(j) ];
end
A(J+1,[J J+1]) = [ 2 * W(J), -(2 * W(J) + alpha(J+1) * dx^2) ];
%spy(A)   % shape of A is tridiagonal
%full(A)  % to see values in *small* J cases

% scale A by rows
scale = max(abs(A),[],2);    % column vector of row maximums
A = spdiags(1./scale,0,size(A,1),size(A,2)) * A;  % divide rows of A
rhs = rhs ./ scale;

% solve by Matlab/Octave default methods
u = A \ rhs;

if nargout>=2, time = toc; end
if nargout>=3, cond = cond(A); end

