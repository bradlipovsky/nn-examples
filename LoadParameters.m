function p = LoadParameters

p.alpha = 1e-3;

p.g = 9.81; p.rho = 916; p.rhow = 1e3;
p.secpera = 31556926;
p.A = 1.0e-13 / p.secpera;
p.n = 3;
p.L = 10e3;
p.Ac = 1e-3;
p.J=300;
dx = p.L / p.J;
p.x = (0:dx:p.L)';
p.b = -p.alpha * p.x;
p.C = 1e4;
p.uleft = 1;
p.initchoice = 2;
p.m=1;

% Good starting geometry
p.Hinit = p.C*p.uleft / (p.rho*p.g*p.alpha);
p.H0 = p.Hinit * (ones(p.J+1,1) -  (1:p.J+1)' * p.alpha);