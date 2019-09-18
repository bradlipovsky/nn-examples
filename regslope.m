  function slope = regslope(dx,f)
  % compute regular grid values of slope  f_x = f'
  % if length(h)=J+1, returns vector with length J+1
  J = length(f) - 1;
  slope = [(f(2)-f(1))/dx; (f(3:J+1)-f(1:J-1))/(2*dx); (f(J+1)-f(J))/dx];
