function [t,x,y]=Kramers_euler_maruyama (omega,gamma, D, x0, y0, sigma0, tmax, P, n, r)

  if ( nargin < 1 )
    omega = 1.0;
  end

  if ( nargin < 2 )
    gamma = 2.0;
  end

  if ( nargin < 3 )
    D = [0.1 0; 0, 0.1];
  end

  if ( nargin < 4 )
    x0 = 2.0;
  end

  if ( nargin < 5 )
    y0 = 3.0;
  end 

  if ( nargin < 6 )
    sigma0 = [0.1 0; 0, 0.1];
  end 
  
  if ( nargin < 7 )
    tmax = 3.0;
  end 

  if ( nargin < 8 )
    P= 100;
  end
  
  if ( nargin < 9 )
    n = 10000;
  end

  if ( nargin < 10 )
    r = 50;
  end
%
%  Set time steps.
%
  dt_large = tmax / n;
  dt_small = tmax / n / r;
%
%  Carry out the Euler-Maruyama approximate integration process.
%
  t = linspace ( 0, tmax, n + 1 );
  x = zeros ( P, n + 1 );
  y = zeros ( P, n + 1 );

  x(:,1) = normrnd(x0,sigma0(1,1),P,1);
  y(:,1) = normrnd(y0,sigma0(2,2),P,1);
  for j = 1 : n
    dwx = sqrt ( dt_small ) * randn ( P, r );
    dwy = sqrt ( dt_small ) * randn ( P, r );
    x(:,j+1) = x(:,j) + dt_large *y(:,j) + sqrt(D(1,1)) * sum ( dwx(:,1:r), 2 );
    y(:,j+1) = y(:,j) + dt_large *( -(omega^2)*x(:,j)-gamma*y(:,j) ) + sqrt(D(2,2)) * sum ( dwy(:,1:r), 2 );
  end

  return
end
