function [ a ] = a_coef( m, n, mu, nu, p )
% a_coef: coeficients for translational addition theorems

  % conditions
  c1 = and( (n + m >= 0)      , (n - m >= 0)    );
  c2 = and( (nu + mu >= 0)    , (nu - mu >= 0)    );
  c3 = and( (p - m - mu >= 0) , (p + m + mu >= 0)    );

  if all([c1 , c2, c3])

    % Computing factors
    f1 = factorial(n + m)      / factorial(n - m);
    f2 = factorial(nu + mu)    / factorial(nu - mu);
    f3 = factorial(p - m - mu) / factorial(p + m + mu);

    % Usage: wigner3j( j1, j2, j3, m1, m2, m3 )
    f4 = wigner3j( n, nu, p, 0, 0, 0);
    f5 = wigner3j( n, nu, p, m, mu, -m -mu);

    % the addition theorem coeficient
    a = power(-1, m + mu ).*(2*p + 1).*sqrt(f1*f2*f3)*f4*f5;

  else
    a = 0;
  end

end  % a_coef
