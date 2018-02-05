function [ y ] = translated_fn( Z_fn, r, th, ph, r0, th0, ph0, n, m )
% translated_fn: Translated function
%  @Z(p, k.*r0)


  k = 1;

  % ----------------------------- Coordinates -----------------------------
  init;

  % --------------------------- new coordinates ---------------------------

  PH1 = ph1(r, th, ph);
  TH1 = th1(r, th, ph);
  R1  = r1(r, th, ph);

% ------------------------------- Summation -------------------------------
  N_nu =  20;

  y = 0;

  for nu = 0 : N_nu
    for mu = - nu : 1 : nu
      for p = abs(n - nu) : 2 : n + nu
      % for p = 1 : 1 : n + nu
      % for p = -N_nu : 1 : N_nu

        % Checking
        a = a_coef(m, n, -mu, nu, p);

        if (a == 0) | (p < m - mu) | (nu < mu)
          break;
        end


        % Summing terms
        y = y + power(-1, mu).*power( 1i, nu + p - n ).*(2*nu + 1)...
        .*a.*sphbes(nu, R1.*k) ...
        .*Z_fn(p, k.*r0) ...
        .*legendreP2( nu, mu, cos(th0) ) ...
        .*legendreP2( p, m - mu, cos(TH1) ) ...
        .*exp( 1i.*( m - mu ).*ph0 + 1i.*mu.*PH1 );
      end
    end
  end


end  % translated_fn
