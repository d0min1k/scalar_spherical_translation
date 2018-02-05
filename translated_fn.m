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
  N_nu =  15;

  y = 0;

  for nu = 0 : N_nu
    for mu = - nu : 1 : nu
      for p = abs(n - nu) : 2 : n + nu % according to Seymour Stein
%       for p = abs(n - nu) : 1 : n + nu 
%       for p = 0 : 1 : n + nu % seems promising
%       for p = 0 : 1 : n + nu % seems promising
%       for p = -N_nu  : 1 : N_nu 

        % Checking
        a = a_coef(m, n, -mu, nu, p);

%         if (a == 0) | (p < m - mu) | (nu < mu)
        if (a == 0) | (a==Inf)
          continue;
%           break;
        end


        % Summing terms
        y_new =  power(-1, mu).*power( 1i, nu + p - n ).*(2*nu + 1)...
        .*a.*sphbes(nu, R1.*k) ...
        .*Z_fn(p, k.*r0) ...
        .*legendreP2( nu, mu, cos(TH1) ) ...
        .*legendreP2( p, m - mu, cos(th0) ) ...
        .*exp( 1i.*( m - mu ).*ph0 + 1i.*mu.*PH1 );
    
%         figure()
%         clf
%         pcolor(r.*sin(th).*cos(ph),r.*sin(th).*sin(ph) ,abs(y_new)); shading interp; hold on
%         scatter(1.5,0, 'r*'); hold off; colorbar;
%         pause
        
        y = y + y_new;
        
        if any(  any( isnan(y) )   )
            disp('problem')
        end
    
      end
    end
  end


end  % translated_fn
