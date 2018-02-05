%% -------------------------------Initialize-------------------------------

% coordinates of system 1 in system 0
th0 = pi/2;
ph0 = 0;
r0 = 1;

% new coordinates in terms of old spherical coordinates
x1 = @(r, th, ph) r.*sin(th).*cos(ph) - r0.*sin(th0).*cos(ph0);
y1 = @(r, th, ph) r.*sin(th).*sin(ph) - r0.*sin(th0).*sin(ph0);
z1 = @(r, th, ph) r.*cos(th) - r0.*cos(th0);

r1  = @(r, th, ph) sqrt( x1(r, th, ph).^2 + y1(r, th, ph).^2 + z1(r, th, ph).^2 );
th1 = @(r, th, ph) acos( z1(r, th, ph) ./ (r1(r, th, ph) + eps) ); % TODO: May include singularity.
% th1 = @(r, th, ph)  pi / 2;
ph1 = @(r, th, ph) atan2( y1(r, th, ph), x1(r, th, ph) );


% --------------------------- The wave function ---------------------------
Z_fn = @(n, x) sphbes(n, x);
