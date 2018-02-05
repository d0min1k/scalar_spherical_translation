clear all
k = 1;
init;

n=2;
m=2;

%% ------ Using addition theorem for spherical scalar wave function -------

fn = @(r, th, ph) Z_fn(n, k.*r).*legendreP2( n, m, cos(th) ).*exp(1i.*m.*ph);


%% ------------------------------- DO plots -------------------------------
T = 10;
x_ = linspace(-T, T, 200);
[x, y] = meshgrid(x_, x_);
z = 0;

R =  sqrt(x.^2 + y.^2 + z.^2);
TH = pi/2;
PH = atan2(y, x);

R1 =  sqrt( (x-r0*cos(ph0)).^2 + (y-r0*sin(ph0)).^2 + z.^2);
TH1 = pi/2;
PH1 = atan2(y, (x-1) );

%% ---------------------------------- I0 -----------------------------------
figure()
I0 = abs( fn(R, TH, PH) );
pcolor( x, y, I0);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])
title('Initial Wave Position')


%% ---------------------------------- I1 -----------------------------------
figure()
I1 = abs( fn(R1, TH1, PH1) );
pcolor( x, y, I1);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])
title('Wanted Wave Position')

hold on
scatter(r0*sin(th0)*cos(ph0), r0*sin(th0)*sin(ph0), '*r')

%% ---------------------------------- I2 -----------------------------------

% I0 = abs( fn(R, TH, PH) );
figure()
I2 = abs( translated_fn( Z_fn, R, TH, PH, r0, th0, ph0, n, m ) );
pcolor(x, y, I2);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])
title('Translated Scalar Spherical Wave')

% hold on
% scatter(r0*sin(th0)*cos(ph0), r0*sin(th0)*sin(ph0), '*r')


%% Difference I2-I1
pcolor(x,y ,I2 - I1)
shading interp
colorbar

%% Difference I2-I0
pcolor(x,y ,I2 - I0)
shading interp
colorbar