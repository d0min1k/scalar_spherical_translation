k = 1;
init;

n=1;
m=1;

%% ------ Using addition theorem for spherical scalar wave function -------

fn = @(r, th, ph) Z_fn(n, k.*r).*legendreP2( n, m, cos(th) ).*exp(1i.*m.*ph);


%% ------------------------------- DO plots -------------------------------
x_ = linspace(-1, 1, 200);
[x, y] = meshgrid(x_, x_);
z = 0;

R =  sqrt(x.^2 + y.^2 + z.^2);
TH = pi/2;
PH = atan2(y, x);

R1 =  sqrt((x-1).^2 + y.^2 + z.^2);
TH1 = pi/2;
PH1 = atan2(y, (x-1) );

%%
figure()
I1 = abs( fn(R1, TH1, PH1) );
pcolor( x+1, y, I1);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])

%%
figure()
I2 = abs( translated_fn( Z_fn, R, TH, PH, r0, th0, ph0, n, m ) );
pcolor(x, y, I2);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])


%%
hold on
scatter(1,0, '*r')
