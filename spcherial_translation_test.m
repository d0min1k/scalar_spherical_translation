k = 1;
init;

n=1;
m=1;

%% ------ Using addition theorem for spherical scalar wave function -------

fn = @(r, th, ph) Z_fn(n, k.*r).*legendreP2( n, m, cos(th) ).*exp(1i.*m.*ph);


%% ------------------------------- DO plots -------------------------------
x_ = linspace(-3, 3, 200);
[x, y] = meshgrid(x_, x_);
z = 0;

R =  sqrt(x.^2 + y.^2 + z.^2);
TH = pi/2;
PH = atan2(y, x);

%%
figure()
I1 = abs(fn(R, TH, PH));
pcolor(x,y,I1);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])

%%
figure()
I2 = abs( translated_fn( Z_fn, R, TH, PH, r0, th0, ph0, n, m ) );
pcolor(x,y,I2);
shading interp;
xlabel('x'); ylabel('y');
colorbar;
pbaspect([1 1 1])


%%
hold on 
scatter(1,0, '*r')