%% Init

n=1;
m=1;

nu=1;
mu=1;


x = linspace(0,10,100);


%% Testing Legendre translation
y1 = legendreP2(n,m,x).*legendreP2(nu,mu,x);

y2=0;
for p=0:20
  y2 = y2 + a_coef( m, n, mu, nu, p ).*legendreP2(p, m + mu, x);
end


%% Plots [Passed]

plot(x, y1);
hold on
plot(x-0.01, y2);
xlabel('x')
