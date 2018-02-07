%%%%%%%%%%%%%%%%%%

% Graf addition theorem

%%%%%%%%%%%%%%%%%%

x=(-1:.01:1)*3;
[X,Y]=meshgrid(x,x);

r1=sqrt(X.^2+Y.^2);
th1=atan2(Y,X);
%reference frame 1 - r1, th1

r0=1;
th0=0;
%displacement between reference frames - r0, th0

r2=sqrt((X-r0).^2+Y.^2);
th2=atan2(Y,X-r0);
%reference frame 2 - r2, th2

n=1;

% lhs=besselj(n,r1).*exp(1i*n*th1);
lhs=bessely(n,r1).*exp(1i*n*th1);
figure, contourf(X,Y,abs(lhs),30), colorbar

rhs=0;
rhs_in=0;
rhs_out=0;
mm=10;
for m=(-mm:1:mm)
%     rhs=rhs+besselj(m+n,r0).*besselj(m,r2).*exp(1i*m*(pi-th2+th0)).*exp(1i*n*th0);
    rhs_in=rhs_in+conj(bessely(m+n,r0).*besselj(m,r2).*exp(1i*m*(pi-th2+th0))).*exp(1i*n*th0); %only works for r2<r0 for non-besselj
    rhs_out=rhs_out+besselj(m+n,r0).*bessely(m,r2).*exp(1i*m*(pi-th2+th0)).*exp(1i*n*th0); %only works for r2>0 for other bessel functions
end

% for non besselj merge inner and outer regions
rhs_out(r2.^2<r0)=-rhs_in(r2.^2<r0);

% figure, contourf(X,Y,abs(rhs),30), colorbar
figure, contourf(X,Y,abs(rhs_out),30), colorbar


