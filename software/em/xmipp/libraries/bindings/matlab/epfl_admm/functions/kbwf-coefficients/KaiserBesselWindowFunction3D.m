function  Y  =  KaiserBesselWindowFunction3D(x,alpha,a,m)

[X1,X2,X3]  =   ndgrid(x);
X           =   sqrt(X1.^2+X2.^2+X3.^2)    ;

Y           =   zeros(size(X));
z           =   alpha*sqrt(1-(abs(X)/a).^2);

Y(abs(X)<a) =  (1/besseli(m,alpha))*((z(abs(X)<a)/alpha).^m).*...
                           besseli(m,z(abs(X)<a));
