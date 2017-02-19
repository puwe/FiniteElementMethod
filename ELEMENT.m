function [K,C] = ELEMENT(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,kx,ky,kz,rho,e)
%
%kx=1;
%ky=1;
%kz=1;
%rho = 0.1;
%e = 2;
%xi = 0;
%eta = 0;
%zeta = 0;

xc = (x1+x2+x3+x4+x5+x6+x7+x8)/8;
yc = (y1+y2+y3+y4+y5+y6+y7+y8)/8;
zc = (z1+z2+z3+z4+z5+z6+z7+z8)/8;

%xi1 = (x1-xc)/abs(x1-xc);
%xi2 = (x2-xc)/abs(x2-xc);
%xi3 = (x3-xc)/abs(x3-xc);
%xi4 = (x4-xc)/abs(x4-xc);
%xi5 = (x5-xc)/abs(x5-xc);
%xi6 = (x6-xc)/abs(x6-xc);
%xi7 = (x7-xc)/abs(x7-xc);
%xi8 = (x8-xc)/abs(x8-xc);
a = 1/8*(abs(x1-xc)+abs(x2-xc)+abs(x3-xc)+abs(x4-xc)+abs(x5-xc)+abs(x6-xc)+abs(x7-xc)+abs(x8-xc));

%eta1 = (y1-yc)/abs(y1-yc);
%eta2 = (y2-yc)/abs(y2-yc);
%eta3 = (y3-yc)/abs(y3-yc);
%eta4 = (y4-yc)/abs(y4-yc);
%eta5 = (y5-yc)/abs(y5-yc);
%eta6 = (y6-yc)/abs(y6-yc);
%eta7 = (y7-yc)/abs(y7-yc);
%eta8 = (y8-yc)/abs(y8-yc);
b = 1/8*(abs(y1-yc)+abs(y2-yc)+abs(y3-yc)+abs(y4-yc)+abs(y5-yc)+abs(y6-yc)+abs(y7-yc)+abs(y8-yc));

%zeta1 = (z1-zc)/abs(z1-zc);
%zeta2 = (z2-zc)/abs(z2-zc);
%zeta3 = (z3-zc)/abs(z3-zc);
%zeta4 = (z4-zc)/abs(z4-zc);
%zeta5 = (z5-zc)/abs(z5-zc);
%zeta6 = (z6-zc)/abs(z6-zc);
%zeta7 = (z7-zc)/abs(z7-zc);
%zeta8 = (z8-zc)/abs(z8-zc);
c = 1/8*(abs(z1-zc)+abs(z2-zc)+abs(z3-zc)+abs(z4-zc)+abs(z5-zc)+abs(z6-zc)+abs(z7-zc)+abs(z8-zc));

N1 = 1/8*(1-xi)*(1+eta)*(1+zeta);
N2 = 1/8*(1-xi)*(1-eta)*(1+zeta);
N3 = 1/8*(1-xi)*(1-eta)*(1-zeta);
N4 = 1/8*(1-xi)*(1+eta)*(1-zeta);
N5 = 1/8*(1+xi)*(1+eta)*(1+zeta);
N6 = 1/8*(1+xi)*(1-eta)*(1+zeta);
N7 = 1/8*(1+xi)*(1-eta)*(1-zeta);
N8 = 1/8*(1+xi)*(1+eta)*(1-zeta);

N = [N1, N2, N3, N4, N5, N6, N7, N8]';
f = (a*b*c).*(rho*e).*N*N';

dN1dxi = 1/8*(-1)*(1+eta)*(1+zeta);
dN2dxi = 1/8*(-1)*(1-eta)*(1+zeta);
dN3dxi = 1/8*(-1)*(1-eta)*(1-zeta);
dN4dxi = 1/8*(-1)*(1+eta)*(1-zeta);
dN5dxi = 1/8*(1)*(1+eta)*(1+zeta);
dN6dxi = 1/8*(1)*(1-eta)*(1+zeta);
dN7dxi = 1/8*(1)*(1-eta)*(1-zeta);
dN8dxi = 1/8*(1)*(1+eta)*(1-zeta);

dN1deta = 1/8*(1-xi)*(1)*(1+zeta);
dN2deta = 1/8*(1-xi)*(-1)*(1+zeta);
dN3deta = 1/8*(1-xi)*(-1)*(1-zeta);
dN4deta = 1/8*(1-xi)*(1)*(1-zeta);
dN5deta = 1/8*(1+xi)*(1)*(1+zeta);
dN6deta = 1/8*(1+xi)*(-1)*(1+zeta);
dN7deta = 1/8*(1+xi)*(-1)*(1-zeta);
dN8deta = 1/8*(1+xi)*(1)*(1-zeta);

dN1dzeta = 1/8*(1-xi)*(1+eta)*(1);
dN2dzeta = 1/8*(1-xi)*(1-eta)*(1);
dN3dzeta = 1/8*(1-xi)*(1-eta)*(-1);
dN4dzeta = 1/8*(1-xi)*(1+eta)*(-1);
dN5dzeta = 1/8*(1+xi)*(1+eta)*(1);
dN6dzeta = 1/8*(1+xi)*(1-eta)*(1);
dN7dzeta = 1/8*(1+xi)*(1-eta)*(-1);
dN8dzeta = 1/8*(1+xi)*(1+eta)*(-1);

dNdxi = [dN1dxi, dN2dxi, dN3dxi, dN4dxi, dN5dxi, dN6dxi, dN7dxi, dN8dxi]';
dNdeta = [dN1deta, dN2deta, dN3deta, dN4deta, dN5deta, dN6deta, dN7deta, dN8deta]';
dNdzeta = [dN1dzeta, dN2dzeta, dN3dzeta, dN4dzeta, dN5dzeta, dN6dzeta, dN7dzeta, dN8dzeta]';

dxidx = 1/a;
dNdx = dNdxi.*dxidx;
detady = 1/b;
dNdy = dNdeta.*detady;
dzetadz = 1/c;
dNdz = dNdzeta.*dzetadz;

Kx = kx.*dNdx*dNdx';
Ky = ky.*dNdy*dNdy';
Kz = kz.*dNdz*dNdz';

k = (a*b*c).*(Kx+Ky+Kz);

%Irons:
%A1F(0,0,0)+B6{F(-b,0,0)+F(b,0,0)+F(0,-b,0)+F(0,b,0)+F(0,0,-b)+F(0,0,b)}+
%C8{F(-c,-c,-c)+F(c,-c,-c)+F(-c,c,-c)+F(-c,-c,c)+F(c,c,-c)+F(c,-c,c)+F(-c,c,c)+F(c,c,c)}
%14Points:B6=0.886426593;C8=0.335180055; b=0.795822426;c=0.758786911;

B6=0.886426593;
C8=0.335180055;
Ib=0.795822426;
Ic=0.758786911;
K = zeros(8);
C = zeros(8);
IB = [
    -Ib, 0,  0;
     Ib, 0,  0;
     0,  -Ib,0;
     0,  Ib, 0;
     0,  0,  -Ib;
     0,  0,  Ib;
    ];
IC = Ic.*[
    -1,  -1, -1;
    1,   -1, -1;
    -1,  1,  -1;
    -1,  -1, 1;
    1,   1,  -1;
    1,   -1, 1;
    -1,  1,  1;
    1,   1,  1;
    ];
for i = 1:6;
    xi = IB(i,1);
    eta = IB(i,2);
    zeta = IB(i,3);
    K = K + B6.*k;
    C = C + B6.*f;
end
for j = 1:8;
    xi = IC(j,1);
    eta = IC(j,2);
    zeta = IC(j,3);
    K = K + C8.*k;
    C = C + C8.*f;
end
