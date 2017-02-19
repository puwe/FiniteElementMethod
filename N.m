function [k f hh pQ pq ph] = N(xi,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE)

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
pQ = (a*b*c).*rho*Q.*N;
hh = (a*b)*h.*N*N';%+z heat convection
pq = (a*b)*q.*N; %-z heat conduction
ph = (a*b)*h*TE.*N; %+z heat convection

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
end
