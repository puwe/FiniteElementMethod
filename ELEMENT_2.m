function [K,C,H,PQ,Pq,PH] = ELEMENT_2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,kx,ky,kz,rho,e,h,Q,q,TE)
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

%Irons:
%A1F(0,0,0)+B6{F(-b,0,0)+F(b,0,0)+F(0,-b,0)+F(0,b,0)+F(0,0,-b)+F(0,0,b)}+
%C8{F(-c,-c,-c)+F(c,-c,-c)+F(-c,c,-c)+F(-c,-c,c)+F(c,c,-c)+F(c,-c,c)+F(-c,c,c)+F(c,c,c)}
%D12{F(-d,-d,0)+F(d,-d,0)+F(-d,0,-d)+F(d,0,-d)+F(0,-d,-d)+...
%14Points:B6=0.886426593;C8=0.335180055; b=0.795822426;c=0.758786911;
%27Points:A1=0.788073483;B6=0.499369002;C8=0.478508449;D12=0.032303742;
%         b=0.848418011;c=0.652816472;d=1.106412899
A1=0.788073483;
B6=0.499369002;
C8=0.478508449;
D12=0.032303742;
Ib=0.848418011;
Ic=0.652816472;
Id=1.106412899;
K = zeros(8);
C = zeros(8);
PQ = zeros(8,1);
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
ID = Id.*[
    -1, -1, 0;
    1,  -1, 0;
    -1, 1,  0;
    1,  1,  0;
    -1, 0, -1;
    -1, 0, 1;
    1,  0, -1;
    1,  0, 1;
    0,  1, -1;
    0,  1, 1;
    0,  -1, -1;
    0,  -1, 1;
    ];
[k,f,~,pQ,~,~] = N(0,0,0,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
K = K + A1.*k;
C = C + A1.*f;
PQ = PQ + A1.*pQ;
for i = 1:6;
    xi = IB(i,1);
    eta = IB(i,2);
    zeta = IB(i,3);
    [k,f,~,pQ,~,~] = N(xi,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
    K = K + B6.*k;
    C = C + B6.*f;
    PQ = PQ + B6.*pQ;
end
for j = 1:8;
    xi = IC(j,1);
    eta = IC(j,2);
    zeta = IC(j,3);
    [k,f,~,pQ,~,~] = N(xi,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
    K = K + C8.*k;
    C = C + C8.*f;
    PQ = PQ + C8.*pQ;
end
for z = 1:12;
    xi = ID(z,1);
    eta = ID(z,2);
    zeta = ID(z,3);
    [k,f,~,pQ,~,~] = N(xi,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
    K = K + D12.*k;
    C = C + D12.*f;
    PQ = PQ + D12.*pQ;
end
% Gauss Integral
% x1 = +-0.861136311594053;H1 = 0.347854845137454
% x2 = +-0.339981043584856;H2 = 0.652145154862546
G1 = 0.861136311594053;
G2 = 0.339981043584856;
H1 = 0.347854845137454;
H2 = 0.652145154862546;
H = zeros(8);
PH = zeros(8,1);
Pq = zeros(8,1);
for eta = [G1,-G1,G2,-G2]
    for xi = [G1,-G1,G2,-G2]
        [~,~,hh,~,~,ph] = N(xi,eta,1,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %+z heat convection
        %[~,~,~,~,~,ph] = N(xi,eta,1,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
        [~,~,hh2,~,pq,ph2] = N(xi,eta,-1,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %-z heat conduction
        if abs(eta) == G1 && abs(xi) == G1;
            H = H + H1*H1.*hh + H1*H1.*hh2;
            PH = PH + H1*H1.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H1.*pq;
        end
        if abs(eta) == G2 && abs(xi) == G2;
            H = H + H2*H2.*hh + H1*H1.*hh2;
            PH = PH + H2*H2.*ph + H1*H1.*ph2;
            Pq = Pq + H2*H2.*pq;
        end
        if abs(eta) == G1 && abs(xi) == G2;
            H = H + H1*H2.*hh + H1*H1.*hh2;
            PH = PH + H1*H2.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
        if abs(eta) == G2 && abs(xi) == G1;
            H = H + H1*H2.*hh + H1*H1.*hh2;
            PH = PH + H1*H2.*ph +  H1*H1.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
    end
end
for zeta = [G1,-G1,G2,-G2]
    for xi = [G1,-G1,G2,-G2]
        [~,~,hh,~,~,ph] = N(xi,1,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %+z heat convection
        %[~,~,~,~,~,ph] = N(xi,eta,1,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
        [~,~,hh2,~,pq,ph2] = N(xi,-1,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %-z heat conduction
        if abs(zeta) == G1 && abs(xi) == G1;
            H = H + H1*H1.*hh + H1*H1.*hh2;
            PH = PH + H1*H1.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H1.*pq;
        end
        if abs(zeta) == G2 && abs(xi) == G2;
            H = H + H2*H2.*hh + H1*H1.*hh2;
            PH = PH + H2*H2.*ph + H1*H1.*ph2;
            Pq = Pq + H2*H2.*pq;
        end
        if abs(zeta) == G1 && abs(xi) == G2;
            H = H + H1*H2.*hh + H1*H1.*hh2;
            PH = PH + H1*H2.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
        if abs(zeta) == G2 && abs(xi) == G1;
            H = H + H1*H2.*hh + H1*H1.*hh2;
            PH = PH + H1*H2.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
    end
end 
for zeta = [G1,-G1,G2,-G2]
    for eta = [G1,-G1,G2,-G2]
        [~,~,hh,~,~,ph] = N(1,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %+z heat convection
        %[~,~,~,~,~,ph] = N(xi,eta,1,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE);
        [~,~,hh2,~,pq,ph2] = N(-1,eta,zeta,a,b,c,kx,ky,kz,rho,e,h,Q,q,TE); %-z heat conduction
        if abs(zeta) == G1 && abs(eta) == G1;
            H = H + H1*H1.*hh + H1*H1.*hh2;
            PH = PH + H1*H1.*ph + H1*H1.*ph2;
            Pq = Pq + H1*H1.*pq;
        end
        if abs(zeta) == G2 && abs(eta) == G2;
            H = H + H2*H2.*hh + H2*H2.*hh2;
            PH = PH + H2*H2.*ph + H2*H2.*ph2;
            Pq = Pq + H2*H2.*pq;
        end
        if abs(zeta) == G1 && abs(eta) == G2;
            H = H + H1*H2.*hh + H1*H2.*hh2;
            PH = PH + H1*H2.*ph + H1*H2.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
        if abs(zeta) == G2 && abs(eta) == G1;
            H = H + H1*H2.*hh + H1*H2.*hh2;
            PH = PH + H1*H2.*ph + H1*H2.*ph2;
            Pq = Pq + H1*H2.*pq;
        end
    end
end
        
