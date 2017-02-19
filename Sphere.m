function [U] = Sphere(U0,V0,a,N,r0,r,t)
E = 0;
for n = 1:N;
    E = E + 1/n*(-1)^n*exp(-n^2*pi^2*a^2*t/r0^2)*sin(n*pi*r/r0);
end
U = U0 + (2*(U0-V0)*r0)/(pi*r)*E;