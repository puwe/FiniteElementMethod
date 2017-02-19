function [T0 M] = Boundary(t0,data)
[m n] = size(data.node);
T0 = 20.*ones(m,1);
M = zeros(m,1);
for j = 1:m;
    if max(abs(data.node(j,2:4))) == 1;
        T0(j) = t0;
        M(j) = 1;
    end
end
%theta = 0.5;
%delta = 0.5;
%KB = CC./delta + theta.*KK;
%QB = (CC./delta - (1-theta).*KK)*T0;
%T1 = inv(KB)*QB;