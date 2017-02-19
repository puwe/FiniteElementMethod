function [T0 M] = Boundary2(t0,data)
[m n] = size(data.node);
T0 = 20.*ones(m,1);
M = zeros(m,1);
for j = 1:m;
    d = abs(data.node(j,2))^2+abs(data.node(j,3))^2+abs(data.node(j,4))^2;
    if abs(d-1)<=0.001;
        T0(j) = t0;
        M(j) = 1;
    end
end