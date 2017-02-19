%function [es ei] = Error(data, T)
[N E R] = BToS(data);
[m d] = size(T);
S = zeros(m,d);
for i = 1:m
    for j = 1:d
        r = R(i);
        t = 0.05*(j-1);
        u = Sphere(50,20,1,2000,1,r,t);
        S(i,j) = u;
    end
end
es = max(max(T-S));
ei = min(min(T-S));
