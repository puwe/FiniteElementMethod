function [Node Element Radius] = BToS(data)
[m ~] = size(data.node);
x = data.node(:,2);
y = data.node(:,3);
z = data.node(:,4);
d = sqrt(x.*x+y.*y+z.*z);
XX = zeros(m,1);
YY = zeros(m,1);
ZZ = zeros(m,1);
N = zeros(m,1);
Radius = zeros(m,1);
for i = 1:m;
    xx = x(i);
    yy = y(i);
    zz = z(i);
    dd = d(i);
    r = max(max(abs(xx),abs(yy)),abs(zz));
    Radius(i) = r;
    if dd == 0;
        XX(i) = 0;
        YY(i) = 0;
        ZZ(i) = 0;
    else
        X = r*xx/dd;
        Y = r*yy/dd;
        Z = r*zz/dd;
        XX(i,1) = X;
        YY(i,1) = Y;
        ZZ(i,1) = Z;
    end
    N(i,1) = i;
end
Node = cat(2,N,XX,YY,ZZ);
Element = data.element;