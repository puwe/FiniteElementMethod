%importdata('data1.txt');
[Data.node Data.element] = BToS(data);
[m n] = size(Data.node); % Node Number; 
kx = 1;
ky = 1;
kz = 1;
rho = 1;
c = 1;
h = 0;
TE = 50;
Q = 0;
q = 0;
t0 = 50;
[T0 Num] = Boundary2(t0,Data);
n = 60;
[KK CC PP] = ASSEMBLE(Data,kx,ky,kz,rho,c,h,Q,q,TE);
for i = 1:m
    if Num(i) == 1;
        PP = PP - KK(:,i).*t0;
    end
end
for d = 1:m
    if Num(d) == 1;
        KK(:,d) = 0;
        KK(d,:) = 0;
        KK(d,d) = 1;
        CC(:,d) = 0;
        CC(d,:) = 0;
        CC(d,d) = 1;
        PP(d) = t0;
    end
end
% for i = 1:m
%     if Num(i) == 1;
%         KK(i,i) = KK(i,i)*10^22;
%         PP(i) = KK(i,i)*10^22*t0;
%     end
% end
    
T = zeros(m,n);
T(:,1) = T0;
theta = 1;
delta = 0.05;
KB = CC./delta + theta.*KK;
PP1 = PP;
for j = 1:n;
    QB = (CC./delta - (1-theta).*KK)*T0 + (1-theta).*PP1 + theta.*PP;
    T1 = KB\QB;
    T(:,j+1) = T1;
%     for k = 1:m;
%         if Num(k) == 1;
%             T1(k) = t0;
%         end
%     end
    T0 = T1;
end
t = 0:delta:n*delta;
plot(t,T)
grid on

% T = zeros(m,n);
% T(1:m,1) = T0;
% hh = 0.0002;
% y0 = T0;
% T(:,1) = y0;
% for k = 1:m
%         if Num(k) == 1;
%             K(k,:) = 0;
%             K(k,k) = 1;
%             P(k,1) = t0;
%         end
% end
% P = PP;
% K = KK;
% C = CC;
% y0 = T0;
% hh = 0.00005;
% for j = 1:n
%    k1 = C\(P-K*y0);
%    k2 = C\(P-K*(y0+1/2*hh.*k1));
%    k3 = C\(P-K*(y0+1/2*hh.*k2));
%    k4 = C\(P-K*(y0+hh.*k3));
%    y1 = y0 + hh/6.*(k1+2.*k2+2.*k3+k4);
%    T(:,j+1) = y1;
%    y0 = y1;
% end

% TT = zeros(m,n);
% TT(:,1) = T0;
% for j = 1:n
%     t = j*hh;
%     phi = PHI(K,C,t,T0);
%     TT(:,j+1) = phi;
% end
% t = 0:hh:n*hh;
% plot(t,T);
% grid on


    %KB = CC./delta + theta.*KK;
    %QB = (CC./delta - (1-theta).*KK)*T0 + (1-theta).*PP + theta.*PP;
    %T1 = KB\QB;
% I = 0;
% for i = 1:m
%     if T(i,2)-T(i,1)<0
%         I = I + 1;
%     end
% end
