function [KK CC PP] = ASSEMBLE(data,kx,ky,kz,rho,c,h,Q,q,TE)
[m,n] = size(data.node);
[i,j] = size(data.element);
%kx = 0.65582;
%ky = 0.65582;
%kz = 0.65582;
%rho = 10^3;
%c = 4.3*10^3;

%KK = zeros(m);
KK1 = zeros(m);
CC = zeros(m);
PPQ = zeros(m,1);
HH = zeros(m);
PPH = zeros(m,1);
PPq = zeros(m,1);

for w = 1:i;
    G = zeros(8,m);
    node1 = data.element(w,2);
    node2 = data.element(w,3);
    node3 = data.element(w,4);
    node4 = data.element(w,5);
    node5 = data.element(w,6);
    node6 = data.element(w,7);
    node7 = data.element(w,8);
    node8 = data.element(w,9);
    x1 = data.node(node1,2);
    y1 = data.node(node1,3);
    z1 = data.node(node1,4);
    x2 = data.node(node2,2);
    y2 = data.node(node2,3);
    z2 = data.node(node2,4);
    x3 = data.node(node3,2);
    y3 = data.node(node3,3);
    z3 = data.node(node3,4);
    x4 = data.node(node4,2);
    y4 = data.node(node4,3);
    z4 = data.node(node4,4);
    x5 = data.node(node5,2);
    y5 = data.node(node5,3);
    z5 = data.node(node5,4);
    x6 = data.node(node6,2);
    y6 = data.node(node6,3);
    z6 = data.node(node6,4);
    x7 = data.node(node7,2);
    y7 = data.node(node7,3);
    z7 = data.node(node7,4);
    x8 = data.node(node8,2);
    y8 = data.node(node8,3);
    z8 = data.node(node8,4);
    
    G(1,node1) = 1; %G = [ 0.......1........0;
    G(2,node2) = 1; %      0
    G(3,node3) = 1; %      0
    G(4,node4) = 1; %      . 
    G(5,node5) = 1; %      .
    G(6,node6) = 1; %      . 
    G(7,node7) = 1; %      . 
    G(8,node8) = 1; %      ...................]
    
    [K,C,~,PQ,~,~] = ELEMENT_2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,kx,ky,kz,rho,c,h,Q,q,TE);
    KK1 = KK1 + G'*K*G;
    CC = CC + G'*C*G;
    PPQ = PPQ + G'*PQ;  % Every Element has Q if Q!= 0
    
    if max(abs([x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8])) == 1;  % Choose the Element of Box;
        %if z1==1||z2==1||z3==1||z4==1||z5==1||z6==1||z7==1||z8==1
            [~,~,H,~,~,PH] = ELEMENT_2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,kx,ky,kz,rho,c,h,Q,q,TE);
            HH = HH + G'*H*G;
            PPH = PPH + G'*PH;
        %end
        %if z1==-1||z2==-1||z3==-1||z4==-1||z5==-1||z6==-1||z7==-1||z8==-1
            [~,~,~,~,Pq,~] = ELEMENT_2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,x5,y5,z5,x6,y6,z6,x7,y7,z7,x8,y8,z8,kx,ky,kz,rho,c,h,Q,q,TE);
            PPq = PPq + G'*Pq;
        %end
    end
end
PP = PPQ + PPq + PPH;
KK = KK1 + HH;
end
