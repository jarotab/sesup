function [x,P,K,s] = dkfstep(x0,th0,P0,s0,t,u,y,R,Q,W,Ks)

nth = size(th0,1);
nx = size(x0,1);
ny = size(y,1);
nu = size(u,1);
if nu<1
    u = 0;
end

x0 = reshape(x0,nx,1);
th0 = reshape(th0,nth,1);
y = reshape(y,ny,1);

w = zeros(nx,1);
e = zeros(ny,1);

C = full(genmod('dgdx',t,x0,u,th0,e));

% Data step

S = C*P0*C' + R;
SumL = Ks*0;
SumR = Ks*0;
for i = 1:nth
    g = C*s0(:,i);
    SumL = SumL + W(:,:,i)*Ks*g*g';
    SumR = SumR +  W(:,:,i)*s0(:,i)*g';
end
Ksv = reshape(Ks,1,nx*ny);
Ksol = vpasolve(Ks*S + SumL == P0*C' + SumR,Ksv);
Kvec = zeros(1,nx*ny);
for i = 1:numel(Ksv)
Kvec(i) = eval(Ksol.(char(Ksv(i))));
end
K = reshape(Kvec,nx,ny);

s = s0*0;
for i = 1:nth
    s(:,i) = (eye(nx) - K*C)*s0(:,i);
end
x = x0 + K*(y-C*x0);
P = (eye(nx)-K*C)*P0*(eye(nx)-K*C)' + K*R*K';

% Time step
A = full(genmod('dfddx',t,x,u,th0,w));
Ath = full(genmod('dfddth',t,x,u,th0,w));

for i = 1:nth
    s(:,i) = A*s(:,i) + Ath(:,i);
end
x = full(genmod('fd',t,x,u,th0,w));
P = A*P*A' + Q;

end