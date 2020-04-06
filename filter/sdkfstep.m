function [x,P,K,S] = sdkfstep(x0,th0,P0,S0,t,u,y,R,Q,W)

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

% Algorithm - single step

A = full(genmod('dfddx',t,x0,u,th0,w));
Ath = full(genmod('dfddth',t,x0,u,th0,w));
C = full(genmod('dgdx',t,x0,u,th0,e));

Xi = A*S0 + Ath;
Ga = C*S0;

K = (A*P0*C' + Xi*W*Ga')/(C*P0*C' + R + Ga*W*Ga');

S = Xi-K*Ga;
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
P = (A-K*C)*P0*(A-K*C)' + Q + K*R*K';

end