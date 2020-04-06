function [x,P,K] = kfstep(x0,th0,P0,t,u,y,R,Q)

nth = size(th0,1);
nx = size(x0,1);
ny = size(y,1);
nu = size(u,1);
if nu<1
    u = 0;
end

w = zeros(nx,1);
e = zeros(ny,1);

x0 = reshape(x0,nx,1);
th0 = reshape(th0,nth,1);
y = reshape(y,ny,1);

A = full(genmod('dfddx',t,x0,u,th0,w));
C = full(genmod('dgdx',t,x0,u,th0,e));

K = (A*P0*C')/(C*P0*C' + R);

x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
P = A*P0*A' + Q - K*C*P0*A';

end