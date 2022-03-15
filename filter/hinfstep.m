function [x,M,K] = hinfstep(x0,th0,M0,t,u,y,gamma)

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
G = full(genmod('dfddw',t,x0,u,th0,w));
f = full(genmod('fd',t,x0,u,th0,w));
g = full(genmod('g',t,x0,u,th0,e));

L = eye(nx);
H = [C;L];
R = blkdiag(eye(ny), -gamma^2*eye(nx));

M = A*M0*A' - A*M0*H'/(H*M0*H' + R)*H*M0*A' + G*G';

Pinv = inv(M0) - gamma^(-2)*(L'*L);

K = A/Pinv*C'/(eye(ny)+C/Pinv*C');

x = f + K*(y-g);

end