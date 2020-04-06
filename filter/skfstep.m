function [x,Px,Kx,P,K] = skfstep(x0,th0,P0,t,u,y,R,Q)

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

z0 = [x0;th0];
Qz = blkdiag(Q, zeros(nth));

% Data step

Cx = full(genmod('dgdx',t,x0,u,th0,e));

C = [Cx,zeros(ny,nth)];

K = (P0*C')/(C*P0*C' + R);

z = z0 + K*(y-C*z0);
P = P0 - K*C*P0;

z(nx+1:end,1) = th0;
P(nx+1:end,nx+1:end) = P0(nx+1:end,nx+1:end);

% Time step
x = z(1:nx,1);
Ax = full(genmod('dfddx',t,x,u,th0,w));
Ath = full(genmod('dfddth',t,x,u,th0,w));
fx = full(genmod('fd',t,x,u,th0,w));

A = [Ax, Ath; 
    zeros(nth,nx), eye(nth)];
z = [fx; th0];
P = A*P*A' + Qz;

x = z(1:nx,1);
Px = P(1:nx,1:nx);
Kx = K(1:nx,1:ny);


end