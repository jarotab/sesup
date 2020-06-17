function dK = K_xdkf_ode(th,K,x,Px,Ps,Psx,t,u,y,R,gamma)

% XDKF algorithm

nx = size(x,1);
ny = size(y,1);
nu = size(u,1);
if nu<1
    u = 0;
end

w = zeros(nx,1);
e = zeros(ny,1);

x = reshape(x,nx,1);
K = reshape(K,nx,ny);

% Additional noise 
V = 10*R;
S = R;

A  = full(genmod('dfddx',t,x,u,th,w));
Ath_col = full(genmod('ddfddxdth',t,x,u,th,w));
C = full(genmod('dgdx',t,x,u,th,e));

alpha = 1-sum(gamma);

% Optimal gain
Athp = reshape(Ath_col(:,1),nx,nx);
Pxs = Psx';
W = alpha*(A*Px*C') + gamma*(Athp*Pxs*C' + A*Ps*C');
U = alpha*(C*Px*C'+R) + gamma*(C*Ps*C'+V);
VrT = gamma*(C*Pxs*C'+S);

dK = - K*U/VrT + W/VrT; 

dK = reshape(dK,nx*ny,1);

end