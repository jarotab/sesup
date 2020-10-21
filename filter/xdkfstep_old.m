function [x,Px,K,Ps,Psx,Kth] = xdkfstep(x0,th0,P0x,P0s,P0sx,t,u,y,R,Q,gamma)

% XDKF algorithm

if numel(th0) > 1
    error('XDKF cannot be used for systems with more than 1 parameter.')
end

nth = 1;
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

A  = full(genmod('dfddx',t,x0,u,th0,w));
Ath_col = full(genmod('ddfddxdth',t,x0,u,th0,w));
Ath = reshape(Ath_col(:,1),nx,nx);
C = full(genmod('dgdx',t,x0,u,th0,e));

% Optimal gain
Kth = zeros(nx,ny);
dAdA_ij_0 = zeros(nx,nx);
for i = 1:nx
    for j = 1:nx
        dAdA_ij = dAdA_ij_0;
        dAdA_ij(i,j) = double(abs(A(i,j))>1e-6);
        dAdth_ij = Ath(i,j);
        Kth = Kth + dKdth_ij_fun(P0x,P0s(:,:,1),dAdA_ij,dAdth_ij,C,R,gamma);
    end
end
K = K_fun(th0,Kth,x0,P0x,P0s(:,:,1),P0sx(:,:,1),t,u,y,A,C,R,gamma);
K = reshape(K,nx,ny);

% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';
% Update error sensitivity covariance
Ps = P0s;
Psx = P0sx;
P0sp  = P0s(:,:,1);
P0sxp = P0sx(:,:,1);
P0xsp = P0sxp';
Psx(:,:,1) = (-Kth*C)*P0x*(A-K*C)' + (A-K*C)*P0sxp*(A-K*C)' ...
    + Kth*R*K';
Ps(:,:,1) = (-Kth*C)*P0x*(-Kth*C)' + (-Kth*C)*P0xsp*(A-K*C)' ...
            + (A-K*C)*P0sxp*(-Kth*C)' + (A-K*C)*P0sp*(A-K*C)' ...
            + Kth*R*Kth';
        
end

function dKdth_ij = dKdth_ij_fun(Px,Ps,dAdA_ij,dAdth_ij,C,R,gamma)

% XDKF algorithm
nx = size(Px,1);
ny = size(C,1);

alpha = 1-sum(gamma);

% Optimal gain
Wij = -alpha*(dAdA_ij*Px*C') - gamma*(dAdA_ij*Ps*C');
U = alpha*(C*Px*C'+R) + gamma*(C*Ps*C');

dKdth_ij = Wij/U*dAdth_ij; 

end

function K = K_fun(th,dK,x,Px,Ps,Psx,t,u,y,A,C,R,gamma)

% XDKF algorithm
nx = size(x,1);
ny = size(y,1);

dK = reshape(dK,nx,ny);

alpha = 1-sum(gamma);

% Optimal gain
Pxs = Psx';
W = alpha*(A*Px*C') + gamma*(A*Ps*C');
U = alpha*(C*Px*C'+R) + gamma*(C*Ps*C');
VrT = gamma*(C*Pxs*C');

K = (-dK*VrT + W)/U; 

K = reshape(K,nx*ny,1);

end