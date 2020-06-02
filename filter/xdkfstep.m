function [x,Px,K,Ps,Psx,Kr] = xdkfstep(x0,th0,P0x,P0s,P0sx,t,u,y,R,Q,Rss,Rys,gamma)

% RES-KF algorithm

if numel(th0) > 1
    error('RES-KF cannot be used for systems with more than 1 parameter.')
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

% Additional noise 
V = Rss;
S = Rys;

% Optimal gain
th_init = th0-1;
ode_th0 = th_init;
ode_A0  = full(genmod('dfddx',t,x0,u,ode_th0,w));
ode_C0 = full(genmod('dgdx',t,x0,u,ode_th0,e));
ode_K0 = ode_A0*P0x*ode_C0'/(ode_C0*P0x*ode_C0'+R);
ode_K0 = reshape(ode_K0,nx*ny,1);
dK_odefun = @(th,K) K_xdkf_ode(th,K,x0,P0x,P0s(:,:,1),P0sx(:,:,1),t,u,y,R,gamma); % define ode
% dK_odefun = @(th,K) K_xdkf_ode_2(th,K,P0x,P0s(:,:,1),P0sx(:,:,1),A,Ath,C,R,gamma); % define ode
[th_sol,K_sol] = ode45(dK_odefun,[th_init,th0],ode_K0);  % solve ode
K = reshape(K_sol(end,:)',nx,ny);
Kr = reshape(dK_odefun(th_sol(end),K),nx,ny);

% Update state and error covariance
x = full(genmod('fd',t,x0,u,th0,w)) + K*(y-C*x0);
Px = (A-K*C)*P0x*(A-K*C)' + Q + K*R*K';
% Update error sensitivity covariance
Ps = P0s;
Psx = P0sx;
P0sp  = P0s(:,:,1);
P0sxp = P0sx(:,:,1);
P0xsp = P0sxp';
Psx(:,:,1) = (Ath-Kr*C)*P0x*(A-K*C)' + (A-K*C)*P0sxp*(A-K*C)' ...
    + Kr*R*K' + K*S*K';
Ps(:,:,1) = (Ath-Kr*C)*P0x*(Ath-Kr*C)' + (Ath-Kr*C)*P0xsp*(A-K*C)' ...
            + (A-K*C)*P0sxp*(Ath-Kr*C)' + (A-K*C)*P0sp*(A-K*C)' ...
            + Kr*R*Kr' + Kr*S*K' + K*S'*Kr' + K*V*K';
        
end