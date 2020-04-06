%% Mandatory definition

sys = struct;

% Propagation equation
sys.fd = @(t,x,u,th,w) [0.9802 0.0196+0.099*th(1); 0  0.9802]*x + eye(2)*w;

% Measurement equation
sys.g = @(t,x,u,th,e) [1 -1]*x + e;

% Dimensions
sys.nx = 2;
sys.ny = 1;
sys.nth = 1;
sys.nu = 0;

% Name
name = 'Sayed (2001)';

%% Optional 

opt = struct;

% Noise covariances
opt.R = 1;
opt.Q = [1.9608 0.0195; 
     0.0195 1.9605];

% Initial conditions
opt.x0 = ones(2,1)*0;
opt.P0 = eye(2);
opt.y0 = 0;

% Parameter prior information
opt.th0 = 0;  % Mean parameter
opt.thBounds = [-1,1];  % Parameter bounds

%% Save system

save('./sayed.mat','sys','opt','name');


