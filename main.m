%% parameters
% matrix parameters
N = 100; % size of matrix
k = 2;   % target rank
p = 1;   % oversampling
ss1 = 1./(1:N).^2; % full-rank singular values (SVs)
ss2 = [1./(1:2*k+p-3) zeros(1,N-2*k-p+3)]; % rank-deficient SVs

% simulation parameters
x_grid = 101;        % resolution of CDF discretization
mc_iter = 1000;      % Monte Carlo iterations for CDF
hist_iter = 10000;   % iterations for histogram
method = "standard"; % method: "standard", "rsi", "rbki"
q = 0;               % iterations for RSI/RBKI 
xx = linspace(0,pi/2,x_grid); % x-axis discretization

%% presets for RBKI; uncomment for RBKI
% method = "rbki";
% q = 2;
% alpha = ss(k);
% gamma = 0.01;
% cheby = [16 0 -20 0 5 0];
% coeffs = (1+gamma)*alpha/polyval(cheby,1+gamma) * cheby .* alpha.^(2*q+1:-1:0);

%% generate histogram
fig1 = make_hist(ss1,k,p,hist_iter,method,q,"#5AB1BB");
fig2 = make_hist(ss2,k,p,hist_iter,method,q,"#5AB1BB");

%% generate CDF
figure(fig1)
hold on
tic
yy = cdf(xx,N,k,p,ss1,mc_iter,method,q);
line(xx,yy,LineWidth=4,Color="#FF453A");
toc

figure(fig2)
hold on
tic
yy = cdf(xx,N,k,p,ss2,mc_iter,method,q);
line(xx,yy,LineWidth=4,Color="#FF453A");
toc