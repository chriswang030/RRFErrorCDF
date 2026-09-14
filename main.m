%% parameters
% matrix parameters
N = 100; % size of matrix
k = 7;   % target rank
p = 7;   % oversampling
assert(N > k+p);

ss_a = 1./(1:N).^2; % full-rank singular values (SVs)
ss_b = [1./(1:2*k+p-3) zeros(1,N-2*k-p+3)]; % rank deficient SVs

% simulation parameters
x_grid = 101;        % resolution of CDF discretization
mc_iter = 1000;      % Monte Carlco iterations for CDF
hist_iter = 10000;   % iterations for histogram
method = "standard"; % method: "standard", "rsi", "rbki"
q = 0;               % iterations for RSI/RBKI 
xx = linspace(0,pi/2,x_grid); % x-axis discretization

%% presets for RBKI; uncomment for RBKI
% method = "rbki";
% q = 2;
% alpha = ss_a(k+1);
% gamma = 0.01;
% cheby = [16 0 -20 0 5 0];
% coeffs = (1+gamma)*alpha/polyval(cheby,1+gamma) * cheby .* alpha.^(-2*q-1:0);

%% generate histogram
fig_a = make_hist(ss_a,k,p,hist_iter,method,q,"#5AB1BB",true);
fig_b = make_hist(ss_b,k,p,hist_iter,method,q,"#5AB1BB",true);

%% generate CDF
figure(fig_a)
hold on
tic
yy = cdf(xx,N,k,p,ss_a,mc_iter,method,q);
line(sin(xx),yy,LineWidth=4,Color="#FF453A");
set(gca,'fontsize',24);
toc
 
figure(fig_b)
hold on
tic
yy2 = cdf(xx,N,k,p,ss_b,mc_iter,method,q);
line(sin(xx),yy2,LineWidth=4,Color="#FF453A");
set(gca,'fontsize',24);
toc