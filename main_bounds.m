%% parameters
N = 100;          % size of matrix
kmin = 1;         % min target rank
kmax = 12;        % max target rank
p = 1;            % oversampling
q = 0;            % subspace iterations
delta = 0.05;     % failure probability tolerance
iter = 20000;     % Monte Carlo iterations

% singular values
ss_slow = 1./(1:N).^2;   % slow: order 2 algebraic decay 
ss_fast = max(0.5.^(1:N),1e-10); % fast: exp decay with bound away from 0

%% plots
% plots for Fig. 2-3 top row and Fig. 4
make_estimate_plots(ss_slow,kmin,kmax,p,delta,iter,q);
make_estimate_plots(ss_fast,kmin,kmax,p,delta,iter,q);

% plots for Fig. 2-3 bottom row
figure
xx = linspace(0,pi/2,101);
intercepts = zeros(1,kmax-kmin+1);
for k = kmin:kmax
    yy = cdf(xx,N,k,p,ss_slow,500);
    [~,idx,~] = unique(yy,'first');
    intercepts(k-kmin+1) = sin(interp1(yy(idx),xx(idx),1-delta));
    line(sin(xx),yy,LineWidth=2,Color='k');
    hold on
end
set(gca,'fontsize',24);
set(gca,'TickLength',[0 0]);
set(gca,'TickLabelInterpreter', 'latex')
ylim([0 1.1]);
xlim([0 1]);
plot([0 1], (1-delta)*ones(1,2), LineWidth=3, Color='#FF453A', ...
    LineStyle="--");
plot(intercepts, (1-delta)*ones(1,kmax-kmin+1), LineStyle='none', ...
    Marker='o', Color='k', MarkerFaceColor='k', MarkerSize=10);
grid on
grid minor