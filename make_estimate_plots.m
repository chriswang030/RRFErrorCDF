function [fig1, fig2] = make_estimate_plots(ss,kmin,kmax,p,delta,iter,q)

ss = ss(:);     % ensure column vector
N = length(ss); % size of matrix
q1 = 2*q+1;     % exponent for subspace iterations
xx = linspace(0,pi/2,501);

% sine bound function
f = @(y) y ./ sqrt(1+y.^2);

%% compute error estimates for various k
plot1 = zeros(3,kmax-kmin+1);
plot2 = zeros(4,kmax-kmin+1);
for k = kmin:kmax
    j = k-kmin+1;

    %% compute relevant constants
    e = exp(1);
    mu = (p-1)/2;
    nu = k*(N-k-p)/2;

    % list of constants, C0-C5 are from our paper
    C = zeros(3);
    C(1) = sqrt(k*(N-k-1) / (2*delta)); % C0
    if mod(p,2) == 1
        C(2) = (local_binom(nu+mu,mu+1) / delta)^(1/(p+1)); % C1
    elseif p > 0
        C(2) = (local_binom(nu+mu,mu+1) / delta)^(1/(p+1)); % C3
    else
        C(2) = k*(N-k) / delta; % C5
    end

    % constant from Eq. (12) in Saibaba (2019)
    C(3) = e*sqrt(k+p)/(p+1) * (2/delta)^(1/(p+1)) ...
        * (sqrt(N-k) + sqrt(k+p) + sqrt(2*log(2/delta)));

    %% compute relevant singular value ratios
    % SV ratio
    rho = (ss(k+1)/ss(k))^q1;
    % Frobenius SV ratio
    xi = sqrt(mean(ss(1:k).^(-2*q1)) * mean(ss(k+1:end).^(2*q1)));
    
    %% compute error estimates
    % Monte Carlo quantities from RSVD
    actual_data = zeros(1,iter);
    ss_approx = zeros(1,N);
    for i = 1:iter
        X = randn(N,k+p);
        K = ss.^q1.*X;     % subspace iterations
        [Y,~] = qr(K,'econ');

        % sin of true largest principal angle
        Y1 = Y(1:k,1:k+p);
        sv = min(1,svds(Y1,1,'smallest')); % adjust for numerical error for acos
        actual_data(i) = sin(acos(sv));
        
        % Monte Carlo for approximate SVs
        ss_approx(1:k+p) = svd(ss'.*Y')'/iter;
    end

    % pad tail for approximate SVs
    ss_approx(k+p+1:end) = ss_approx(k+p);

    % error estimate from CDF inversion + exact/approx. SVs
    zz = cdf(xx,N,k,p,ss_approx.^q1,500);
    zz = real(zz(~isnan(zz)));
    xx_temp = xx(~isnan(zz));
    [~,idx,~] = unique(zz,'first');
    bd_cdf = sin(interp1(zz(idx),xx_temp(idx),1-delta));

    % get (1-delta)th percentile from sample data
    actual = prctile(actual_data,(1-delta)*100);

    plot1(:,j) = [ ...
         actual; ...        % empirical
         f(xi*C(1)); ...    % bound from Thm 4.7
         f(rho*C(3)) ];     % bound from Saibaba (2019) Thm 6
    plot2(:,j) = [ ...
         actual; ...
         f(xi*C(1)); ...    % bound from Thm 4.7
         f(xi*C(2)); ...    % bound from Conj 5.1
         bd_cdf ];          % bound from RSVD approximation
end

%% plot
% colors
black  = 'k';
blue   = '#0072BD';
green  = '#77AC30';
orange = '#D95319';

p1_colors = {black, blue, orange};
p2_colors = {black, blue, blue, green};
p1_styles = {"-", "-", "-"};
p2_styles = {"-", "-", "--", "-"};

fig1 = figure;
for j = 1:3
    semilogy(kmin:kmax,plot1(j,:), ...
        LineWidth=2, ...
        LineStyle=p1_styles{j}, ...
        Marker='o', ...
        MarkerFaceColor=p1_colors{j}, ...
        MarkerSize=8, ...
        Color=p1_colors{j});
    hold on
end

set(gca,'fontsize',24);
set(gca,'TickLength',[0 0]);
set(gca,'TickLabelInterpreter', 'latex')
ylim([min(plot1,[],'all')/1.2 1.1]);
yticks(0:0.1:1)
xlim([0 kmax]);
xticks(0:5:kmax);
grid on

fig2 = figure;
for j = 1:4
    semilogy(kmin:kmax,plot2(j,:), ...
        LineWidth=2, ...
        LineStyle=p2_styles{j}, ...
        Marker='o', ...
        MarkerFaceColor=p2_colors{j}, ...
        MarkerSize=8, ...
        Color=p2_colors{j});
    hold on
end

set(gca,'fontsize',24);
set(gca,'TickLength',[0 0]);
set(gca,'TickLabelInterpreter', 'latex')
ylim([min(plot2,[],'all')/1.2 1.1]);
yticks(0:0.1:1)
xlim([0 kmax]);
xticks(0:2:kmax);
grid on
end

function x = local_binom(a,b)
if b >= 0 && mod(b,1) == 0
    if a >= b && mod(a,1) == 0
        x = nchoosek(a,b);
    elseif b > 0
        x = prod(a-b+1:a ./ 1:b);
    else
        x = 1;
    end
else
    x = gamma(a+1) / gamma(b+1) / gamma(a-b+1);
end
end