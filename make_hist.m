function fig = make_hist(ss,k,p,iter,method,q,color,sin_scale)
%% MAKE_HIST
%   Plot histogram of empirical subspace approximation error from samples
%   of RRF applied to a matrix A, in CDF format.
%
%   Input:
%   * ss   : array of singular values of A
%   * k    : target rank >= 0
%   * p    : oversampling >= 0; require odd for exactness
%   * iter : number of RRF samples
%   * (optional) method : "standard", "rsi" (randomized subspace
%                 iteration), or "rbki" (randomized block Krylov iteration)
%   * (optional) q      : number of iterations for RSI or RBKI
%   * (optional) color  : color of histogram bars
%   * (optional) sin_scale : whether to plot CDF of angular error theta
%                 itself or scale to sin(theta)
%
%   Output:
%   * fig : MATLAB figure with histogram

if (~exist('method', 'var'))
    method = "standard";
end
if (~exist('q','var') || method == "standard")
    q = 0;
end
if (~exist('color','var'))
    color = "white";
end
if (~exist('sin_scale','var'))
    sin_scale = false;
end

fig = figure;
n = length(ss);
data = zeros(1,iter);
for i = 1:iter
    X = randn(n,k+p);
    K = ss'.*X;
    if method == "standard" || method == "rsi"
        for j = 1:q
            K = ss.^2'.*K; % equivalent to diag(ss)' * diag(ss) * X
        end
        [Y,~] = qr(K,'econ');
    elseif method == "rbki"
        K_hat = [K zeros(n,q*(k+p))];
        for j = 1:q
            K = ss.^2'.*K;
            K_hat(:,j*(k+p)+1:(j+1)*(k+p)) = K;
        end
        [Q,~] = qr(K_hat,'econ');
        [U,~,~] = svds(Q'.*ss,k+p); % equivalent to Q'*diag(ss)
        Y = Q*U;
    end

    Y1 = Y(1:k,1:k+p);
    sv = svds(Y1,1,'smallest');
    data(i) = acos(min(sv,1));
end

if sin_scale
    binwidth = 1/30;
    binlim = [0 1];
    xts = 0.2*(0:5);
    data = sin(data);
else
    binwidth = pi/60;
    binlim = [0 pi/2];
    xts = pi/16*(0:8);
end

histogram(data, ...
    Normalization='cdf', ...
    BinWidth=binwidth, ...
    BinLimits=binlim, ...
    LineWidth=1, ...
    FaceColor=color, ...
    FaceAlpha=0.8);

set(gca,'fontsize',24);
set(gca,'TickLength',[0 0]);
set(gca,'TickLabelInterpreter', 'latex')
ylim([0 1.1]);
xlim(binlim);
xticks(xts);
if ~sin_scale
    xticklabels({'0','$\frac\pi{16}$','$\frac\pi{8}$','$\frac{3\pi}{16}$', ...
        '$\frac{\pi}{4}$','$\frac{5\pi}{16}$','$\frac{3\pi}{8}$', ...
        '$\frac{7\pi}{16}$','$\frac{\pi}{2}$'})
end
end