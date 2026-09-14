function yy = cdf(xx,n,k,p,ss,iter,method,q,coeffs)
%% CDF
%   Compute the CDF of the largest principal angle between the k-dominant
%   left singular subspace of a matrix A against the (k+p)-dimensional
%   approximation from the randomized range finder algorithm.
%
%   Input:
%   * xx   : discretization points on the interval [0,pi/2]
%   * n    : size of matrix (square, without loss of generality)
%   * k    : target rank >= 0
%   * p    : oversampling >= 0; require odd for exactness
%   * ss   : array of singular values of A
%   * iter : number of Monte Carlo iterations
%   * (optional) method : "standard", "rsi" (randomized subspace
%                 iteration), or "rbki" (randomized block Krylov iteration)
%   * (optional) q      : number of iterations for RSI or RBKI
%   * (optional) coeffs : coefficients of polynomial to test for RBKI;
%                 require odd-degree only, degree <= 2q+1, ordered highest 
%                 degree to lowest
%   * (optional) M      : for even p, set truncation point for the
%                 non-terminating hypergeometric series
%
%   Output:
%   * yy : values of CDF at points specified in xx

if ~exist('method', 'var')
    method = "standard";
end
if ~exist('q','var')
    q = 0;
end
if ~exist('coeffs','var')
    coeffs = [1 zeros(1,2*q+1)];
end

assert(length(ss) == n);
assert(all(ss(1:k) >= 0)); % leading singular values must be non-singular

% address an edge case by embedding in larger space
if n < 2*k+p
    ss = [ss zeros(1,2*k+p-n)];
    n = 2*k+p;
end

if method == "standard" % usual randomized range finder
    ss1 = ss(1:k);
    ss2 = ss(k+1:n);
elseif method == "rsi"  % randomized subspace iteration
    ss1 = ss(1:k).^(2*q+1);
    ss2 = ss(k+1:n).^(2*q+1);
elseif method == "rbki" % randomized block Krylov iteration
    assert(length(coeffs) == 2*q+2); % degree must be <= 2q+1
    assert(all(coeffs(2:2:length(coeffs)) == 0)); % no even-degree terms
    ss1 = polyval(coeffs,ss(1:k));
    ss2 = polyval(coeffs,ss(k+1:end));
end

yy = zeros(size(xx));

% Monte Carlo iteration
for i = 1:iter
    [H1,~] = qr(randn(n-k,k+p),'econ');
    [Q1,~] = qr(randn(k+p,k),'econ');

    % full-rank case
    if sum(ss2 > eps) >= k+p
        [C,r] = chol(Q1'/(H1'.*ss2.^2*H1)*Q1);
        if r == 0
            Xh = diag(ss1.^(-1))/C;
            X_eig = eig(Xh*Xh');
        else
            Xh = diag(ss1.^(-2))/(Q1'/(H1'.*ss2.^2*H1)*Q1);
            X_eig = eig((Xh+Xh')/2);
        end
        
    % rank-deficient case
    else
        A1 = H1'.*ss2;
        A = A1*A1';
        [N,~] = qr(Q1'*null(A),'econ');
        P = null(N');
        R = P*P';
        X = pinv(R'*Q1'*pinv(A)*Q1*R)./ss1./ss1';
        X_eig = eig((X+X')/2);
    end

    if mod(p,2) == 1
        y = arrayfun(@(x) integrand_odd(x,n,k,p,X_eig), xx);
    else
        y = even_constant(n,k,p) ...
            * arrayfun(@(x) integrand_even(x,n,k,p,X_eig), xx);
    end
    yy = yy + y/iter;
end
end

% helper function for computing integrand
function z = integrand_odd(x,n,k,p,X_eig)
if x == 0 % avoid Inf * 0 situations for cot(0)
    z = 0;
    return
end

% set matrix argument
S = 1./(1+cot(x)^2 .* X_eig);

% use Eq. (3.18) for odd p with terminating series
z = prod(S)^((n-k-p)/2) * mhg([k*(p-1)/2,(p-1)/2],2,(n-k-p)/2,[],1-S);
return
end

% use Eq. (3.1) for even p
function z = integrand_even(x,n,k,p,X_eig)
if x == 0 % avoid Inf * 0 situations for cot(0)
    z = 0;
    return
end

% set matrix argument
S = 1./(1+cot(x)^2 .* X_eig);
z = prod(S)^((n-k-p)/2) * mhg(100,2,[(-p+1)/2,(n-k-p)/2],(n-p+1)/2,S);
end

function K = even_constant(n,k,p)
K = 1;
for j = 1:k
    K = K * gamma((n-j+1)/2) * gamma((k-j+2)/2) ...
        / gamma((n-p-j+2)/2) / gamma((k+p-j+1)/2);
end
end
