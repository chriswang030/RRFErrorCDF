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
%
%   Output:
%   * yy : values of CDF at points specified in xx

if (~exist('method', 'var'))
    method = "standard";
end
if (~exist('q','var'))
    q = 0;
end
if (~exist('coeffs','var'))
    coeffs = [1 zeros(1,2*q+1)];
end

assert(all(ss(1:k) >= 0)); % leading singular values must be non-singular
if method == "standard"
    ss1 = ss(1:k);
    ss2 = ss(k+1:n);
elseif method == "rsi"
    ss1 = ss(1:k).^(2*q+1);
    ss2 = ss(k+1:n).^(2*q+1);
elseif method == "rbki"
    assert(length(coeffs) == 2*q+2); % degree <= 2q+1
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
        R = eye(k) - N*N';
        X = pinv(R'*Q1'*pinv(A)*Q1*R)./ss1./ss1';
        X_eig = eig((X+X')/2);
    end

    y = arrayfun(@(x) J(x,n,k,p,X_eig), xx);
    if any(isnan(y))
        disp();
    end
    yy = yy + y/iter;
end
end

% helper function for computing integrand
function val = J(x,n,k,p,X_eig)
assert(mod(p,2) == 1) % require odd p for terminating series
if x == 0 % avoid Inf * 0 situations for cot(0)
    val = 0;
else
    S = 1./(1+cot(x)^2.*X_eig);
    val = prod(S)^((n-k-p)/2) * mhg([k*(p-1)/2,(p-1)/2],2,(n-k-p)/2,[],1-S);
end
end