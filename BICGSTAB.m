function [ x ] = BiCGSTAB_unprec( A, b, maxIter, tol )
tic;
n = size(A, 2);
x = rand(n, 1);
%x = ones(n, 1);
%x = zeros(n, 1);
r = b - A * x;
r_hat = r;
rho = 1;
alpha = 1;
w = 1;
v = zeros(n, 1);
p = zeros(n, 1);
for i = 1 : maxIter
    rho_old = rho;
    rho = dot(r_hat, r);
    beta = (rho / rho_old) * (alpha / w);
    p = r + beta * (p - w * v);
    v = A * p;
    alpha = rho / dot(r_hat, v);
    h = x + alpha * p;
    if norm(h-x, inf) < tol
        x = h;
        break;
    end
    s = r - alpha * v;
    t = A * s;
    w = dot(t, s) / dot(t, t);
    x_new = h + w * s;
    if norm(x_new-x, inf) < tol
        x = x_new;
        break;
    end
    x = x_new;
    r = s - w * t;
end
fprintf('Algorithm converged in %d iterations\n', i);
toc;
end
