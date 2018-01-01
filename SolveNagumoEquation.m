tic;
h1 = 0.01; % for x,y
h2 = 0.25; % for t
 
m = round(1 / h1);
n = (m - 1)^2;
r = h2 / (h1^2);
e = 0.1;
a = 0.5;    % normally in typical Nagumo equation a should be in [0, 1]
Tmax = 2;
 
I = eye(n);
 
A = -4 * eye(n);
 
for i = 1 : (n-1)
    if (mod(i, m - 1) == 0)
        A(i, i + 1) = 0;
        A(i + 1, i) = 0;
    else
        A(i, i + 1) = 1;
        A(i + 1, i) = 1;
    end    
end
 
for i = 1 : (n - m + 1)
    A(i, m - 1 + i) = 1;
    A(m - 1 + i, i) = 1;
end
    
A = I - (e * r * A);
 
b1 = zeros(n, 1);
k = 1;
for i = 0 : m : m
    if (i == m)
        k = n - m + 2;
    end
    for j = 1 : (m-1)
        x = i * h1;
        y = j * h1;
        o = u(x, y);
        b1(k) = o;
        k = k + 1;
    end
end
 
b1 = e * r * b1;
 
b2 = zeros(n, 1);
step = m - 2;
 
k = 1;
for i = 1 : (m-1)
    for j = 0: m: m
        x = i * h1;
        y = j * h1;
        o = u(x, y);
        b2(k) = o;
        if (j == 0)
            k = k + step;
        else
            k = k + 1;
        end
    end
end
    
b2 = e * r * b2;
 
b3 = zeros(n, 1);
k = 1;
for i = 1 : (m - 1)
    for j = 1 : (m - 1)
        x = i * h1;
        y = j * h1;
        o = u(x, y);
        b3(k) = o * (1 - o) * (o - a);
        k = k + 1;
    end
end
    
b3 = h2 * b3;
b = b1 + b2 + b3;
 
K1 = diag(diag(A));
K2 = inv(K1);
A_prec = inv(K1) * A * K1;
b_prec = K2 * b;
 
maxIter = 500;
tol = 0.01;
tNew = BiCGSTAB_unprec( A_prec, b_prec, maxIter, tol );
 
m_t = Tmax/h2;
solution = zeros(n, m_t - 2);
solution(:,1) = tNew;
for k = 2 : (m_t-1);
    b = b1 + b2 + tNew;
    b_prec = K2 * b;
    tNew = BiCGSTAB_unprec( A_prec, b_prec, maxIter, tol );
    solution(:, k) = tNew;
end
toc;
