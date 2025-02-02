clear all;
close all;
% Parameters
n = 200; % Number of grid points in x
m = 200; % Number of grid points in y
x = linspace(0, 2*pi, n);
y = linspace(0, 1, m);
dx = x(2) - x(1);
dy = y(2) - y(1);
[X,Y]=meshgrid(x,y);
w=1.7; %relaxation factor

% Initialize temperature field and apply boundary conditions
T = cos(2*X); % initial condition for all points
T(:, 1) = 1; % x=0
T(:, n) = 1; % x=2*pi
T_old=T;

% Convergence criteria
tol = 1e-6;
crit = 1;

% %if matrix A
% a_extracted = diag(A, -1);  % sub-diagonal
% b_extracted = diag(A,  0);  % main diagonal
% c_extracted = diag(A,  1);  % super-diagonal





% Construct the tridiagonal system for row j
        a = -1/dx^2 * ones(1, n-3); % Sub-diagonal
        b = (2/dx^2 + 2/dy^2) * ones(1, n-2); % Main diagonal
        c = -1/dx^2 * ones(1, n-3); % Super-diagonal
        A = spdiags([a b c], [-1 0 1], n-2, n-2); % Sparse tridiagonal matrix

d=ones(1,n-2); % The solutions vector (b)
counter=0;

while crit > tol
    counter=counter+1;
    for j = 2:m-1 % Loop over internal rows
        
        d = (T(j-1 ,2:n-1 ) + T(j+1,2:n-1 )) / dy^2; % RHS
        
        % Adjust RHS for boundary conditions
        d(1) = d(1) + T(j, 1) / dx^2; % Add left boundary condition
        d(end) = d(end) + T(j, n) / dx^2; % Add right boundary condition
        
        % Solve tridiagonal system using the provided Thomas algorithm
        T_new = thomas_algorithm_vectorized(a,b,c,d);
        T(j, 2:n-1)=(1-w)*T(j, 2:n-1)+w*T_new.'; %SOR method
    end

    % Check for convergence
    crit = max(max(abs(T - T_old)));
    T_old=T;
    %disp(crit);
end

%disp(counter);
figure(1);
surf(X,Y,T);
axis([0 2*pi 0 1]);
colorbar;
shading interp;
xlabel('x');
ylabel('y');
title('Temperature Distribution');

function x = thomas_algorithm_vectorized(a, b, c, d)
    n = length(b); % Number of equations
    
    % Forward sweep (elimination step)
    c(1) = c(1) / b(1);
    d(1) = d(1) / b(1);
    
    for i = 2:n-1
        factor = 1 / (b(i) - a(i-1) * c(i-1)); % Compute the scaling factor
        c(i) = c(i) * factor;
        d(i) = (d(i) - a(i-1) * d(i-1)) * factor;
    end
    
    % Last row special case (no super-diagonal element)
    d(n) = (d(n) - a(n-1) * d(n-1)) / (b(n) - a(n-1) * c(n-1));

    % Backward substitution
    x = zeros(n, 1);
    x(n) = d(n);
    
    % Vectorized backward substitution
    for i = n-1:-1:1
        x(i) = d(i) - c(i) * x(i+1);
    end
end


