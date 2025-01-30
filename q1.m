close all;
clear all;
% Parameters
n = 200; % Number of grid points in x
m = 200; % Number of grid points in y
x = linspace(0, 2*pi, n);
y = linspace(0, 1, m);
dx = x(2) - x(1);
dy = y(2) - y(1);
[X,Y]=meshgrid(x,y);

% Initialize temperature field and apply boundary conditions
T = cos(2*X); % initial condition for all points
T(2:m-1, 1) = 1; % x=0
T(2:m-1, n) = 1; % x=2*pi


% Convergence criteria
tol = 1e-6;
crit = 1;
T_old = 1;


% Construct the tridiagonal system for row j
        a = -1/dx^2 * ones(1, n-3); % Sub-diagonal
        b = (2/dx^2 + 2/dy^2) * ones(1, n-2); % Main diagonal
        c = -1/dx^2 * ones(1, n-3); % Super-diagonal
        %A = spdiags([a b c], [-1 0 1], n-2, n-2); % Sparse tridiagonal matrix
        A = diag(b) + diag(a, 1) + diag(c, -1);
d=ones(n-2); % The solutions vector (b)
tic;
while crit > tol
    
    for j = 2:m-1 % Loop over internal rows
        
        d = (T(j-1 ,2:n-1 ) + T(j+1,2:n-1 )) / dy^2; % RHS
        
        % Adjust RHS for boundary conditions
        d(1) = d(1) + T(j, 1) / dx^2; % Add left boundary condition
        d(end) = d(end) + T(j, n) / dx^2; % Add right boundary condition
        
        % Solve tridiagonal system using the provided Thomas algorithm
        T(j, 2:n-1) = tomas(A,d);
        

    end
    
    % Check for convergence
    crit = max(max(abs(T - T_old)));
    T_old=T;
    %disp(crit);
end
toc;
 figure(1);
 contourf(X,Y,T);
 axis([0 2*pi 0 1]);
 colorbar;
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

function [m]=tomas(A,d)
%A is sparse
%solve linear equation system, A is a matrix and d is vector of results and
%m is the answer vector
I=length(d); %get the length of the vector of results
m=zeros(1,I);%preallocation
%fix the matrix for 3- diagonal matrix
d(I)=d(I)-d(I-2)*(A(I,I-3)/(A(I-2,I-3)));
A(I,:)=A(I,:)-A(I-2,:).*(A(I,I-3)/(A(I-2,I-3)));
d(I)=d(I)-d(I-1)*(A(I,I-2)/(A(I-1,I-2)));
A(I,:)=A(I,:)-A(I-1,:).*(A(I,I-2)/(A(I-1,I-2)));
A(I,:)=A(I,:)-A(I-2,:).*(A(I,I-3)/(A(I-2,I-3)));
A(I,:)=A(I,:)-A(I-1,:).*(A(I,I-2)/(A(I-1,I-2)));

E=zeros(1,I);%preallocation
F=E;%preallocation
E(1)=A(1,2)/A(1,1); F(1)=d(1)/A(1,1);%intial condition
for L=2:I %run the algorithm
    if(L~=I)
       E(L)=(A(L,L+1))/(A(L,L)-E(L-1)*A(L,L-1)); 
    end
    F(L)=(d(L)-F(L-1)*A(L,L-1))/(A(L,L)-E(L-1)*A(L,L-1));
end
%applied the solutions:
m(I)=F(I);
for l=I-1:-1:1
    m(l)=F(l)-E(l)*m(l+1);
end
end