clc; clear all; close all;

%% parameters

c=50; %grid size
x=linspace(0,2*pi,c);
y=linspace(0,1,c);
dx=x(2)-x(1);
dy=y(2)-y(1);
dt = 0.0001;

%grid coordinates
[X,Y]=meshgrid(x,y);

%cylinder coordinates
phi = linspace(0,2*pi,100);
R = 0.2;
cylinder = [pi+R*cos(phi);1/2+R*sin(phi)];
eps = cylinder(1,:);
eta = cylinder(2,:);
%boundary condition on the cylinder surface
T_cyl = 2;


% %initial condition on grid
% TT = [ones(c*c,1)];
% for ii=1:c
%     TT = [ones(c,1)*cos(x(ii));TT];
% end

%% initialazing mat A and vec b, including IC

%initialazing the matrix
b=ones(c*c,1); %IC of one in whole domain
q=zeros(c*c,1);
alpha=2*(1/dx^2+1/dy^2)+1/dt;
beta=-1/dx^2;
gamma=-1/dy^2;

% for row i theres:
% alpha on column i
% beta on column i-1,i+1
% gamma on column i+n,i-n
vec_alpha=alpha*ones(c*c-c,1);
vec_beta=beta*ones(c*c-c,1);
vec_gamma=gamma*ones(c*c-c,1);
A1=spdiags([vec_gamma vec_beta vec_alpha vec_beta vec_gamma], [1 c c+1 c+2 2*c+1 ],c*c-2*(c+1) , c*c);
A2=spdiags([ones(c+1,1)], [0], c+1,c*c);
A3=spdiags([ones(c+1,1)], [c*c-c-1], c+1,c*c);
A=[A2;A1;A3];

for ii=1:c 
    if ii>=2
        A(ii*c:ii*c+1,:)=0;
        A(ii*c,ii*c)=1;   
        A(ii*c+1,ii*c+1)=1; 
    end
    % defining the parts of b that are based on the boundry conditions
    % where x=0 and x=2*pi
    b((ii-1)*c+1)=1;
    b(ii*c)=1;
    disp(b);
end
% defining the parts of b that are based on the boundry conditions where
% y=0 and y=1
b(1:c)=cos(2*x(ii))/dt;
b(c*c-c+1:c*c)=cos(2*x(ii))/dt;

%% calculating transient proccess

for n=0:dt:1
    %step 1: calculate intermidiate temp b
    %step 2: interpolation temp to body
    %step 3: calculate het flux q
    %step 4: regularize q to grid 
    %step 5: solve for new temp field
end

