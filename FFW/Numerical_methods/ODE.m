%% SOLVE y'(x) = cos( (x^2 + x^3)^0.25 )
%        y(0)  = 0
%
%  Objective: Solve for y(x)
%             Demonstrate convergence by varying step-size h

clc
clear
clf

h = 0.1; % Step-size 0.1 0.5 0.05

% Initial condition y(x=0) = 0;
y(1) = 0;

% Setup Plot
hold on
xlabel('x')
ylabel('y(x)')
axis([0 12 -3.5 1])

% Iterate
N = 12/h; % total number of iterations to perform
n = 0;    % current number of iterations performed
while n < N
    
    n = n + 1;
    
    x(n) = (n-1)*h;
    f(n) = cos( (x(n)^2 + x(n)^3)^0.25 );
    y(n+1) = y(n) + h*f(n);
    
    plot(x(n),y(n),'r.',x,y(1:n),'r')
    pause(0.01)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SOLVE y'(x) = -4*y(x)  
%        y(0)  = 1
%  Exact solution: y(x) = exp(-4x)
%
%  Objective: Demonstrate dependence of step-size h on numerical stability

clc
clear
clf

h = 0.1; % Step-size 0.1 0.5 1

% Initial condition y(x=0) = 0;
y(1) = 1;

% Setup Plot
hold on
xlabel('x')
ylabel('y(x)')

% Iterate
N = 8/h; % total number of iterations to perform
n = 0;   % current number of iterations performed
while n < N
    
    n = n + 1;
    
    x(n) = (n-1)*h;
    f(n) = -4*y(n);
    y(n+1) = y(n) + h*f(n);
    
    plot(x(n),y(n),'r.',x,y(1:n),'r',x,exp(-4*x),'b')
    pause(0.01)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%