clc
clear
clf

% SOLVE y'(x) = cos( (x^2 + x^3)^0.25 )
%       y(0)  = 0

h = 0.1; % Step-size 

% Initial condition y(x=0) = 0;
y(1) = 0;

% Setup Plot
hold on
xlabel('x')
ylabel('y(x)')
axis([0 12 -3.5 1])

% Iterate
N = 12/h; 

n = 0;
while n < N
    
    n = n + 1;
    
    x(n) = (n-1)*h;
    f(n) = cos( (x(n)^2 + x(n)^3)^0.25 );
    y(n+1) = y(n) + h*f(n);
    
    plot(x(n),y(n),'b.')
    pause(0.01)
end