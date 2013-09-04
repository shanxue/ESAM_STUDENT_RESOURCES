%Introduction to MATLAB: tips, tricks, how to?s, and resources
%Karna Gowda
%ESAM First Year Fundamentals Workshop 2013

%Contents:

%Intro: What is MATLAB, how does it differ from other programming languages
% Obtaining MATLAB, installing it, creating a script, running a script.
% Talk about command window, workspace, directories, commenting, functions 
% vs. scripts, semicolons, spaces, tabs, version control.

%Simple scripts: Run through a simple Forward Euler scheme

%Plotting: Plotting in MATLAB, saving files and images

%Protips: Measuring script runtime, profiler, using built-in functions,
% vectorizing, indexing

%Resources: Useful links, see accompanying file. Include a bunch of useful
% functions

%% Simple script: Forward Euler
%This is the most basic numerical integration method out there.
%You can derive it from a Taylor expansion

%Diff eq: x'(t) = f(x,t), solution: x(t), IC: x(t_0) = x_0
%Taylor expand solution about t_0: x(t_0 + dt) = x(t_0) + x'(t_0) dt + ...
% (x(t_0 + dt) - x(t_0))/dt -> x'(t_0) = f(x(t_0), t_0) as dt -> 0
%Then: x(t_0 + dt) = x(t_0) + dt*f(x(t_0), t_0)
%So FE step is: x_n+1 = x_n + dt*f(x(t_n), t_n)

%Variables: we want to define a few variables to hold details like the IC,
% the time step size dt, the number of time steps taken (or the max time),
% etc. In MATLAB, you don't have to declare variable types.

dt=0.01;    %this gets automatically typed as as float (decimal) variable
t0=0;       %these next couple also get typed as floats, even though they 
tMax = 10;  % look like integers
x0 = 0.1;

%A variable can contain references to other variables. Note that we should 
% make sure this is an integer.
tSteps = tMax/dt;

%If we want to store all the output of our method, i.e. the entire solution
% trajectory, we can do that with a vector. This is basically an array. You
% can append to a vector on the fly, but MATLAB likes it a bit more if you
% preallocate a vector. We do this by creating a giant vector of zeros,
% with a space for each point resolved on our trajectory. We'll use the
% MATLAB zeros(rows,cols) function to create this.

x = zeros(1,tSteps); %protip: zeros(N) will create an NxN matrix

%We want the first element of the trajectory to be the IC.

x(1) = x0; %notice that indexing starts with 1 in MATLAB. Annoying.

%In addition, we might want a vector of time values. MATLAB's linspace
% function can do that.
t = linspace(t0,tMax,tSteps);

%Finally, we need our f(x,t) function. We can define this inline with the
% function handle formalism. 

f = @(x,t) x - x.^3; %dependent variables go in parens after the @ symbol

%If we want to call this function, we would simply say "f(0,0)" for
% example. Notice that even though we've already defined x and t, this
% function treats x and t as symbolic variables.

%That's pretty much all the setup that we need. Time for the loop. We have
% some choices here, but I'm just going to do a simple "for" loop. We'll 
% create a counter variable "i" in the for loop statement. i will take
% values in the range 2 to tSteps. The contents of the for loop must go
% between the "for" statement and "end".
for i=2:tSteps
    %Now we just implement the Forward Euler step.
    x(i) = x(i-1) + dt*f(x(i-1),t(i-1));
end

%Great! Running this section now will numerically solve the system 
% x'(t) = x - x^3 with ICs t=t0, x(t_0)=x_0. Now what?

%% Plotting: Basic 2D lineplot
%A simple 2D lineplot is the "plot(t,x)" function. This can accept vector
% inputs. By running the section above, you've generated numerical solution
% vectors x and t already, and they are sitting in your "workspace." That
% is to say, we can plot them right now. It's easy to add labels as well.

plot(t,x)
xlabel('t') %Ok, "xlabel" means independent variable label
ylabel('x') %dependent variable label

%% Plotting: Getting fancy with a lineplot
%We can change the appearance of the plot by assigning it a variable name
% so we can change its object properties. The main box containing the plot
% is a "figure" and we'll name that too, so we can save it later.

close all %this command closes all previously open plots
f = figure
p = plot(t,x)
xlabel('t')
ylabel('x')
set(p,'Color','red','LineWidth',2)

%% Plotting: Saving your output
%Now that you've generated a nice plot, you want to save it and put it in a
% paper or something. My favorite format to save in is "eps", which saves
% your plot as a vectorized image. This basically means that the file has
% instructions for your computer telling it how to render the curves rather
% than information about what color goes in what pixel like a jpeg. This
% means that your image will not get pixelated if it's rescaled in any way.
% LaTeX knows how to deal with EPS files, and so should LyX.

print(f, '-depsc', 'figure1.eps')

%This saves the file in the active directory. '-depsc' tells MATLAB to save
% as an EPS color file. There are lots of other commands to save as other
% types of files that you can find on the MATLAB documentation for "print".
% Alternatively, you can save manually by selecting File>Save in the
% figure. You can choose image formats from a dropdown.

%% Plotting: Other types of plots
%MATLAB has loads of other built in imaging options. Some cool/useful ones
% for visualizing the contents of matrices are "spy" and "imagesc". First
% we'll generate a random matrix and see what's in it.

clear %this gets rid of all the variables in our workspace
M = magic(20); %this generates a 20x20 "magic" square where columns and rows have equal sums.
imagesc(M)
colorbar %this draws the color scale so you can see what the colors mean.

%This can help you visualize what a matrix is made of. Perhaps you wouldn't
% notice the structure of a matrix if it wasn't in living color. This is
% also useful for plotting solutions to 2D PDEs. Each element in a matrix
% could represent the height of a solution surface.
%%
%Another tool is spy. Spy just tells you where nonzero elements of a matrix
% are. This is helpful if you wish to check whether your matrix is
% tridiagonal, or something.

clear
M = eye(20); %this creates an eye-dentity matrix.
spy(M)
%% 
%If you have a 2D surface, the "surf" and "mesh" functions are good to use.
% Meshgrid creates a matrix of x and y values so that (x(n),y(n)) is a
% coordinate point in a 2D space. We'll use the sqrt and sin functions to
% create something that looks funky.

[x1,y1] = meshgrid(-8:.1:8);
r = sqrt(x1.^2 + y1.^2);  %notice the '.^2'. This means that squaring should
z1 = sin(r)./r;          % be done element wise on the matrix/vector. NOT a 
                        % matrix power or norm. The ".operator" convention
                        % everywhere implies elementwise operations.
mesh(x1,y1,z1);

%Try rotating the plot around once it's generated.
%%
%Surf basically works the same way and does something similar, just filling
% in the spaces between mesh lines and making the mesh dark. Surf looks
% pretty when you have a low resolution grid. Mesh is better for higher
% resolution grids.

[x2,y2] = meshgrid(-8:.5:8);
r = sqrt(x2.^2 + y2.^2);
z2 = sin(r)./r;
surf(x2,y2,z2);

%% Plotting: Multiple plots
%If you want multiple plots to show up in the same figure, for the sake of
% comparison or whatever, you can use "subplot." subplot(n,m,selected)
% creates a figure with n rows, m columns, and starts with the 3rd plot if
% selected = 3. It counts from the upper left and moves right, then down,
% like reading English. Here I'm going to create a plot consisting of the
% mesh/surf comparison.

subplot(1,2,1)
mesh(x1,y1,z1)
title('This is a mesh plot') %the title function works similarly to xlabel.
subplot(1,2,2)
surf(x2,y2,z2)
title('This is a surf plot')

%This sets the plot's position and fixes the dimensions (1000x400 pixels)
% so that everything looks nice.
set(gcf,'Position',[100 100 1000 400])

%This sets the "colormap" of the 2D figure. MATLAB has a bunch of them
% built in, which can be found by searching "colormap" in documentation.
colormap('Summer')
%% Protips: Built-in functions
%Because of how MATLAB is designed, it is often computationally faster to
% do things in a certain way. MATLAB is an interpreted language, which
% means that coding in it is mostly intuitive and flexible. However, this
% comes at the cost of pure computational speed. MATLAB wastes time trying
% to figure out what types variables are etc. Luckily, a lot of really
% common/useful tasks are included in built-in functions that can be
% accessed and run very quickly, because they are in the form of compiled
% code. Let's check out the difference between interpreted and compiled
% functions.

clear %this empties the workspace
clc %this clears out all the output in your command window
randomNumbers = rand(10^7,1); %creates a row of U(0,1) random numbers

%Let's try summing the elements of this vector manually.
tic %this starts a timer.
total = 0;
for i=1:length(randomNumbers)
    total = total + randomNumbers(i);
end
total
toc %this ends the time and spits out the run time.

%Now let's try summing the elements using MATLAB's sum function, which
% probably looks a lot like ours in terms of source code but is compiled and
% should run more quickly.

tic %restart the timer
sum(randomNumbers)
toc

%Now check out your workspace. There is about an order of magnitude
% difference in compute time here! The lesson here is that if you want to
% do a simple task, chances are that MATLAB has a built in function that
% does it (faster). Knowing what to google may be tricky, so feel free to
% ask around the department!

%Also, tic/toc are kind of simple ways for checking performance. The
% profiler is a much more sophisticated tool. If running MATLAB 2012b or
% newer, select "Run and Time" when running a function to see how long the
% basic components of your program take to run. Search the documentation 
% for "Profiler" if you have an earlier version of MATLAB. The profiler can 
% help you identify trouble spots/bottlenecks.

%% Protips: Vector arithmetic > loops
%If you want to do something to a vector/matrix, like sum its elements or
% performe a finite difference (discretization method for PDEs), there is
% usually a way to do with built-in functions and vector operations.
% Wherever possible, avoid for loops because they are slow. You need a for
% loop for your time integration in the Forward Euler method, because each
% step depends on the step before it. But when this is not the case, let's
% try to use vector operations.

close all
clc
clear

x = -2:0.00001:2; %discretized x values on the interval [-2,2]

tic
%create the discrete function y(x) = x^2 one element at a time
for i = 1:length(x)
   y1(i) = x(i)^2;
end
toc

%Now do the same thing using the vector exponent operator
tic
y2 = x.^2;
toc

subplot(2,1,1)
plot(x,y1)
subplot(2,1,2)
plot(x,y2)
 
%Check your command window. You should see at least an order of magnitude 
% difference here. If the vectors are small, the effect of "overhead" in
% the form of startup/background processes in MATLAB may cause the
% improvement to seem more negligible. In any case, using vector operators
% is a good practice, especially since it is more elegant and will help if
% you decide to scale up your program.

%% Protips: Indexing
%You can do neat things with logical operators and MATLAB's vector
% indexing. Say you want to find all elements in a vector that are
% bigger than some number.

clear
clc
close all

randomNumbers = rand(100,1);
biggerThanIndices = randomNumbers > 0.9;

%"biggerThanIndeces" contains a 100x1 logical array, basically an array of
% true/false answers to the question, "is this element bigger than
% 0.9". We can feed this back into the randomNumbers array to see which
% values these are, specifically.

randomNumbers(biggerThanIndices)

%Check your command window to see all the values in the randomNumbers array
% that are bigger than 0.9. Note that we could do this all in one line:
% randomNumbers(randomNumbers > 0.9)

%%
%Another nifty indexing feature is the "find" function. By itself,
% find(vector) spits out the indices of nonzero elements in the form of a
% logical array. Using logical operators, we can use find to tell us where
% certain values are in a vector. Say we want to find out where the number
% 0.6 is in a linearly spaced vector from 0 to 1.

clear
clc

x = 0:0.01:1;
indices = find(x==0.6)

%Since 0.6 occurs only once, the command window should report that 0.6 only
% occurs in index 61.

x(indices)
