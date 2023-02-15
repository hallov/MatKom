clc; clear; close all
%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1;
Tfinal = 1000;
J = 100;

dx = L/(J + 1);
dt = 1e-3/4;

t = 1;
xx = linspace(0,L,J+2 ).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PHYSICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1/7;
v = 1;
xxx = linspace(0,L,J ).';
weigt=exp(-2*xxx()-0.8);
f = @(t,u) 0; %400*abs((sin(t*10/dt))).*weigt; Bara min lek
r =@(c) 0 ; %% negativ reaction ger mer c?? -c/(250*dt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(J,1); 
for i=1:J+2

    if ((xx(i)>=0.45)&&((xx(i)<=0.55)))
        c(i)=1;
   
    end
end
%%c = c0(xx(2:end-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CONSTRUCTING FINITE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = full(spdiags(bsxfun(@times,ones(J,1),[-1 0 1]),[-1 0 1],J,J))/(2*dx);
C=zeros(J-2,1);
E=[1;C;0]/(2*dx)*-v; %V från bondires
B=[1;C;0]/(dx*dx)*d; %D från bondries
D = (full(spdiags(bsxfun(@times,ones(J,1),[1 -2 1]),[-1 0 1],J,J)))/(dx*dx);

M = d*D - v*A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PREPARING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = axes;
p = plot(xx,[1;c;0],'.-');
p.DisplayName = "Pollutant concentration";
ax.YLimMode = 'manual';
xlabel('m')
ylabel('kg/m^3')
ylim([0 2])
title(sprintf("t = %.3e",t));
legend show
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% SOLVING THE PROBLEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while t < Tfinal
    c = c + dt*(M*c+B+E - r(c) + f(t,c)); %B o E är för boundires
    p.YData = [1;c;0]; %% tvingar kanterna är 0 tror inte den som fuckar
    t = t + dt;
    ax.Title.String = sprintf("t = %.3e",t);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clc; clear; close all