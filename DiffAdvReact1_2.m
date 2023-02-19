clc; clear; close all
%%%%%%%%%% SIMULATION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
L = 1;
Tfinal = 100;
J = 100;

dx = L/(J + 1); 
dt = 1e-3/4;

t = 1;
xx = linspace(0,L,J + 2).';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PHYSICAL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
d = 1/10;
v = 2;
f = @(t,u) 0*u; %% Ingen nytt ? %%Vi skickar in c så vad skiljer f från r?
r =@(c) 0; %% negativ reaction ger mer c?? %%Vid positiv r går partiklarna in i reaktioner och 
% försvinner (tolkningsfråga, annars ändra tecken i uppdateringen där nere för omvänd tolkning)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% INITIAL CONDITION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=zeros(J,1); 
for i=1:J+2

    if ((xx(i)>=0)&&((xx(i)<=0.2)))
        c(i)=1;
   
    end
end
%%c = c0(xx(2:end-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% CONSTRUCTING FINITE DIFFERENCES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = full(spdiags(bsxfun(@times,ones(J,1),[-1 0 1]),[-1 0 1],J,J))/(2*dx);
C=zeros(J-2,1);
B=[1;C;0]; %%Tror har med ska hjälpa med boundires
D = (full(spdiags(bsxfun(@times,ones(J,1),[1 -2 1]),[-1 0 1],J,J)))/(dx*dx);

M = d*D - v*A;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PREPARING PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax = axes;
p = plot(xx,[1;c;0],'.-'); %%Gränsvärdena här ändras direkt vid första beräkningen men måste finnas för att fylla grafen?
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
    c = c + dt*(M*c+B - r(c) + f(t,c));
    p.YData = [1;c;0]; %% tvingar kanterna är 0 tror inte den som fuckar %%Tror inte detta påverka beräkningarna utan bara ploten
    t = t + dt;
    ax.Title.String = sprintf("t = %.3e",t);
    drawnow
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%clc; clear; close all