% Networked Control Systems Project

clear, close all

%% Define adjacency graphs and vigilancy and susceptability rate functions

n = 46;
p = 0.2;
d = 10;

% USING ERDOS-ROGAN FOR PHYSICAL GRAPH FOR NOW (NEED TO FIGURE OUT THE
% JAPANESE POPULATION DATABASE)
% Beta = ER(n,p)+eye(n); % Epidemics connectivity graph
% degree = sum(Beta,2);
% Beta = diag(1./(2*degree))*ones(n).*Beta; % Normalize such that Beta*1n = 0.5*1n

% Import adjacency matrix from Japanese database
load("Beta.mat")
Beta = 0.5*adj;

% Using Watts-Strotgatz model for opinion graph
Omega = WS(n,d,p); % Opinions connectivity graph
degree = sum(Omega,2);
Omega = diag(1./degree)*ones(n).*Omega; % Normalize such that Omega*1n = 1n

LW = Laplacian(Omega);

% Control parameters
epsilon = 0;
tau = 0;

Gamma = @(xOk) diag((0.4-epsilon)*ones(length(xOk),1) - (0.4-epsilon)*xOk); % Rate of a vigilant person to become susceptible
Theta = @(xOk) diag((0.2+epsilon)*ones(length(xOk),1) + (0.3-epsilon)*xOk); % Rate of a susceptible person to become vigilant

Phi = 0.1*eye(n); % Impact of becoming infected on the opinion

% Load recovery rate matrix from japanese hospitals data
% Delta = 0.2*eye(n); % Recovery rate of infected population (CAN BE DERIVED FROM JAPANESE MINISTRY OF HEALTH DATA)
load('Delta.mat')
Delta = diag(Delta);

PsiFun = @(o) Theta(o)./(Gamma(o)+Theta(o)); % Calculate Psi for Condition numbers

ob = 0:0.01:0.8; % Opinion state boundaries

for j = 1:length(ob)
    PSI(j) = PsiFun(ob(j));
end

% Find psi and psihat
PsiHat = max(PSI);
Psi = min(PSI);

R0Vbar = eigs(eye(n) - Delta + Beta - PsiHat*Beta, 1, 'lm');
R0V = eigs(eye(n) - Delta + Beta - Psi*Beta, 1, 'lm');

if R0Vbar > 1
    disp('R0V_bar > 1, there is at least one endemic equilibrium')
else
    disp('R0V_bar <= 1, there may or may not be an endemic equilibrium')
end

if R0V <= 1
    disp('R0V <= 1, there is a healthy globally stable equilibirum')
else
    disp('R0V > 1, there is no healthy globally stable equilibrium')
end

ROVSIS = eigs(eye(n) - Delta + Beta, 1, 'lm');


%% Simulation

T = 49; % Final time

% Intitial conditions
xOk = 0.8*rand(n,1);
xSk = 0.1+0.3*rand(n,1);
xIk = 0.1+0.5*rand(n,1);
xVk = ones(n,1) - xIk - xSk;

xk = [xSk;
      xIk;
      xVk;
      xOk];

% Sim
for i = 1:T
    xk1 = SIVO(xk(:,i),n,Omega,Beta,Phi,Delta,Gamma,Theta,tau);
    xk = [xk xk1];
end

xSk = xk(1:n,:); % Susceptible
xIk = xk(n+1:2*n,:); % Infected
xVk = xk(2*n+1:3*n,:); % Vigilant
xOk = xk(3*n+1:4*n,:); % Opinion

%% Plot

% des_prefs = pref_names';
des_prefs = {'Kanagawa','Aichi','Hyogo','Miyagi','Tottori'};

choice = find(matches(pref_names,des_prefs));

linethick = 1;

% plot S
figure
hold on
plot(xSk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xSk,1)','k--','LineWidth',linethick);
title('Susceptible Population Evolution')
if length(des_prefs) <= 5
    legend([des_prefs {'Average'}])
else
    legend(p2,'Average')
end
xlabel('Time Step (k)')
ylabel('Susceptible Population')
movegui('northwest')

% plot I
figure
hold on
plot(xIk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xIk,1)','k--','LineWidth',linethick);
title('Infected Population Evolution')
if length(des_prefs) <= 5
    legend([des_prefs {'Average'}])
else
    legend(p2,'Average')
end
xlabel('Time Step (k)')
ylabel('Infected Population')
movegui('northeast')

% plot V
figure
hold on
plot(xVk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xVk,1)','k--','LineWidth',linethick);
title('Vigilant Population Evolution')
if length(des_prefs) <= 5
    legend([des_prefs {'Average'}])
else
    legend(p2,'Average')
end
xlabel('Time Step (k)')
ylabel('Vigilant Population')
movegui('southwest')


% plot O
figure
hold on
plot(xOk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xOk,1)','k--','LineWidth',linethick);
title('Population Opinion Evolution')
if length(des_prefs) <= 5
    legend([des_prefs {'Average'}])
else
    legend(p2,'Average')
end
xlabel('Time Step (k)')
ylabel('Population Opinion ')
movegui('southeast')

% All four
figure
subplot(2,2,1)
hold on
plot(xSk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xSk,1)','k--','LineWidth',linethick);
xlabel('Time Step (k)')
ylabel('Susceptible Population')
subplot(2,2,2)
hold on
plot(xIk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xIk,1)','k--','LineWidth',linethick);
xlabel('Time Step (k)')
ylabel('Infected Population')
subplot(2,2,3)
hold on
plot(xVk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xVk,1)','k--','LineWidth',linethick);
xlabel('Time Step (k)')
ylabel('Vigilant Population')
subplot(2,2,4)
hold on
plot(xOk(choice,:)','LineWidth',linethick)
p2 = plot(mean(xOk,1)','k--','LineWidth',linethick);
xlabel('Time Step (k)')
ylabel('Population Opinion')
if length(des_prefs) <= 5
    legend([des_prefs {'Average'}],'Location','northeast')
else
    legend(p2,'Average','Location','northeast')
end
sgtitle(['SIVO Epidemics Model w/ ROV: ' num2str(R0V) ' and ROV_{bar}: ' num2str(R0Vbar)])


%% Functions

function [xk1] = SIVO(xk,n,Omega,Beta,Phi,Delta,Gamma,Theta,tau)

I = eye(n);
One =  ones(n,1);

xSk = xk(1:n); % Susceptible
xIk = xk(n+1:2*n); % Infected
xVk = xk(2*n+1:3*n); % Vigilant
xOk = xk(3*n+1:4*n); % Opinion

xIk1 = xIk + (diag(One - xVk - xIk))*Beta*xIk - Delta*xIk; % Update infected
xVk1 = xVk + Delta*xIk + Theta(xOk)*(One - xVk - xIk) - Gamma(xOk)*xVk; % Update vigilant
% xOk1 = Phi*xIk + (I - Phi)*(xOk + (I - I*xOk)*-Laplacian(Omega)*xOk); % Update opinion
xOk1 = tau+(1-tau)*(Phi*xIk + (I - Phi)*((I - I*xOk)*Omega + xOk)*xOk); % Update opinion
xSk1 = One - xVk1 - xIk1; % Update susceptible

xk1 = [xSk1;
       xIk1;
       xVk1;
       xOk1];

end

function LW = Laplacian(W)

% Laplacian =  I - W

n = length(W);

LW = eye(n) - W;

end

function A = ER(n,p)
%A = ER(n,p) generates an Erdos–Renyi random graph with network size n and 
% link probability p

% When you write functions, it is often a good idea to put some checks,
% just to verify that inputs are correctly defined (e.g., that p is between
% 0 and 1, n positive integer,...)

A=zeros(n); 

% Remark: How to simulate an event (in our case, generate a link) that occurs with a certain probability p? Generate a random number and if the random number is less than or equal to p, than the event occurs

for i=1:n
    for j=i+1:n
        if rand<=p
            A(i,j)=1;
            A(j,i)=1;
        end
    end
end

% We can do it more efficiently as: A=(rand(n)<=p); A=triu(A,1)+triu(A,1)'

end

function A = ring(n,m)
%A = ring(n,m) generates a ring with network size n and m number of
%connections (per side)

A=zeros(n); 

for i=1:m
    A=A+diag(ones(n-i,1),i)+diag(ones(n-i,1),-i)+diag(ones(i,1),n-i)+diag(ones(i,1),-n+i); %we create the diagonal bands
end

end

function A = WS(n,d,p)
m=d/2;
A=ring(n,m); 
for i=1:n
    for j=i+1:i+m 
        %I create an additional function to take care for links of the form (n,1)
        if A(i,modn(j,n))==1
            %check if we rewire it
            if rand<=p
                % rewire
                flag=0;
                while flag==0
                    S=randi(n);
                    if A(i,S)==0 && S~=i  %arc not present and not self loop
                        A(i,modn(j,n))=0; %cancel existing arc
                        A(modn(j,n),i)=0; %in both directions
                        A(i,S)=1; %create new arc
                        A(S,i)=1;
                        flag=1; %exit the loop
                    end  
                end
                %no instructions in the "else" are needed, since we keep the original arc in A
            end
        end
    end
end
end

function y = modn(x,n)
%modn(x,n) computes x modulus n (the reminder of the division x:n), but it returns n instead of z is the
%reminder is 0
y=mod(x,n);
if y==0
    y=n;
end
end