% Compute R0 = beta/gamma of US,FR,JP,TW
% from SIR model we get:
%
% dI/dt = beta*I*(N-I-R)/N - gamma*I
% dR/dt = gamma*I
% 
% Let dI/dt = I_1 - I_0 = beta*I_0*(N-I_0-R_0)/N-gamma*I_0 
%     dR/dt = R_1 - R_0 = gamma*I
% We can use fsolve to solve
% F1(beta, gamma) = beta*I_0*(N-I_0-R_0)/N-gamma*I_0 - dI/dt = 0
% F2(beta, gamma) = gamma*I - dR/dt = 0
%
% US
tic
N = 329227746; % population (constant)
% import data
DATA = csvread('US.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
days = length(I);

%   
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = 0;
    else
        R_0(k) = X(1)/X(2);
    end    
end
d = 1:days-1;
subplot(4,1,1)
plot(d, R_0)
grid on
legend('US')

% FR
N = 68809894; % population (constant)
DATA = csvread('FR.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
   
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = 0;
    else
        R_0(k) = X(1)/X(2);
    end    
end
subplot(4,1,2)
plot(d, R_0,'c')
grid on
legend('FR')

% JP
N = 124271318; % population (constant)
DATA = csvread('JP.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
days = length(I);
   
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = 0;
    else
        R_0(k) = X(1)/X(2);
    end    
end
subplot(4,1,3)
plot(d, R_0,'y')
grid on
legend('JP')

% TW
N = 23563356; % population (constant)
DATA = csvread('TW.csv');
I = flip(DATA(:,1));
R = flip(DATA(:,2));
   
for k = 1:days-1
    dIdt = I(k+1)-I(k);
    dRdt = R(k+1)-R(k);
    X = fsolve(@(X) IR_equations(X,N,I(k),R(k),dIdt,dRdt),[1,0.01]);
    if X(2) == 0
        R_0(k) = 0;
    else
        R_0(k) = X(1)/X(2);
    end    
end
subplot(4,1,4)
plot(d, R_0,'m')
grid on
legend('TW')
toc
function F = IR_equations(X, N, I, R, dIdt,dRdt)
    F(1,1) = X(1)*I*(N-I-R)/N-X(2)*I-dIdt;
    F(2,1) = X(2)*I-dRdt;
end