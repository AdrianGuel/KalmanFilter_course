%% MSD using UKF
clearvars; close all;
%% Solution to real model to get measurement data
Ts=1e-2;
t=0:Ts:10;
rng(0);
%F=zeros(1,length(t));
u=0.1*sin(t);
w=1e-4; v=1e-3;
kr=1.9; b=0.5; m=1;
x=zeros(2,length(t));
y=zeros(1,length(t));
x(:,1)=[1;0.1];
for j=2:length(t)
    x(:,j)=x(:,j-1)+Ts*model(x(:,j-1),u(j-1),m,b,kr,w);
    y(j)=x(1,j)+v*randn;
end

%% UKF params
L=2;
Wm=zeros(1,2*L+1);
Wc=zeros(1,2*L+1);
alpha=0.99; beta=2; lambda=L*(alpha^2-1); eta=sqrt(L+lambda);
Wm(1)=lambda/(L+lambda); Wc(1)=(lambda/(L+lambda))+1-alpha^2+beta;
Wm(2:end)=1/(2*(L+lambda));
Wc(2:end)=1/(2*(L+lambda));
ffactor=0.003;

P=1e-6*eye(L,L);
Q=1e-3*eye(L,L);
R=1e-2;

x_pred=zeros([2,length(t)]);
x_est=zeros([2,length(t)]);
y_pred=zeros([1,length(t)]);
x_est(:,1)=[0.1,0.2]';

%% UKF
for k=2:length(t)
    sP = chol(P,'lower'); % square root of P
    X = x_est(:,k-1).*ones([1,2*L+1]) + eta*[zeros([L,1]), sP, -sP];% Sigma points

    % prediction
    Y = F(X,(k-1),m,b,kr,Ts,L); %Dynamics
    x_pred(:,k)=Y*Wm'; % Dynamics mean value
    Xs = (X(:,2:end) - x_est(:,k-1).*ones([1 2*L]))*sqrt(Wc(2));
    Xs1 = X(:,1) - x_est(:,k-1);
    P_pred = Xs*Xs' + Wc(1)*(Xs1*Xs1')+Q; %Dynamics Covariance
    Z= H(X,L); % system output
    y_pred(k)=Z*Wm';% system mean output
    Zs = (Z(:,2:end) - y_pred(k)*ones([1 2*L])) * sqrt(Wc(2));
    Zs1 = Z(:,1) - y_pred(k);
    P_y = Zs*Zs' + Wc(1)*(Zs1*Zs1')+R; %system output covariance
    
    % estimation/update equations
    P_xy= Xs*Zs' + Wc(1)*(Xs1*Zs1');% correlation between dynamics and output
    K_k=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+K_k*(y(k)-y_pred(k));
    P=P_pred-L*P_y*L';

    % Q-adaptive
    %% Q learning
    Q=(1-ffactor)*Q+ffactor*K_k*(y(k)-y_pred(k))*(y(k)-y_pred(k))'*K_k';    
end

%% Fuctions
function aux=model(x,F,m,b,k,w)
aux=[x(2)+w*randn;
    -(k/m)*x(1)-(b/m)*x(2)+F+w*randn];
end

function Y=F(X,u,m,b,k,Ts,L)
Y=zeros([L,2*L+1]);
Ad=eye(2)-Ts*[[0,1];[-k/m,-b/m]]; Bd=Ts*[0,1/m]';
    for k=1:2*L+1
        Y(:,k)=Ad*X(:,k)+Bd*u;
    end
end

function Z=H(X,L)
Z=zeros([1,2*L+1]);
    for k=1:2*L+1
        Z(k)=X(1,k);
    end
end
