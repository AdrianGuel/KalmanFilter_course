%% MSD using UKF
clearvars; close all;
%% Solution to real model to get measurement data
Ts=1e-2;
t=0:Ts:50;
rng(0);
u=zeros(1,length(t));
%u=0.1*sin(t);
w=1e-4; v=1e-3;
kr=1.5; b=0.2; m=1;
x=zeros(2,length(t));
y=zeros(1,length(t));
x(:,1)=[4;0.1];
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
ffactor=0.03;

P= chol(diag([0.5 0.18]));
Q=chol(diag([1e-2,1e-4]));
R=chol(1e-3);

x_pred=zeros([2,length(t)]);
x_est=zeros([2,length(t)]);
y_pred=zeros([1,length(t)]);
x_est(:,1)=[1,0.04]';

%% UKF
for k=2:length(t)
    X = x_est(:,k-1).*ones([1,2*L+1]) + eta*[zeros([L,1]), P, -P];% Sigma points

    % prediction
    Y = F(X,u(k-1),m,b,kr,Ts,L); %Dynamics
    x_pred(:,k)=Y*Wm'; % Dynamics mean value
    % prediction Cov update
    Xs=(Y(:,2:end)-x_pred(:,k))*sqrt(Wc(2));
    C=[Xs,Q];
    S_minus = qr(C',"econ");
    x_dev = Y(:,1) - x_pred(:,k);
    Sx=cholupdate(S_minus,sqrt(Wc(1))*x_dev);
    % system output
    Z= H(X,L); % system output
    y_pred(k)=Z*Wm';% system mean output
    %output cov update
    Zs = sqrt(Wc(2))*(Z(2:end) - y_pred(k));
    N = [Zs, R];
    Sy = qr(N',"econ"); 
    y_dev = Z(1) - y_pred(k);
    Sy = cholupdate(Sy, sqrt(Wc(1))*y_dev);  

    P_xy= Xs*Zs' + Wc(1)*(x_dev*y_dev');% correlation between dynamics and output
    %% Kalman Gain
    K = (P_xy/Sy')/Sy;   % (dim: 3x1) Kalman gain calculation [Ref1, Eqn. 27]
    x_est = x_pred + K*(y(k) - y_pred(k));      % (dim: 3x1) State estimation [Ref1, Eqn. 28]
    U = K*Sy;  % (dim: 3x1) [Ref1, Eqn. 28]
    S = cholupdate(Sx, U, "-");  % (dim: 3x3) State covariance estimation [Ref1, Eqn. 29]

    % Q-adaptive
    %% Q learning
    %Q=(1-ffactor)*Q+ffactor*K*(y(k)-y_pred(k))*(y(k)-y_pred(k))'*K';    
end


%% Plotting

subplot(2,1,1)
yyaxis left
plot(t,y,'k')
hold on
plot(t,x_est(1,:),'b.-')
ylabel('$x$','Interpreter','latex',FontSize=12)
yyaxis right 
plot(t,abs(y-x_est(1,:)),'r')
legend('y','x1 estimated','error')

subplot(2,1,2)
yyaxis left
plot(t,x(2,:),'k')
hold on
plot(t,x_est(2,:),'b.-')
ylabel('$\dot{x}$','Interpreter','latex',FontSize=12)
yyaxis right 
plot(t,abs(x(2,:)-x_est(2,:)),'r')
legend('x2 ','x2 estimate ','error')

%% Fuctions
function aux=model(x,F,m,b,k,w)
aux=[x(2)+w*randn;
    -(k/m)*x(1)-(b/m)*x(2)+F+w*randn];
end

function Y=F(X,u,m,b,k,Ts,L)
Y=zeros([L,2*L+1]);
Ad=eye(2)+Ts*[[0,1];[-k/m,-b/m]]; Bd=Ts*[0,1/m]';
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
