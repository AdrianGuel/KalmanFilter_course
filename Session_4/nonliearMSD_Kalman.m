%% LKF Guel cortez 2024
% Constant acceleration model
% Ref https://towardsdatascience.com/tuning-q-matrix-for-cv-and-ca-models-in-kalman-filter-67084185d08c

clearvars;
close all;
Ts=1e-2; L=50;
nonlinearMSD;

P=1e5*eye(2);
Q=1e-6*eye(2);
R=1e-5;
alpha=0.03;

x_pred=zeros([2,length(t)]);
x_est=zeros([2,length(t)]);
y_pred=zeros([1,length(t)]);
C=[0,1];

x_pred(:,1)=[0.1,0.2]';
for k=2:length(t)
    x_pred(1,k)=x_pred(1,k-1)+Ts*(x_pred(2,k));
    x_pred(2,k)=x_pred(2,k-1)+(Ts/m)*(-(w^2)*x_pred(1,k-1)^3-gamma*x_pred(2,k-1)+u);
    y_pred(k)=x_pred(1,k);

    Ak=[[1,Ts];[-3*Ts*(w^2)*(x_pred(1,k)^2)/m,-Ts*gamma/m+1]];
    P_pred=Ak*P*Ak'+Q;
    P_y=C*P_pred*C'+R;
    P_xy=P_pred*C';

    L=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+L*(y(k)-y_pred(k));
    P=P_pred-L*P_y*L';

    %% Q-learning
    Q=(1-alpha)*Q+alpha*L*(y(k)-y_pred(k))*(y(k)-y_pred(k))'*L';
end

figure
subplot(2,1,1)
yyaxis left
plot(t,x(1,:),'b')
hold on
plot(t,x_est(1,:),'k--')
ylabel("position")
yyaxis right
plot(t,abs(x(1,:)-x_est(1,:)))


subplot(2,1,2)
yyaxis left
plot(t,y,'b')
hold on
plot(t,x_est(2,:),'k--')
ylabel("velocity")
yyaxis right
plot(t,abs(y-x_est(2,:)))
legend("truth", "estimate", "error")

