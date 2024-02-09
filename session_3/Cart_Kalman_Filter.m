%% LKF Guel cortez 2024
% LKF for cart model


clearvars;
close all;

b=1; m=2; k_m=5; Ts=1e-3; L=100;
Cart_generative_model;

Q=1e-4*eye(2);
R=v^2;
P=1e4*eye(2);
alpha=0.1;

x_pred=zeros([2,length(t)]);
x_est=zeros([2,length(t)]);
y_pred=zeros([1,length(t)]);

x_pred(:,1)=[6,0]';
for k=2:length(t)
    x_pred(:,k)=Ad*x_pred(:,k-1)+Bd*u;
    y_pred(k)=C*x_pred(:,k);

    P_pred=Ad*P*Ad'+Q;
    P_y=C*P_pred*C'+R;
    P_xy=P_pred*C';

    L=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+L*(y(k)-y_pred(k));
    P=P_pred-L*P_y*L';
    
    % if k<3
    % Q=(1-alpha)*Q+alpha*L*(y(k)-y_pred(k))*(y(k)-y_pred(k))'*L';
    % end
end

subplot(2,1,1)
yyaxis left
plot(t,y,'b')
hold on
plot(t,x_est(1,:),'k--')
ylabel("position")
yyaxis right
plot(t,abs(x(1,:)-x_est(1,:)))
ylabel("error")
legend("true","estimation","error")

subplot(2,1,2)
yyaxis left
plot(t,x(2,:),'b')
hold on
plot(t,x_est(2,:),'k--')
ylabel("velocity")
yyaxis right
plot(t,abs(x(2,:)-x_est(2,:)),'k')
ylabel("error")
