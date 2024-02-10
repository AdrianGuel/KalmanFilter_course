%% LKF Guel cortez 2024
% Constant acceleration model
% Ref https://towardsdatascience.com/tuning-q-matrix-for-cv-and-ca-models-in-kalman-filter-67084185d08c

clearvars;
close all;

Ts=1e-2; L=1000;
CA_generative_model;

P=1e4*eye(3);
Q=1e-3*eye(3);
R=1e-2;
alpha=0.001;

x_pred=zeros([3,length(t)]);
x_est=zeros([3,length(t)]);
y_pred=zeros([1,length(t)]);
v_sum=zeros([1,15]);
eps_max=0.1;
Q_scale_factor=10; count=0;

x_pred(:,1)=[3,-1,10]';
for k=2:length(t)
    x_pred(:,k)=Ad*x_pred(:,k-1);
    y_pred(k)=C*x_pred(:,k);

    P_pred=Ad*P*Ad'+Q;
    P_y=C*P_pred*C'+R;
    P_xy=P_pred*C';

    L=P_xy/P_y;
    x_est(:,k)=x_pred(:,k)+L*(y(k)-y_pred(k));
    P=P_pred-L*P_y*L';
    
    %% Innovation based adaptive method 
    % v_sum(end)=y(k)-y_pred(k);
    % v_sum=circshift(v_sum,-1);
    % vk=mean(v_sum);
    % Q=vk*(L*L');

    %% Q-learning
    %Q=(1-alpha)*Q+alpha*L*(y(k)-y_pred(k))*(y(k)-y_pred(k))'*L';
    
    %% scaling factor
    %v_sum(end)=y(k)-y_pred(k);
    %v_sum=circshift(v_sum,-1);
    %alpha_k=(mean(v_sum)-R)/trace(C*P_pred*C');

    %% Bar Shalom 
    %eps= (y(k)-y_pred(k))/P_y;
    %if eps > eps_max
    %    Q = Q*Q_scale_factor;
    %    count = count + 1;
    %elseif count > 0
    %    Q = Q/Q_scale_factor;
    %    count = count- 1;
    %end

    %% Z 
   %if abs(y) > std_scale * std_dev
   %     Q= Q+phi;
   %     count = count + 1;
   %elseif count > 0
   %     Q = Q - phi;
   %     count = count- 1;
   %end
end

subplot(3,1,1)
yyaxis left
plot(t,x(1,:),'b')
hold on
plot(t,x_est(1,:),'k--')
ylabel("position")
yyaxis right
plot(t,abs(x(1,:)-x_est(1,:)))
legend("truth", "estimate", "error")


subplot(3,1,2)
yyaxis left
plot(t,x(2,:),'b')
hold on
plot(t,x_est(2,:),'k--')
ylabel("velocity")
yyaxis right
plot(t,abs(x(2,:)-x_est(2,:)))
legend("truth", "estimate", "error")


subplot(3,1,3)
yyaxis left
plot(t,y,'b')
hold on
plot(t,x_est(3,:),'k--')
ylabel("acceleration")
yyaxis right
plot(t,abs(x(3,:)-x_est(3,:)))
legend("truth", "estimate", "error")
