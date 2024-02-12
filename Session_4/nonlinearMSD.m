%% Guel Cortez 2024
% Constant acceleration model
% Model simulating the real system

rng(0)

t=0:Ts:L;
x=zeros([2,length(t)]);
y=zeros([1,length(t)]);
u=0;
x(:,1)=[1,0.5]';
R=1e-5;
w=3;m=1; gamma=0.8;
xi_1=1e-6; xi_2=1e-6;

for k=2:length(t)
    x(1,k)=x(1,k-1)+Ts*(x(2,k-1)+sqrt(xi_1)*randn);
    x(2,k)=x(2,k-1)+(Ts/m)*(-(w^2)*x(1,k-1)^3-gamma*x(2,k-1)+u+sqrt(xi_2)*randn);
    y(k)=x(2,k)+sqrt(R)*randn;
end
figure
subplot(2,1,1)
plot(t,x(1,:))
subplot(2,1,2)
plot(t,y)
