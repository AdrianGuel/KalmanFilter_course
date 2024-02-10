%% Guel Cortez 2024
% Constant acceleration model
% Model simulating the real system

rng(0)

Ad=[[0,Ts,0.5*Ts^2];[0,1,Ts];[0,0,1]];
Bd=0;
C=[0,0,1];
D=0;


t=0:Ts:L;
x=zeros([3,length(t)]);
y=zeros([1,length(t)]);
u=0;
x(:,1)=[0,-50,10]';
Q=1e-6*eye(3);
R=1e-3;

for k=2:length(t)
    x(:,k)=Ad*x(:,k-1)+sqrt(Q)*randn([3,1]);
    y(k)=C*x(:,k)+sqrt(R)*randn;
end
