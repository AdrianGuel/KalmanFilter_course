%% Guel Cortez 2024
% Simple cart model in one dimension
% Generative model
A=[[0,1];[-k_m/m,-b/m]];
B=[0,1/m]';
C=[1,0];
D=0;

Ad=eye(2)+Ts*A;
Bd=Ts*B;
t=0:Ts:L;
x=zeros([2,length(t)]);
y=zeros([1,length(t)]);
u=0;
x(:,1)=[6,3]';
w=3e-5;
v=1e-3;

for k=2:length(t)
    x(:,k)=Ad*x(:,k-1)+Bd*u+w*randn([2,1]);
    y(k)=C*x(:,k)+v*randn;
end
