maxT=200;
x=zeros(4,maxT);
x(:,1)=[0;0.1;0;0.1];
T=0.1;
A=[1 T 0 0; 0 1 0 0;...
    0 0 1 T; 0 0 0 1];
B= [T^2/2 0; T 0; 0 T^2/2; 0 T];
covar = [2, 0.75; 0.75, 1];
corr = chol(covar,'lower');

for i=1:100
    for k=2:maxT
        x(:,k)=A*x(:,k-1)+B*corr*randn(2,1);
    end
    
    plot(x(1,:),x(3,:));
    hold on
end
