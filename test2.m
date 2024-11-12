N=10;
I=eye(N);
H=zeros(N);
Q=zeros(N);
for j = 1 : N
    for k = 1 : N
        H(j,k)=1/(j+k-1);
        if k < j+1
            Q(j,k)=H(j,k);
        end
    end
end
M=I-inv(Q)*H;
p=poly(M);
x=(0:1200)/1000;
y=polyval(p,x);
plot(x,y,'o-',LineWidth=0.5);
grid minor
hold on
plot([x(1),x(1201)],[0,0]);
xlabel('x');
ylabel('y');
legend('p(x)');