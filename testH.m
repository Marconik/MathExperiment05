num=[10,50,100,200,500,1000];
for i = 1 : length(num)
    N=num(i);
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
    T=eig(M);
    fprintf("%.18f\n",max(abs(T)));
end