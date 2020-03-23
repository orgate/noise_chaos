N = 10000;
R = 3:0.01:4;
x = zeros(N,length(R));
x(1,:) = 0.2;
y = zeros(1,N);
y(1) = 0.01;
d = 0.1;
X = zeros(100,length(R));

for r=1:1%length(R)
    for i=1:N-1
        %x(i+1,r) = R(r)*x(i,r)*(1-x(i,r));
        y(i+1) = y(i)+0.01*(y(i)*sqrt(d) + (d*y(i)*y(i)) - (y(i)^3));
    end
    X(:,r) = hist(x(:,r),100);
end

hist(y,100);
%surf(X);