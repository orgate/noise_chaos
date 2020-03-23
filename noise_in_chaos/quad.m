r = [];
n=100000;
%num = 1;

for i=1:n
    r = [r,roots([rand,rand,rand])];
end

length(find(imag(r)==0))/(2*n)

%plot(r(find(imag(r)==0)));
hist(r(find(imag(r)==0)),100);