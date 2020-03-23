function dx = F(T, X, a, b, c)
% Evaluates the right hand side of the Lorenz system
% x' = sigma*(y-x)
% y' = x*(rho - z) - y
% z' = x*y - beta*z
% typical values: rho = 28; sigma = 10; beta = 8/3;

    d = length(a);
    dx = zeros(d,1);
    for i=1:d
        dx_1=0;
        for j=1:d
            dx_2=0;
            for k=1:d
                dx_3=0;
                for l=1:d
                    dx_3 = dx_3 + c(i,j,k,l)*X(j)*X(k)*X(l);
                end
                dx_2 = dx_2 + dx_3 + b(i,j,k)*X(j)*X(k);
            end
            dx_1 = dx_1 + dx_2 + a(i,j)*X(i);
        end
        dx(i) = dx(i) + dx_1 - X(i).^5;
    end
%    dx = a*X + (b*X)*X;
%    dx = dot(a,X) + dot(dot(b,X),X) + dot(dot(dot(c,X),X),X) - X.^5;

    return
end