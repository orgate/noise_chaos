%rho = [0:50];
rho = logspace(0,3,50);
sigma = [10:10];
beta = [8/3:8/3];
X_array = zeros(length(rho), 100);
%Y_array = zeros(1, 100000);
%Z_array = zeros(1, 100000);
xmax = zeros(length(rho),1);
xmin = zeros(length(rho),1);

for rhoi=1:length(rho)
    for sigmai=1:length(sigma)
        for betai=1:length(beta)
            [X Y Z] = lorenz(rho(rhoi), sigma(sigmai), beta(betai), [0 1 1.05], [0 2500], 0.000001);
            xmax(rhoi) = max(X);
            xmin(rhoi) = min(X);
            X_array(rhoi,:) = hist(X,100);
%            [hx, x] = hist(X,100);
%            [hy, y] = hist(Y,100);
%            [hz, z] = hist(Z,100);
%            for i=1:100
%                indx = round(x(i)*1000+50000)+1;
%                X_array(indx) = X_array(indx)+hx(i);
%            end
        end
    end
end

surf(X_array, del2(X_array));
colormap;
legend('show');
%hist(X_array);