% LORENZ Function generates the lorenz attractor of the prescribed values
% of parameters rho, sigma, beta
%
%   [X,Y,Z] = LORENZ(RHO,SIGMA,BETA,INITV,T,EPS)
%       X, Y, Z - output vectors of the strange attactor trajectories
%       RHO     - Rayleigh number
%       SIGMA   - Prandtl number
%       BETA    - parameter
%       INITV   - initial point
%       T       - time interval
%       EPS     - ode solver precision
%
% Example.
%        [X Y Z] = lorenz(28, 10, 8/3);
%        plot3(X,Y,Z);

%if nargin<3
%  error('MATLAB:lorenz:NotEnoughInputs','Not enough input arguments.'); 
%end

%if nargin<4
%  eps = 0.000001;
%  T = [0 25];
%  initV = [0 1 1.05];
%end

d=20;
a = randn(d,d);
b = randn(d,d,d);
c = randn(d,d,d,d);

%rho = 28;
%sigma = 10;
%beta = 8/3;
eps = 0.000001;
T = [0 25];
%initV = [0 1 1.05];
init = randn(d,1);
histx = zeros(1,20000);

[T,X] = ode45(@(T,X) F(T, X, a, b, c), T, init);

%for i=1:50
%    init = randn(d,1);
%    %options = odeset('RelTol',eps,'AbsTol',eps*ones(1,d));
%    [T,X] = ode45(@(T,X) F(T, X, a, b, c), T, init);
%    for j=1:10
%        for k=1:length(X)
%            ind = round(X(k,j)*1000+10000);
%            histx(i,ind) = histx(i,ind) + 1;
%        end
%    end
%end
%plot3(X(:,1),X(:,2),X(:,3));
%axis equal;
%grid;
%title('Lorenz attractor');
%xlabel('X'); ylabel('Y'); zlabel('Z');

%x = X(:,1);
%y = X(:,2);
%z = X(:,3);

