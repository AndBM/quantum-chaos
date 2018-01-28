% Note that total spin is both F and j.

% Check the kappa=0 results, maybe fine if they don't match article since
% article includes experimental stuff.

% Just make the figures and finish the report.

% kappa=0 results don't match old results. It seems like the push is acting
% differently?? Also the chaos measure seems to be off for the old data.
% You could check the amplitude in f(t). Or maybe just accept that the new
% stuff is correct?

% Allow interpolation and free time resolution?

% Check if kappa != 0 gives mixed states

% Sometimes Tr(rho^2) > 1, this is numerical effect. It's a pure state all
% the way.

% I think dW now better implemented, but still pretty non-differentiable.
% Maybe smoothen out or something?

% To do:
% X TEST MEASUREMENT IMPLEMENTATION
% - Experiment with different gauss widths?
% X Implement strict chaos signatures?
% X Implement weak measurement
% X Test relation between measurement strength and chaos
% X Check with husimi dists that the spin length measure is correct.
    %- Husimi checks out, but lengths don't...
% - Stuff that Mølmer said
% - randomData0 is strangely homogeneous.
% - Spin-½ not constant for non-zero kappa
% - Spin coherence test shows states pretty different from coherent states

function propagator()
set(0, 'DefaultTextInterpreter', 'latex');
set(0, 'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultAxesFontSize',11);
set(0, 'defaultAxesTitleFontSizeMultiplier',1.5);
set(0, 'defaultLegendInterpreter','latex');

tic
close all;

%------------------------- Settings
% Hamiltonian dynamics
F = 3;                  % Spin
p = 0.99;               % Angle
k = 2;                  % Degree of chaos
tau = 1;                % "Pushing" period
eta = 1;                % Detector efficiency
kappa = 0;           % Probe field strength

% Numerics
endt = 41;              % Timespan (or number of periods if tau = 1)
res = 100;              % Resolution of plot
tres = 'Free'; %endt*100+1;       % Time resolution of ODE
                        % (if tres is 'Free' ODE45 chooses its own resolution)
                        % (this can lead to strange animations)
                        % (this messes with dW)
calcNew = 0;            % Calculate new data? 1 yes, 0 no.
folder = 'randomData0/';% Load from folder if no new calculation
N = 100;    % Number of husimi-distributions to calculate.
                        % If (0 < N < tres) is false, N becomes tres (all are
                        % calculated)
dWstep = 1e-5;
seed = 4;               % Set random seed.
rng(seed,'twister')

% Set initial state
Fs1 = [0.7, 0.7, -0.16];         % Order islands
Fs2 = [-0.94, -0.31, -0.16];   % Chaos sea 2
% Fs = [-0.94, 0.31, -0.16];    % Chaos sea (not)
Fs = [0, -0.99, -0.16];       % Major island
rho = setInit(F, Fs(1), Fs(2), Fs(3));
rho01 = setInit(F, Fs1(1), Fs1(2), Fs1(3));
rho02 = setInit(F, Fs2(1), Fs2(2), Fs2(3));

%-------------------------

% Propagate according to the relevant master eq.

%TEST IF SPIN-½ IS CONSTANT
rho0 = setInit(F, 0, 0, 1);
disp('Propagating...')
[rhos, t] = propagate(F, rho01, tres, p, k, tau, endt, eta, kappa, dWstep);
toc
wienerCheck(t,dWstep);
spinvec(rhos, endt)
%testCoherence(rhos,t)
% Calculate relevant Husimi distributions and save everything
disp('Calculating Husimi distributions...')
Qs = husimi(rhos, res, t, F, N);
animator(Qs,res,endt,tau);
%meanplot(Qs,res);
% toc
% save(datestr8601);

  
% else
%     % Load data from a previous calculation
%     figure
%     hold on
%     files = dir(folder);
%     for i = 1:10
%         filename = [folder files(i+2).name];
%         load(filename)
%         %spinvec(rhos, endt)
%     end
%     legend('\kappa = 0', '\kappa = 0.05', '\kappa = 0.25', '\kappa = 0.5',...
%            '\kappa = 0.75', '\kappa = 1', 'Location', 'Northeast')
%     title('Purity measure, chaos sea')
% end
% else
%     % Load data from a previous calculation
%     figure
%     for i = [1 2 6 11 16 21]
%         filename = ['kappa' num2str((i-1)*5)];
%         load(filename)
%         spinvec(rhos, endt, kappa(i))
%     end
%     legend('\kappa = 0', '\kappa = 0.05', '\kappa = 0.25', '\kappa = 0.5',...
%            '\kappa = 0.75', '\kappa = 1', 'Location', 'East')
%     title('Purity measure, order islands')
% end


toc
%------------------------- Main function end.

function testCoherence(rhos,t)
% Subtract ket*bra from rho, should be zero if spin coherent state

F = (size(rhos,1)-1)/2;
N = size(rhos,3);
rhoDiff = zeros(1,N);
cohRho = zeros(size(rhos));

for i = 1:N
    X = trace(rhos(:,:,i)*Sx(F));
    Y = trace(rhos(:,:,i)*Sy(F));
    Z = trace(rhos(:,:,i)*Sz(F));
    [theta, elev] = cart2sph(X, Y, Z);
    phi = pi/2 - elev;
    ket = expand(theta, phi, F);
    bra = ket';
    cohRho(:,:,i) = ket*bra;
    rhoDiff(i) = sum(sum(abs(rhos(:,:,i) - ket*bra)))/sum(sum(rhos(:,:,i)));
end

cohQ = husimi(cohRho, 100, t, F, 0);

figure
animator(cohQ, 100, 41, 1)
figure
plot(rhoDiff)

function plotJx(seed, rhos, t, F)
% Reset random number generator and recreate dWs
rng(seed,'twister');
dWs = randn(size(t));

Jxs = zeros(size(t));
Sxs = Jxs;

% Define Jx
for i = 1:length(t)
    avSx = trace(rhos(:,:,i)*Sx(F));
    Sxs(i) = avSx;
    Jxs(i) = avSx + dWs(i);
end

figure; hold on;
plot(t, Sxs, '--')
plot(t, Jxs)
xlabel('t')
legend('<S_x>','J_x')

function covariance(F, N, rhos)
covars = zeros(3, 3, N);
eigvals = covars;
eigvecs = covars;

Sops = Sx(F);
Sops(:,:,2) = Sy(F);
Sops(:,:,3) = Sz(F);

for k = 1:N
    for i = 1:3
        for j = 1:3
            avgOp = 1/2*( Sops(:,:,i)*Sops(:,:,j)...
                        + Sops(:,:,j)*Sops(:,:,i) );
            covars(i,j,k) = trace( rhos(:,:,k)*avgOp );
        end
    end
    [eigvecs(:,:,k), eigvals(:,:,k)] = eig(covars(:,:,k));
end
disp('lol')

function spinvec(rhos, endt)
F = (size(rhos,1)-1)/2;
N = size(rhos,3);
lengths = zeros(1,N);
spinvecs = zeros(3,N);

for i = 1:N
    spinvec = [trace(rhos(:,:,i)*Sx(F))
               trace(rhos(:,:,i)*Sy(F))
               trace(rhos(:,:,i)*Sz(F))];
    spinvecs(:,i) = spinvec;
    lengths(i) = norm(spinvec);
end

if max(max(imag(spinvecs))) > 1e-12
    disp('Warning: "spinvecs" has significant imaginary parts (>1e-12)')
end
spinvecs = real(spinvecs);
vars = sum(sum(var(spinvecs')));
disp(vars)

for i=1:length(rhos)
    toplot(i)=trace(rhos(:,:,i)^2);
end

figure; hold on
X = linspace(0,endt,N);
plot(X,lengths)
xlabel('Period number')
ylabel('Measure')
%ylabel('|(\langle F_x \rangle, \langle F_y \rangle, \langle F_z \rangle) |')
%axis([0 endt 0 F+0.1])

figure; hold on;
plot(X,spinvecs(1,:),'b-')
plot(X,spinvecs(2,:),'k-.')
plot(X,spinvecs(3,:),'r--')

legend('Sx','Sy','Sz')

figure; hold on
plot(real(toplot))
xlabel('t')
ylabel('Tr($\rho^2$)')
axis([0 length(toplot) 0 max(real(toplot))+0.1]);

function animator(Qs,res,endt,tau)
% Animate Q on a sphere
figure
set(gcf, 'Position', [53, 274, 1155, 423])
for i = 1:size(Qs,3)
    % Show the pushes
    subplot(1,3,1)
    t = linspace(-tau,0);
    plot(t,f(t + i*endt/size(Qs,3), tau))
    axis([-tau, tau, -0.1, 20]);
    
    % One side of the sphere
    subplot(1,3,2)
    sphereplot(Qs(:,:,i), res, 0);
    xlabel('$F_x/F$'); ylabel('$F_y/F$'); zlabel('$F_z/F$');
    titlestr = ['Period ' num2str(round(i/(size(Qs,3)/endt)))];
    title(titlestr)
    
    % The other side of the sphere
    subplot(1,3,3)
    sphereplot(Qs(:,:,i), res, 180);
    xlabel('$F_x/F$'); ylabel('$F_y/F$'); zlabel('$F_z/F$');
    titlestr = ['Period ' num2str(round(i/(size(Qs,3)/endt)))];
    title(titlestr)
    pause(0.1)
end
    
function rho0 = setInit(F, Fx, Fy, Fz)
% Set initial state as spin coherent state centered at (Fx, Fy, Fz)/F
% Precalculate for speed:
binomial = zeros(2*F+1, 1);
for m = -F:F
    binomial(index(F,m)) = nchoosek(2*F, F+m);
end

[phi, elevation, ~] = cart2sph(Fx,Fy,Fz);
theta = pi/2 - elevation;
ket = expand(theta, phi, F, binomial);
bra = ket';
rho0 = ket*bra;

function meanplot(Qs, res)
Q = mean(Qs,3);

figure
%set(gcf,'position',[3102 1143 662 316])

subplot(1,2,1)
sphereplot(Q, res, 0);
title('$F_y < 0$'); xlabel('$F_x/F$'); ylabel('$F_y/F$'); zlabel('$F_z/F$');
subplot(1,2,2)
sphereplot(Q, res, 180);
title('$F_y > 0$'); xlabel('$F_x/F$'); ylabel('$F_y/F$'); zlabel('$F_z/F$');

function multiplot(Qs, res, tres, N)
sample = floor(linspace(1,tres,N));

figure
for i = 1:N
    Q = Qs(:,:,sample(i));
    subplot(2,N,2*i-1)
    sphereplot(Q, res, 0);
    subplot(2,N,2*i)
    sphereplot(Q, res, 180);
end

function h = sphereplot(Q, res, viewangle)
% Plot Q on a sphere
[thetas, phis] = makeangles(res);
elevation = pi/2-thetas;
[Phi, Elev] = meshgrid(phis, elevation);
[x, y, z] = sph2cart(Phi, Elev, 1);
h = surf(x, y, z, Q);
view(viewangle,0);
daspect([1 1 1]); axis tight; shading interp; colormap(jet);

function Q = husimi(rhos, res, t, F, N)
% Calculate Husimi distribution
[thetas, phis] = makeangles(res);

% Make sure N has proper value, otherwise calculate all
if not(0 < N && N < length(t))
    N = length(t);
end
Q = zeros(length(thetas), length(phis), N);
sampleInd = floor(linspace(1, length(t), N));

% Precalculate for speed:
binomial = zeros(2*F+1, 1);
for m = -F:F
    binomial(index(F,m)) = nchoosek(2*F, F+m);
end

% Run through all coordinates and get a Q-value for each rho
for k = 1:N
    disp(k/N)
    for i = 1:length(thetas)
        for j = 1:length(phis)
            ket = expand(thetas(i), phis(j), F, binomial);
            bra = ket';
            rho = rhos(:,:,sampleInd(k));
            Q(i,j,k) = (2*F + 1)/(4*pi)* bra*rho*ket;
        end
    end
end

% Q is complex, should be real
% Ignore if negligible, otherwise make warning.
if max(max(max(imag(Q)))) > 1e-12
    disp('Warning: Q has significant imaginary parts (>1e-12)')
end
Q = real(Q);

function [thetas, phis] = makeangles(res)
% Get angles of correct resolution
thetas = linspace(0, pi, ceil(res/2));
phis = linspace(0, 2*pi, res);

function [rhos, t] = propagate(F, rho0, tres, p, k, tau, endt, eta, kappa, dWstep)
% Allow both set and free time resolution
if strcmp(tres,'Free')
    tspan = [0 endt];
else
    tspan = linspace(0,endt,tres);
end

% Generate Wiener process
Wtimes = 0:dWstep:endt;
dW = dWstep^2*randn(size(Wtimes)); % Multiply by dt^2 and divide by dt.
Wt = griddedInterpolant(Wtimes, [0 cumsum(dW(2:end))]);
%[Wt; 0 cumsum(dW(2:end))]';
% W_0 = 0, so start at i = 2

% Solve the differential equation for the density
odefun = @(t,y) dRhodt(t, y, F, p, k, tau, eta, kappa, Wt);
opts = odeset('RelTol',1e-6,'AbsTol',1e-6);

y0 = [rho0(:)]; % Shape into vector...
[t, y] = ode45(odefun, tspan, y0, opts);
rhos = reshape(y.',2*F+1,2*F+1,[]); %... and back into matrix

function result = dRhodt(t, y, F, p, k, tau, eta, kappa, Wt)
% Differential equation for density

% Interpolate Wiener process
dW = Wt(t);

% Reshape y into matrix
rho = reshape(y, 2*F+1, 2*F+1);
rho = repair(rho); % Fix residue from approximations

simple = -1i*( H(p,t,tau,F,k) * rho - rho * H(p,t,tau,F,k) );
               
decoherence = kappa*( operatorD( Sy(F), rho )...
                    + operatorD( Sy(F), rho )...
                    + operatorD( Sx(F), rho ));
                
backaction = sqrt(eta*kappa)*( operatorH( Sz(F), rho)...
                             + operatorH( Sy(F), rho)...
                             + operatorH( Sx(F), rho))*dW;
               
master = simple + decoherence + backaction;

% Shape the result back into a matrix
result = master(:);

function result = H(p,t,tau,F,k)
% The article has the same sign on each. Why does it only work with
% different signs?
result = p*f(t,tau)*Sy(F) + k/(2*F*tau)*Sx(F)^2;

function rho = repair(rho)
% Handle inconsistencies in density after approximations in ODE45
n = length(rho);
rho(1:n+1:n*n) = real(diag(rho)); % Remove imaginary residue
rho = 1/trace(rho) * rho; % Renormalize rho

function result = operatorD(a, rho)
result = a*rho*a' - 1/2*( a'*a * rho + rho * a'*a );

function result = operatorH(a, rho)
result = a*rho + rho*a' - trace( a*rho + rho*a' )*rho;

function result = generatedW(t)
global t_prev;
dt = t - t_prev;
t_prev = t;
result = dt*randn(); % Multiply by dt^2 and divide by dt.

function result = f(t,tau)
% "Push" the top with f
a = 1; % Amplitude
c = tau/40; % Width
x = mod(t,tau);

% Gaussian is in the middle of periodic interval
%result = 1/(c*sqrt(2*pi)) * exp(-x.^2/(2*c^2));
result = 1/(c*sqrt(2*pi)) * exp(-(x-tau/2).^2/(2*c^2));
%result = a * exp(-x.^2/(2*c^2));
%result = a * exp(-(x-tau/2).^2/(2*c^2));

function state = expand(theta, phi, j, binomial)
% Expand spin coherent state in spin eigenstates
% Return ket
state = zeros(2*j+1,1);
spineigs = eye(2*j+1);

%nchoosek(2*j, j+m)
for m = -j:j
    term = sqrt(binomial(index(j,m))) * cos(theta/2)^(j+m)... 
            * sin(theta/2)^(j-m) * exp(1i*(j-m)*phi)...
            * spineigs(:,index(j,m));
    state = state + term;
end

function result = stepUp(j)
% Step up operator
result = zeros(2*j+1,2*j+1);
for m = -j:j-1
    mm = index(j,m);
    result(mm,mm+1) = sqrt(j*(j+1) - m*(m+1));
end

function result = stepDown(j)
% Step down operator
result = zeros(2*j+1,2*j+1);
for m = -j+1:j
    mm = index(j,m);
    result(mm,mm-1) = sqrt(j*(j+1) - m*(m-1));
end

function result = Sx(j)
% S_x operator
result = 1/2 * ( stepUp(j) + stepDown(j) );

function result = Sy(j)
% S_y operator
result = -1i/2 * ( stepUp(j) - stepDown(j) );

function result = Sz(j)
% S_z operator
% CHANGED SIGN 6/12
result = -diag(-j:j);

function index = index(j, m)
% Convert quantum number m to vector index
index = m + j + 1;

function passed = rhoSanity(rho)
% Density matrix sanity check
passed = 1; % Assume all OK

if any(imag(diag(rho)) > 1e-12)
    disp('Warning! Rho has complex diagonal.')
    passed = 0;
end

if trace(rho)-1 > 1e-12
    disp('Warning! Rho is not normalized.')
    passed = 0;
end
    
if rho' ~= rho
    disp('Warning! Rho is not hermitian.')
    passed = 0;
end

if trace(rho^2)-1 > 1e-3
    disp('Warning! State significantly exceeds unity purity.')
    passed = 0;
end

function fsanity(t, tau)
plot(t,f(t,tau),'b-')
xlabel('t')
ylabel('f(t)')

function wienerCheck(t,dWstep)
if min(diff(t)) <= dWstep
    disp('Warning! ODE-solver step smaller than pseudo-Wiener increment, results may be compromised')
end




























