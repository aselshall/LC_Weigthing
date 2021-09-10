clear *; clc
%#ok<*SAGROW>

files={'R3210','R321X','R32XX','R3XXX'};
files={'RXXX0'};
for file=files
    clear err err1 arfitness MW
    
for run=1:3
%% [0] General simulation optimization opitions 
strfitnessfct='opt_a2Case04';                  % Simulation model function

%% [1] Set decision variables
folder='';
nm=6;
[LCA,KB,row,ZOSN,ZOSS,NR]=opt_ESM_AVISO_data(file{1},folder);
nd=length(NR);
PS=nd*100;
IP=[PS 200];                    % Population zie, Number of iterationsround(PS/10)

%% [2] CMA-ES setting and model selection
N =nd;                  % number of objective variables/problem diminsion
sigma = 0.5;            %(0.5) coordinate wise standard deviation (step size)
stopfitness = -1e9;     % stop if fitness < stopfitness (minimization)
lambda =IP(1);          % population size, offspring number 4+floor(3*log(N));
ITR=IP(2);              %Number of iterations
stopeval =ITR*lambda;   %stop after stopeval number of function evaluations 1e1*N^4;
rangemin=-4;                            %Normal distribution range
rangemax=4;                            %Normal distribution range
PrMsMin=0;
PrMsMax=1;
%disp(['Number of iterations ' num2str(ITR) '    Iteration size ' num2str(lambda)%])



%% [4] Initialize CMA-ES
%(4.1)  --------------------  Initialization --------------------------------
if contains(file,'R3210')
    xmean=[-3.78010101029446;-3.82631300294148;-0.366665990595674;-2.55419490924260;-4;
        -0.910087690538266;-1.13000654506987;4;-3.20792954394902;-2.49623053893893;-1.86535129564843];
elseif contains(file,'R321X')
    xmean=[-3.78010101029446;-3.82631300294148;-0.366665990595674;
        -0.910087690538266;-1.13000654506987;4;-3.20792954394902;-2.49623053893893;-1.86535129564843];
elseif contains(file,'R32XX')
    xmean=[-3.82631300294148;
        -0.910087690538266;-1.13000654506987;4;-3.20792954394902;-2.49623053893893;-1.86535129564843];
elseif contains(file,'R3XXX')
    xmean=[-1.13000654506987;4;-2.49623053893893;-1.86535129564843];

end
xmean=xmean+(rand(N,1)*1e-6);
SOF=7;
disp([file{1} '| Population size ' num2str(PS) ' | run ' num2str(run) ' | Decision variables ' num2str(nd) ' |  OF ' num2str(SOF)])

% Strategy parameter setting: Selection
mu = lambda/2;                          % number of parents/points for recombination
weights = log(mu+1/2)-log(1:mu)';       % muXone array for weighted recombination
mu = floor(mu);
weights = weights/sum(weights);         % normalize recombination weights array
mueff=sum(weights)^2/sum(weights.^2);   % variance-effectiveness of sum w_i x_i

% Strategy parameter setting: Adaptation
cc = (4 + mueff/N) / (N+4 + 2*mueff/N); % time constant for cumulation for C
cs = (mueff+2) / (N+mueff+5);           % t-const for cumulation for sigma control
c1 = 2 / ((N+1.3)^2+mueff);             % learning rate for rank-one update of C
cmu = 2 * (mueff-2+1/mueff) / ((N+2)^2+mueff);      % and for rank-mu update
damps = 1 + 2*max(0, sqrt((mueff-1)/(N+1))-1) + cs; % damping for sigma
% usually close to 1

% Initialize dynamic (internal) strategy parameters and constants
pc = zeros(N,1); ps = zeros(N,1);   % evolution paths for C and sigma
B = eye(N,N);                       % B defines the coordinate system
D = ones(N,1);                      % diagonal D defines the scaling
C = B * diag(D.^2) * B';            % covariance matrix C
invsqrtC = B * diag(D.^-1) * B';    % C^-1/2
eigeneval = 0;                      % track update of B and D
chiN=N^0.5*(1-1/(4*N)+1/(21*N^2));  % expectation of

%(4.2) -------------------- Generation Loop --------------------------------
%Generate solutions and obtain fitness (OF) of each solution
counteval = 0;  % the next 40 lines contain the 20 lines of interesting code
CurrentFitness=inf;

Start=0;
RErr=1e15;

%% [5] Start CMA-ES by generating and evaluationg solutions (in [1]Paralell or [2]Serial)
while counteval < stopeval   %&& CurrentFitness >= stopfitness
% Generate and evaluate lambda offspring
arx=zeros(N,lambda);            %Initialize solutions
arfitness=zeros(1,lambda);      %Initialize fintess (OF)

for k=1:lambda
    arx(:,k) = xmean + sigma * B * (D .* randn(N,1));   % m + sig * Normal(0,C)
    [err(k), err1(k),arfitness(k), MW(k,:)] = feval(strfitnessfct, arx(:,k),LCA,ZOSN,ZOSS,NR,KB,nm,row,rangemin,rangemax,PrMsMin,PrMsMax,SOF); % Objective function call 
    counteval = counteval+1;                            %Count solutions on by one
end

[Err,IDX]=min(err);
if Err<RErr
    RErr1=err1(IDX);
    RErr=Err;
    Rfitness=arfitness(IDX);
    RMW=MW(IDX,1:N);
end



%% [6] Analyize fitness values to propose new solutions
% Sort by fitness and compute weighted mean into xmean
[arfitness, arindex] = sort(arfitness);  % minimization
xold = xmean;
xmean = arx(:,arindex(1:mu)) * weights;  % recombination, new mean value
xmin = arx(:, arindex(1)); % Return best point of last iteration.
% Cumulation: Update evolution paths
ps = (1-cs) * ps ...
    + sqrt(cs*(2-cs)*mueff) * invsqrtC * (xmean-xold) / sigma;
hsig = sum(ps.^2)/(1-(1-cs)^(2*counteval/lambda))/N < 2 + 4/(N+1);
pc = (1-cc) * pc ...
    + hsig * sqrt(cc*(2-cc)*mueff) * (xmean-xold) / sigma;
% Adapt covariance matrix C
artmp = (1/sigma) * (arx(:,arindex(1:mu)) - repmat(xold,1,mu));  % mu difference vectors
C = (1-c1-cmu) * C ...                      % regard old matrix
    + c1 * (pc * pc' ...                    % plus rank one update
    + (1-hsig) * cc*(2-cc) * C) ...         % minor correction if hsig==0
    + cmu * artmp * diag(weights) * artmp'; % plus rank mu update
% Adapt step size sigma
sigma = sigma * exp((cs/damps)*(norm(ps)/chiN - 1));
%sigma=0.6;
% Update B and D from C
if counteval - eigeneval > lambda/(c1+cmu)/N/10 % to achieve O(N^2)
    eigeneval = counteval;
    C = triu(C) + triu(C,1)';                   % enforce symmetry
    [B,D] = eig(C);                             % eigen decomposition, B==normalized eigenvectors
    D = sqrt(diag(D));                          % D contains standard deviations now
    invsqrtC = B * diag(D.^-1) * B';
end
CurrentFitness=arfitness(1);

%% [7] Print to screen and save some outputs per iteration 

%(7.1) Print iteration progrss on screen 
more off;  % turn pagination off in Octave
ITR=counteval/lambda;
REM = mod(ITR,1);
if ITR==1 || REM==-0
disp([ 'Iteration#' num2str(ITR) ' :  Fitness ' num2str(arfitness(1)) ' : Min Err ' num2str(RErr) ' : Min Err1 ' num2str(RErr1)])
end

% if RErr1>2
%     break
% end

end


end
end