function [PtfRisk, PtfReturn, MDD, OptimWeights, PI, ER] = StaticBLAlloc(Dates, returns, UsdMktCaps, EstStartDate, EstEndDate, lambda, P, Q, nbport, n, ExecDel, plotflag)
%%
tic
[IdxEstStartDate,~] = find(Dates == datenum(EstStartDate,'dd.mm.yyyy'));
[IdxEstEndDate,~] = find(Dates == datenum(EstEndDate,'dd.mm.yyyy'));

nbassets = size(returns,2);

%% STATIC OPTIMIZATION:  in sample (02.01.2004 to 29.09.2008) 

%% Michaud Resampled Efficiency with implied equilibrium returns incorporating investor's views (sim 3)

%% Bulding PI : the Market Equilibrium Implied Expected Return
% Compute the market cap weights for each asset and for each date
UsdMktCaps = UsdMktCaps(IdxEstStartDate:IdxEstEndDate,:);
TotMktCap = sum(UsdMktCaps,2);
MktCapWts = zeros(size(UsdMktCaps,1),size(UsdMktCaps,2));
for j = 1:size(UsdMktCaps,2)
    for i = 1:size(UsdMktCaps,1)
        MktCapWts (i,j) = UsdMktCaps(i,j) ./ TotMktCap(i);
    end
end

% Skewness of returns
alpha = skewness(returns(IdxEstStartDate:IdxEstEndDate,:));

%% Beware!... Black-Litterman model under construction...

% Building inputs

% Estimation of the variance covariance matrix 
cov_matrix = cov(returns(IdxEstStartDate:IdxEstEndDate,:)) + .0001 * eye(size(returns(IdxEstStartDate:IdxEstEndDate,:),2)); % in order to have it positive semi definite

nbviews = size(Q,1);

% scaling factor used in the Black-Litterman model
tau = 1/(IdxEstEndDate - IdxEstStartDate + 1);

% Omega : covariance matrix of the error's for the views
Omega = zeros(nbviews, nbviews);
for j = 1:nbviews
Omega(j,j) = (P(j,:) * cov_matrix * P(j,:)');
end

% Building PI and the Black-Litterman new combined return vector ER
PI = (lambda * cov_matrix * MktCapWts(IdxEstEndDate,:)')';
ER = (inv((inv(tau*cov_matrix)+P'*Omega*P))*inv(tau*cov_matrix)*PI'+P'*Omega*Q)';

% Simulation of 1000 scenarios of random returns drawn from a student distribution

df = 7; % kurto_sim1 = kurtosis(returns_sim1) : used to check that the simulated returns have approx the same distr than the historical ones
cases  = IdxEstEndDate - IdxEstStartDate + 1;
returns_sim3 = zeros(cases,nbassets,n);
mu3 = zeros(n,nbassets);
covmat3 = zeros(nbassets,nbassets,n);
PortRisk_sim3 = zeros(n,nbport);
PortReturn_sim3 = zeros(n,nbport);
PortWts_sim3 = zeros(nbport,nbassets,n);

parfor l = 1:n
    returns_sim3(:,:,l) = st_rndmvt( ER, cov_matrix, df, cases );
    mu3 = mean(returns_sim3(:,:,l));
    covmat3 = nancov(returns_sim3(:,:,l));
    [PortRisk_sim3(l,:), PortReturn_sim3(l,:), PortWts_sim3(:,:,l)] = portopt(mu3, covmat3, nbport);
end

% Average of the 1000 simulated optimal portfolio weights

StocksPortWts_RE3 = mean(PortWts_sim3,3);
OptimWeights = StocksPortWts_RE3;

%% Portfolios OOS Performance Comparison

% Out of sample returns

OOS_Stat1 = (OptimWeights * returns(IdxEstEndDate+ExecDel:end,:)')';
OOS_Stat1Eqw = ((1/nbassets) * ones(nbport,nbassets) * returns(IdxEstEndDate + ExecDel:end,:)')';

% NAVs
equity0 = 100;
NAV_Stat1 = nan(size(OOS_Stat1)+1);
NAV_Stat1Eqw = nan(size(OOS_Stat1Eqw)+1);

for j = 1:size(OOS_Stat1,2)
    
    for i = 2:size(OOS_Stat1,1)+1
        NAV_Stat1(1,j) = equity0;
        NAV_Stat1(i,j) = NAV_Stat1(i-1,j) * (1 + OOS_Stat1(i-1,j));
    end
end

for j = 1:size(OOS_Stat1Eqw,2)
    
    for i = 2:size(OOS_Stat1Eqw,1)+1
        NAV_Stat1Eqw(1,j) = equity0;
        NAV_Stat1Eqw(i,j) = NAV_Stat1Eqw(i-1,j) * (1 + OOS_Stat1Eqw(i-1,j));
    end
end

% Outputs

PtfRisk = sqrt(portvar(returns(IdxEstEndDate+ExecDel:end,:),vertcat(OptimWeights, (1/nbassets) * ones(1,nbassets)))) .* (sqrt(252) * ones(nbport+1,1));
PtfReturn = portror(mean(returns(IdxEstEndDate+ExecDel:end,:),1),vertcat(OptimWeights, (1/nbassets) * ones(1,nbassets))) .* (252*ones(1,nbport+1));
PtfReturn = PtfReturn';
MDD = maxdrawdown([NAV_Stat1(:,1:nbport) NAV_Stat1Eqw(:,1)]);

% Plot NAVs
if plotflag == 1
f = figure;
p1 = plot(Dates(IdxEstEndDate + ExecDel-1:end),NAV_Stat1); hold on; p2 = plot(Dates(IdxEstEndDate + ExecDel-1:end),NAV_Stat1Eqw); hold on; h = [p1(1:nbport);p2]; legend(h,{'Port1','Port2','Port3','Port4','Port5','Equally Weighted Port'},'Location','westoutside'); title('Black-Litterman MV Optimized Portfolios vs Equally Weighted Portfolio');
datetick('x','mmm-yy','keepticks','keeplimits')
end
toc
end
