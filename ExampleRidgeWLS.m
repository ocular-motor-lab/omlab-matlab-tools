
 mber of observations
p = 2; % number of variables

MuPrior = zeros(p,1); % mean of the prior
SigmaPrior = eye(p); % covariance matrix of the prior
SigmaLikelihood = eye(n); % covariance matrix of the observations

X = randn(n,p); % design matrix

beta = chol(SigmaPrior) * randn(p,1)+ MuPrior; % draw betas from the prior
y = X*beta + chol(SigmaLikelihood)*randn(n,1); % generative model


% estimate beta and its varaince
betahat = (X'/SigmaLikelihood*X + eye(p)/SigmaPrior) \ (X'/SigmaLikelihood*y + SigmaPrior\MuPrior);
betahatVar = (X'/SigmaLikelihood*X + eye(p)/SigmaPrior) \ X'/SigmaLikelihood*X  /(X'/SigmaLikelihood*X + eye(p)/SigmaPrior)
