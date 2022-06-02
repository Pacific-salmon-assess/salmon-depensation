data
{

   real R;         // Number of recruits
   real S;         // Spawning stock biomass

   real alphaLow;           // supports of the parameters
   real alphaUp;
   real sigmaLow;
   real sigmaUp;

}

transformed data
{

   real sigma2Low;
   real sigma2Up;

   sigma2Low = sigmaLow^2;
   sigma2Up = sigmaUp^2;

}

parameters
{

   real<lower=alphaLow, upper=alphaUp> alpha;
   real<lower=sigma2Low, upper=sigma2Up> sigma2;


}

transformed parameters
{

   real sigma;

   real mu;
   real ER;

   sigma = sqrt(sigma2);

   ER = alpha * S;
   mu = log(ER) - 0.5*sigma2;

}

model
{

   alpha ~ uniform(alphaLow, alphaUp);
   sigma2 ~ uniform(sigma2Low, sigma2Up);

}

generated quantities
{

   real log_lik;

   log_lik = lognormal_log(R, mu, sigma);

}
