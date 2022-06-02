data
{

   real R;         // Number of recruits
   real S;         // Spawning stock biomass

   real kLow;      // supports of the parameters
   real kUp;
   real SkLow;
   real SkUp;
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

   real<lower=kLow,upper=kUp> k;
   real<lower=SkLow,upper=SkUp> Sk;
   real<lower=sigma2Low,upper=sigma2Up> sigma2;

}

transformed parameters
{

   real sigma;

   real mu;
   real ER;

   sigma = sqrt(sigma2);

   ER = k * S/Sk * exp(1-S/Sk);
   mu = log(ER) - 0.5*sigma2;

}

model
{

   k ~ uniform(kLow,kUp);
   Sk ~ uniform(SkLow,SkUp);
   sigma2 ~ uniform(sigma2Low,sigma2Up);

   R ~ lognormal(mu, sigma);

}

generated quantities
{

   real log_lik;

   log_lik = lognormal_log(R, mu, sigma);

}
