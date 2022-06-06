data
{

   int T;               // Number of time steps

   vector[T] R;         // Number of recruits
   vector[T] S;         // Spawning stock biomass

   real kLow;           // supports of the parameters
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

   vector[T] mu;
   vector[T] ER;

   sigma = sqrt(sigma2);

   for (t in 1:T) {
      ER[t] = k * S[t]/Sk * exp(1-S[t]/Sk);
   }
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

   vector[T] log_lik;

   for (t in 1:T) {
      log_lik[t] = lognormal_lpdf(R[t]| mu[t], sigma);
   }

}
