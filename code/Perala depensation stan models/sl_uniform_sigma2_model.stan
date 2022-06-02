functions {
   real fair_prior_log(real c) {
      return log( 1 + 0.62269*( ( 1-2001^(2*c-2) ) / ( 1+2001^(2*c-2) ) ) );
   }
}

data
{

   int T;               // Number of time steps

   vector[T] R;         // Number of recruits
   vector[T] S;         // Spawning stock biomass

   real kLow;           // supports of the parameters
   real kUp;
   real SkLow;
   real SkUp;
   real cLow;
   real cUp;
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
   real<lower=cLow,upper=cUp> c;
   real<lower=sigma2Low,upper=sigma2Up> sigma2;

}

transformed parameters
{

   real sigma;

   vector[T] mu;
   vector[T] ER;

   sigma = sqrt(sigma2);

   for (t in 1:T) {
      ER[t] = k * pow(S[t]/Sk,c) * exp(c*(1-S[t]/Sk));
   }
   mu = log(ER) - 0.5*sigma2;

}

model
{

   k ~ uniform(kLow,kUp);
   Sk ~ uniform(SkLow,SkUp);
   c ~ fair_prior();
   sigma2 ~ uniform(sigma2Low,sigma2Up);

   R ~ lognormal(mu, sigma);

}

generated quantities
{

   vector[T] log_lik;

   for (t in 1:T) {
      log_lik[t] = lognormal_log(R[t], mu[t], sigma);
   }

}
