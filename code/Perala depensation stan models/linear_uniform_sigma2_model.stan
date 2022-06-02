data
{

   int T;               // Number of time steps

   vector[T] R;         // Number of recruits
   vector[T] S;         // Spawning stock biomass

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

   vector[T] mu;
   vector[T] ER;

   sigma = sqrt(sigma2);

   for (t in 1:T) {
      ER[t] = alpha * S[t];
   }

   mu = log(ER) - 0.5*sigma2;

}

model
{

   alpha ~ uniform(alphaLow, alphaUp);
   sigma2 ~ uniform(sigma2Low, sigma2Up);

   R ~ lognormal(mu, sigma);

}

generated quantities
{

   vector[T] log_lik;

   for (t in 1:T) {
      log_lik[t] = lognormal_log(R[t], mu[t], sigma);
   }

}
