data
{

   real R;         // Number of recruits
   real S;         // Spawning stock biomass

   real q;

   real RinfLow;           // supports of the parameters
   real RinfUp;
   real SqLow;
   real SqUp;
   real sigmaLow;
   real sigmaUp;

}

transformed data
{

   real sigma2Low;
   real sigma2Up;
   real p;

   sigma2Low = sigmaLow^2;
   sigma2Up = sigmaUp^2;
   p = (1-q)/q;

}

parameters
{

   real<lower=RinfLow,upper=RinfUp> Rinf;
   real<lower=SqLow,upper=SqUp> Sq;
   real<lower=sigma2Low,upper=sigma2Up> sigma2;

}

transformed parameters
{

   real sigma;

   real mu;
   real ER;

   sigma = sqrt(sigma2);

   ER = Rinf / (p * Sq / S + 1);
   mu = log(ER) - 0.5*sigma2;

}

model
{

   Rinf ~ uniform(RinfLow,RinfUp);
   Sq ~ uniform(SqLow,SqUp);
   sigma2 ~ uniform(sigma2Low,sigma2Up);

}

generated quantities
{

   real log_lik;

   log_lik = lognormal_log(R, mu, sigma);

}
