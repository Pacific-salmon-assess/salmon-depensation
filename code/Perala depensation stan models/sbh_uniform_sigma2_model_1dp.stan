functions {
   real fair_prior_log(real c) {
      return log( 1 + 0.62269*( ( 1-2001^(2*c-2) ) / ( 1+2001^(2*c-2) ) ) );
   }
}

data
{

   real R;         // Number of recruits
   real S;         // Spawning stock biomass

   real q;

   real RinfLow;           // supports of the parameters
   real RinfUp;
   real SqLow;
   real SqUp;
   real cLow;
   real cUp;
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

   real<lower=RinfLow, upper=RinfUp> Rinf;
   real<lower=SqLow, upper=SqUp> Sq;
   real<lower=cLow, upper=cUp> c;
   real<lower=sigma2Low, upper=sigma2Up> sigma2;

}

transformed parameters
{

   real sigma;

   real mu;
   real ER;

   sigma = sqrt(sigma2);

   ER = Rinf / (p * pow(Sq / S, c) + 1);
   mu = log(ER) - 0.5*sigma2;

}

model
{

   Rinf ~ uniform(RinfLow, RinfUp);
   Sq ~ uniform(SqLow, SqUp);
   c ~ fair_prior();
   sigma2 ~ uniform(sigma2Low, sigma2Up);

   R ~ lognormal(mu, sigma);

}

generated quantities
{

   real log_lik;

   log_lik = lognormal_log(R, mu, sigma);

}
