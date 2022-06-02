functions {
   real fair_prior_log(real c) {
      return log( 1 + 0.62269*( ( 1-2001^(2*c-2) ) / ( 1+2001^(2*c-2) ) ) );
   }
}

data
{

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

   sigma2Low = sigmaLow^2;
   sigma2Up = sigmaUp^2;

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
   sigma = sqrt(sigma2);

}

model
{

   Rinf ~ uniform(RinfLow, RinfUp);
   Sq ~ uniform(SqLow, SqUp);
   c ~ fair_prior();
   sigma2 ~ uniform(sigma2Low, sigma2Up);

}