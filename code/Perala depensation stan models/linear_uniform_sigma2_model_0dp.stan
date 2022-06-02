data
{

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

   sigma = sqrt(sigma2);

}

model
{

   alpha ~ uniform(alphaLow, alphaUp);
   sigma2 ~ uniform(sigma2Low, sigma2Up);

}