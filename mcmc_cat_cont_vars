log_pdf_uniform <- function(x, lower, upper){

        return(sum(dunif(x,lower, upper, log=TRUE)))

}

​

log_pdf_normal <- function(x, mu, sig){

        return(sum(dnorm(x,mu,sig, log=TRUE)))

}

​

log_pdf_gamma <- function(x, alpha, beta){

        return(sum(dgamma(x,alpha,rate=beta, log=TRUE)))

}

​

log_pdf_exp <- function(x, l){

        # convert input into vector

        x = c(x)

        loglik = log(l)-l*x

        # sum up loglik values to get lik of whole data

        loglik_data = sum(loglik)

        return(loglik_data)

}

​

setwd("/Users/tobias/GitHub/bayesian_course_material/")

temp_tbl = read.table("data/global_temperature_NASA_scaled_time_start_1900.txt",h=T)

predictor_tbl = read.table("data/co2_ozone_data.txt",h=T)

temperature = temp_tbl$Temperature

CO2 = predictor_tbl$co2

O3 = predictor_tbl$ozone_depleting_substances

​

logfile <- "output/r_mcmc_samples_white_noise_multi_regression_indicators.txt"

cat(c("iteration","posterior","likelihood","prior","mu0","sig","a","b","sigma_hp","i_a","i_b","\n"),file=logfile,sep="\t")

​

#_____________PRIORS________________

# boundaries of the uniform prior on mu0

minMu = -10

maxMu = 40

# shape and rate parameters of the gamma prior on sigma

shapeSig = 2

rateSig = 1

# normal prior on a

muA = 0

# normal prior on b

muB = 0

# add the prior on-off switch rate, low values lead to fewer on-switch events, e.g. a rate of 0.1 will turn the switch on with 10% probability

rate_binomial = 0.5

# we use an exponential prior on the sigma hyperprior parameter

rate_hyp_sig = 1

#___________________________________

​

​

# init paprameters

current_mu0 = 14

current_sigma0 = 1

current_a = 0 # effect size CO2

current_b = 0 # effect size O3

current_sigma_hp = 1

current_ia = 0

current_ib = 0

​

mu_t = current_mu0 + current_a*current_ia*CO2 + current_b*current_ib*O3

current_likelihood = log_pdf_normal(temperature, mu_t, current_sigma0)

current_prior = log_pdf_uniform(current_mu0,minMu,maxMu) + log_pdf_gamma(current_sigma0,shapeSig,rateSig) + 

                log_pdf_normal(current_a, muA, current_sigma_hp) +

                log_pdf_normal(current_b, muB, current_sigma_hp) +

                log_pdf_exp(current_sigma_hp,rate_hyp_sig) + # <- hyper-prior!

                dbinom(current_ia,1,rate_binomial,log=TRUE) + # <- prior for on-off switch parameter

                dbinom(current_ib,1,rate_binomial,log=TRUE)

current_posterior = current_likelihood+current_prior # apply Bayes theorem to calculate the posterior

​

window_size_update_hyp_sig = 1

window_size_update_mu = 0.1

window_size_update_sig = 0.1

window_size_update_a = 0.05

window_size_update_b = 0.001

​

        

for (i in 1:500000){

        # reset new params

        new_mu0 =    current_mu0

        new_sigma0 = current_sigma0 

        new_a =      current_a 

        new_b =      current_b 

        new_sigma_hp = current_sigma_hp

        new_ia = current_ia

        new_ib = current_ib   

        

        # no hastings here becasue I'm using only symmetric proposals

        rand <- runif(1)

        if (rand < 0.05){

                new_ia = abs(current_ia - 1)

        } else if(rand < 0.1){

                new_ib = abs(current_ib - 1)

        } else if (rand < 0.4){

                new_mu0 = new_mu0 + rnorm(1,0,window_size_update_mu)

                new_sigma0 = abs(new_sigma0 + rnorm(1,0,window_size_update_sig))

        }else if (rand < 0.6){

                new_a = new_a + rnorm(1,0,window_size_update_a)

        }else if (rand < 0.8){

                new_b = new_b + rnorm(1,0,window_size_update_b)

        }else{

                new_sigma_hp = abs(new_sigma_hp + rnorm(1,0,window_size_update_hyp_sig))

        }

        hastings_ratio = 0 # in this case it's always 0 because we are using symmetric proposal functions

        

        mu_t = new_mu0 + new_a*new_ia*CO2 + new_b*new_ib*O3

        new_likelihood = log_pdf_normal(temperature, mu_t, new_sigma0)

        new_prior = log_pdf_uniform(new_mu0,minMu,maxMu) + log_pdf_gamma(new_sigma0,shapeSig,rateSig) + 

                log_pdf_normal(new_a, muA, new_sigma_hp) +

                log_pdf_normal(new_b, muB, new_sigma_hp) +

                log_pdf_exp(new_sigma_hp,rate_hyp_sig) + # <- hyper-prior!

                dbinom(new_ia,1,rate_binomial,log=TRUE) + # <- prior for on-off switch parameter

                dbinom(new_ib,1,rate_binomial,log=TRUE)

        new_posterior = new_likelihood+new_prior # apply Bayes theorem to calculate the new posterior

        

        post_ratio <- new_posterior-current_posterior

        # calculate acceptance probability as product of the posterior ratio and Hastings ratio (in log space!)

        r = post_ratio + hastings_ratio

        

        if (log(runif(1)) < r){

                current_posterior = new_posterior

                current_prior <- new_prior

                current_likelihood <- new_likelihood

                current_mu0 <- new_mu0  

                current_sigma0 <- new_sigma0

                current_a <- new_a    

                current_b <- new_b    

                current_sigma_hp <- new_sigma_hp

                current_ia = new_ia

                current_ib = new_ib                

        }

        if (i %% 100 == 0){

                cat(c(i,current_posterior, current_likelihood, current_prior, current_mu0, current_sigma0,

                         current_a, current_b, current_sigma_hp,current_ia,current_ib,"\n"),file=logfile,sep="\t",append=T)

        }

        if (i %% 5000 == 0){

                cat(c(i,current_posterior, current_likelihood, current_prior, current_mu0, current_sigma0,

                      current_a, current_b, current_sigma_hp,current_ia,current_ib,"\n"),sep="\t",append=T)

        }

}

​
