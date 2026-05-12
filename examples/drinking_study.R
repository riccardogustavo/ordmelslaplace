


data <- ordmels::drinking_study


# for Testpurposes only use study2
#data <- data[ data$study == "study2", ]


str(data)
head(data)

data$age_z <- scale(data$age)
data$dealt.stress_z <- scale(data$dealt.stress)

#evtl nur study 1 oder 2 ?
fixed <- alc.intox ~ age_z + gender + dealt.stress_z
random <- ~ 1 + dealt.stress_z | PID
fixed_scale <- ~ age_z + gender + dealt.stress_z
latent_error_dist <- "logit"
nPoints <- 5

fit <- ordmels::ord_mels( fixed = fixed,
          random = random,
          fixed_scale = fixed_scale,
          latent_error_dist = latent_error_dist,
          data = data,
          nPoints = nPoints )

start <- c( -2.763517, 1.300459, -0.07035313, -0.06603048, 0.2877098, 0.8611179,
            0.0750252, 0.05442614, 0.00731079, 0.009606506, 0.1551004, -0.0009395244,
            -0.2166048, 0.004689187, 1.064253, -6.438353, -0.0647763,1 -3.768646 )


fit <- ordmelspureR::ord_mels( fixed = fixed,
                 random = random,
                 fixed_scale = fixed_scale,
                 latent_error_dist = latent_error_dist,
                 data = data,
                 nPoints = 5,
                 start = start )
