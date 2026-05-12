
set.seed(1234)
data <- ordmelslaplace::sim_ordinal_mels()



fit <- ordmelslaplace::ord_mels( fixed = y ~ x,
                 random = ~ 1 + x | id,
                 fixed_scale = ~ x,
                 data = data,
                 latent_error_dist = "probit",
                 nPoints = 5 )

# Mit startwerten
start <- c(-6.544256, 1.706173, 1.910515, 2.691738, 0.1763091, 1.020261,
           -0.2647777, -0.3425024, 0.9378529, -0.1859619, -0.3282065 )



fit <- ordmelslaplace::ord_mels( fixed = y ~ x,
                 random = ~ 1 + x | id,
                 fixed_scale = ~ x,
                 data = data,
                 latent_error_dist = "probit",
                 nPoints = 8,
                 start = start )

# fixed <- y ~ x
# random <- ~ 1 + x | id
# fixed_scale <- ~ x
#
# latent_error_dist <- "probit"
# nPoints <- 5
# SEs <- F
# starting_values <- NULL
# maxiter <- 100
# conv_tol <- 1e-4
# trace <- T
# nlminb_ctrl <- list( eval.max = 500, iter.max = 350 )
