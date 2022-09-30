source('package_code/1_load_all.R')
library( tictoc )

binom_innov <- function(n, x) {
  fam <- binomial()
  rbinom(n, prob = fam$linkinv(x), size = 1)
}

conj_tab <- function( x, y ) {
  
  fits <- x$glm_diagnostics$fitted.values >= 0.5
  table( fits, y )
  
}

b <- c( -2, 0.5, 0.5, 0, 0 )
nvar <- length(b) - 1
set.seed(123)
test <- irts_midas_sim( n_y = 1500, n_num = nvar, seasonal_adjust = TRUE,
                        beta = b, innov_fun = binom_innov,
                        DLF_parameters = list( c( 14, -1 ), c( 0.5, -0.1 ),
                                               c( 14, -1 ), c( 0.5, -0.1 ) ) )$Data

mrtest <- midas_regimes(value = test[1:nvar], time = test$INDEX, DLF_parameters = list(c(0, -0.1)))
rtest <- response_regime(value = test$y, time = test$INDEX)
form <- rtest ~ mrtest

tic()
fit_non_gr <- irts_midas( form, data = test, family = "binomial", 
                          method = "BFGS", gr = NULL )
toc() 

# 2 covars: ~19 seconds
# 4 covars: ~65 seconds

non_gr_tab <- conj_tab( fit_non_gr, rtest$All$value )
sum( diag(non_gr_tab) ) / sum( non_gr_tab )

tic()
fit <- irts_midas( form, data = test, family = "binomial", method = "BFGS" )
toc() 

# 2 covars: ~9 seconds
# 4 covars: ~8 seconds

gr_tab <- conj_tab( fit, rtest$All$value )
sum( diag(gr_tab) ) / sum( gr_tab )


