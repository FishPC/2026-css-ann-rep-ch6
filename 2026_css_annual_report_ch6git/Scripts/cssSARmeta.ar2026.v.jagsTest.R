# setup -------------------------------------------------------------------
# 1---source functions
source(
  'Functions/packFontHandler.R'
  )

# 2---call/install packages and reconcile fonts
## function can be modified @ Scripts/Functions/packHandler.R
packFontHandler()

# data steps --------------------------------------------------------------
## Not run:
# 1---create 'notin' operator
# '%notin%' <- Negate(
#   '%in%'
#   )
# End(**Not run**)

# 2---define vector for wild fish
target.wild <- c(
  'ROSA', # Yakima River wild Chinook
  'JDAC', # John Day River wild Chinook
  'YAKS', # Yakima River wild steelhead
  'JDAS', # John Day River wild steelhead
  'AGCW', # Snake River wild spring-summer Chinook
  'AGWS', # Snake River wild steelhead aggregate
  'EMCR', # Entiat and Methow rivers wild Chinook
  'EMWS'  # Entiat and Methow rivers wild steelhead
)

## Not run:
# 3---define vector for wild and hatchery stocks
# target.tot <- c(
#   'AGGR',
#   'AGGA',
#   'AGGB'
# )
# End(**Not run**)

# 4---set target grouping (should be one of the above: 'target.wild' or 'target.tot')
targ.grp <- 'target.wild'

# 5---read-in and manipulate data
sar.meta.dat <- readRDS(
  'Data/sarMeta.rds'
) %>% 
  filter(
    css.grp
    %in%
      !! sym(targ.grp)
  ) %>%
  filter(
    mig.yr >= 1994
  ) %>% 
  ## add fields for adult returns, row id and distance to CRM
  mutate(
    cases = (sar.est/100)*juv.pop,
    obs.id = 1:n(),
    dist_rkm = ifelse(
      css.grp == 'EMCR' | css.grp == 'EMWS',
      529,
      ifelse(
        css.grp == 'AGCW' | css.grp == 'AGWS',
        461,
        ifelse(
          css.grp == 'ROSA' | css.grp == 'YAKS',
          236,
          ifelse(
            css.grp == 'JDAC' | css.grp == 'JDAS',
            113,
            NA
          )
        )
      )
    ),
    orig_grp = ifelse(
      css.grp == 'EMCR' | css.grp == 'EMWS',
      'E_M',
      ifelse(
        css.grp == 'AGCW' | css.grp == 'AGWS',
        'SNK',
        ifelse(
          css.grp == 'ROSA' | css.grp == 'YAKS',
          'YAK',
          ifelse(
            css.grp == 'JDAC' | css.grp == 'JDAS',
            'JDA',
            NA
          )
        )
      )
    )
  ) %>%  
  group_by(
    css.grp
  ) %>% 
  mutate(
    es.id = seq_along(css.grp)
  )

# 6---calculate summary effect sizes
sar.meta.es <- escalc (
  xi = cases, # adult returns
  ni = juv.pop, # juvenile population
  data = sar.meta.dat , # specify data (from above)
  measure = 'PLO',  # data transformation (logit in this case)
  add = 1/100,  # add constant to account for '0' cell entries
  to = 'only0' # specify which records to apply constant to
)
  drop_na(yi)




#######################JAGS INTERPRETATION######################################
dat <- list(
  N = nrow(sar.meta.es),
  yi = sar.meta.es$yi,
  vi = sar.meta.es$vi,
  
  orig_grp = as.integer(factor(sar.meta.es$orig_grp)),
  obs_id   = as.integer(factor(sar.meta.es$obs.id)),
  spp_code = as.integer(factor(sar.meta.es$spp.code)),
  mig_yr   = as.integer(factor(sar.meta.es$mig.yr)),
  
  N_orig_grp = length(unique(sar.meta.es$orig_grp)),
  N_obs_id   = length(unique(sar.meta.es$obs.id)),
  N_spp      = length(unique(sar.meta.es$spp.code)),
  N_mig_yr   = length(unique(sar.meta.es$mig.yr))
)

inits <- function() list(
  beta0 = rnorm(1, 0, 1),
  sigma_orig_grp = runif(1, 0.1, 1),
  sigma_obs_id   = runif(1, 0.1, 1),
  sigma_spp      = runif(1, 0.1, 1),
  sigma_mig_yr   = runif(1, 0.1, 1),
  
  u_orig_grp = rnorm(dat$N_orig_grp, 0, 0.1),
  u_obs_id   = rnorm(dat$N_obs_id, 0, 0.1),
  u_spp      = rnorm(dat$N_spp, 0, 0.1),
  u_mig_yr   = rnorm(dat$N_mig_yr, 0, 0.1)
)

params <- c(
  "beta0",
  "sigma_orig_grp",
  "sigma_obs_id",
  "sigma_spp",
  "sigma_mig_yr",
  "u_orig_grp", 
  "u_obs_id", 
  "u_spp",
  "u_mig_yr"
)

model_string_ml_like <- "
model {

  for (i in 1:N) {

    #### Latent sampling variance (always defined)
    vi_latent[i] ~ dlnorm(0, 1.0E-6)   # extremely weak prior, always positive

    #### Link observed vi to latent vi
    vi[i] ~ dnorm(vi_latent[i], 1.0E6)  # ignored when vi[i] = NA

    #### Precision uses latent vi
    tau_obs[i] <- 1 / vi_latent[i]

    #### Latent true effect
    y_latent[i] ~ dnorm(mu[i], tau_obs[i])

    #### Link observed yi to latent yi
    yi[i] ~ dnorm(y_latent[i], 1.0E6)   # ignored when yi[i] = NA

    #### Linear predictor
    mu[i] <- beta0 +
             u_orig_grp[orig_grp[i]] +
             u_obs_id[obs_id[i]] +
             u_spp[spp_code[i]] +
             u_mig_yr[mig_yr[i]]
  }

  #### Fixed intercept (flat ≈ ML)
  beta0 ~ dnorm(0, 1.0E-12)

  #### Random effects: orig_grp
  for (j in 1:N_orig_grp) {
    u_orig_grp[j] ~ dnorm(0, tau_orig_grp)
  }
  tau_orig_grp <- 1 / (sigma_orig_grp * sigma_orig_grp)
  sigma_orig_grp ~ dunif(0, 1000)

  #### Random effects: obs.id
  for (j in 1:N_obs_id) {
    u_obs_id[j] ~ dnorm(0, tau_obs_id)
  }
  tau_obs_id <- 1 / (sigma_obs_id * sigma_obs_id)
  sigma_obs_id ~ dunif(0, 1000)

  #### Random effects: spp.code
  for (j in 1:N_spp) {
    u_spp[j] ~ dnorm(0, tau_spp)
  }
  tau_spp <- 1 / (sigma_spp * sigma_spp)
  sigma_spp ~ dunif(0, 1000)

  #### Random effects: mig.yr
  for (j in 1:N_mig_yr) {
    u_mig_yr[j] ~ dnorm(0, tau_mig_yr)
  }
  tau_mig_yr <- 1 / (sigma_mig_yr * sigma_mig_yr)
  sigma_mig_yr ~ dunif(0, 1000)

}

"
writeLines(model_string_ml_like, "meta_model_ml_like.jags")



library(rjags)

m <- jags.model(
  file = "meta_model.jags",
  data = dat,
  inits = inits,
  n.chains = 3,
  n.adapt = 2000
)

update(m, 5000)   # burn-in

samples <- coda.samples(
  model = m,
  variable.names = params,
  n.iter = 20000,
  thin = 10
)

summary(samples)
plot(samples)

posterior_means <- summary(samples)$statistics[, "Mean"]
posterior_means

################################################################################



# 7---generate data summary table
table.6.1 <- sar.meta.es %>% 
  dplyr::select(
    zone,
    # css.grp,
    sar.reach,
    spp.code
  ) %>%  
  unique() %>%  
  mutate(
    zone=replace(zone, zone=='SNAK', 'Snake R.'),
    zone=replace(zone, zone=='MCOL', 'Mid. Col. R'),
    zone=replace(zone, zone=='UCOL', 'Upp. Col. R.'),
    zone=replace(zone, zone=='UCOL', 'Upp. Col. R.'),
    spp.code = replace(spp.code, spp.code=='CH', 'Chinook'),
    spp.code = replace(spp.code, spp.code=='ST', 'steelhead')
  ) %>%  
  mutate(
    riv = c(
      rep('Snake R.',2),
      rep('Yakima R.',2),
      rep('John Day R.',2),
      rep('Entiat and Methow R.', 2)
    )
  ) %>%  
  mutate(
    n = c(
      nrow(filter(sar.meta.es, css.grp=='AGCW')),
      nrow(filter(sar.meta.es, css.grp=='AGWS')),
      nrow(filter(sar.meta.es, css.grp=='ROSA')),
      nrow(filter(sar.meta.es, css.grp=='YAKS')),
      nrow(filter(sar.meta.es, css.grp=='JDAC')),
      nrow(filter(sar.meta.es, css.grp=='JDAS')),
      nrow(filter(sar.meta.es, css.grp=='EMCR')),
      nrow(filter(sar.meta.es, css.grp=='EMWS'))
    )
  ) %>%
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    zone  = 'Zone',
    # css.grp = 'CSS group ID',
    spp.code = 'Species',
    riv = 'Origin',
    sar.reach = 'SAR est. reach',
    n = 'No. estimates'
  ) %>% 
  set_formatter(
    n = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f', big.mark = ','))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>%
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  width(
    j = c(
      1,
      3,
      5
    ), 
    1.1
  ) %>%  
  width(
    j = c(
      2
    ),
    3.0
  ) %>% 
  width(
    j = c(
      4
    ),
    1.7
  ) %>%  
  mk_par( 
    j = 2,
    i = 1:2,
    value = as_paragraph(
      'Lwr. Granite',
      as_sub('juv.'),
      '\u2013',
      'Bonn.',
      as_sub('ad.')
    )
  ) %>% 
  mk_par( 
    j = 2,
    i = 3:4,
    value = as_paragraph(
      'McNary',
      as_sub('juv.'),
      '\u2013',
      'Bonn.',
      as_sub('ad.')
    )
  ) %>%  
  mk_par( 
    j = 2,
    i = 5:6,
    value = as_paragraph(
      'John Day',
      as_sub('juv.'),
      '\u2013',
      'Bonn.',
      as_sub('ad.')
    )
  ) %>% 
  mk_par( 
    j = 2,
    i = 7:8,
    value = as_paragraph(
      'Rocky Reach',
      as_sub('juv.'),
      '\u2013',
      'Bonn.',
      as_sub('ad.')
    )
  ) %>%
  print()

# save_as_image(
#   'Output\\Tables\\Table 6.1.png',
#   res = 2000
# )

# specify models ----------------------------------------------------------
# 1---assess basic model structural components
# a---fit models
# i---model 0 -- basic multi-level model (i.e., null model):
## intercept = TRUE;
## random effects = NA; 
## fixed effects = NA;
## var-cov matrix = unstructured
mod0 <- rma.mv(
  yi = yi, 
  V = vi,
  slab = css.grp,
  data = sar.meta.es,
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

# ii---model 1:
## intercept = TRUE; 
## random effect(s) = group-level; 
## fixed effects = NA;
## var-cov matrix = unstructured
mod1 <- rma.mv(
  yi = yi, 
  V = vi,
  slab = orig_grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(orig_grp)
  ), 
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

# iii---model 2:
## intercept = TRUE;
## random effects = observation-level within css.grp; 
## fixed effects = NA;
## var-cov matrix = unstructured
mod2 <- rma.mv(
  yi = yi, 
  V = vi,
  slab = orig_grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(orig_grp)/factor(obs.id)
  ), 
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

# iv---model 3:
## intercept = TRUE;
## random effects = observation-level within css.grp; migration year; 
## fixed effects = NA;
## var-cov matrix = unstructured
mod3 <- rma.mv(
  yi = yi, 
  V = vi,
  slab = css.grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(orig_grp)/factor(obs.id),
    ~ 1 | mig.yr
  ),
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

# v---model 4:
## intercept = TRUE;
## random effects = observation-level within css.grp; mig. yr. within spp. 
## fixed effects = NA;
## var-cov matrix = unstructured
mod4 <- rma.mv(
  yi = yi, 
  V = vi,
  slab = css.grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(orig_grp)/factor(obs.id),
    ~ 1 | factor(spp.code)/mig.yr
  ),
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

# b---create table of overall (structural) model comparisons
# i---manipulate data
bsFitStats <- fitstats(
  mod0,
  mod1,
  mod2,
  mod3,
  mod4
) %>% 
  t() %>% 
  as.data.frame()

table.6.3 <- bsFitStats %>% 
  dplyr::select('logLik:',
                'AICc:'
  ) %>%  
  rename(
    log.likliehood = 1,
    aicc = 2
  ) %>% 
  mutate(
    'Intercept' = c('+', '+','+','+','+'),
    'Fixed.Effects' = c('-','-','-','-','-'),
    'Random.Effects' = c('-','1|origin.grp.','1|origin.grp/obs.id','1|origin.grp/obs.id; 1|mig. yr.', '1|origin.grp/obs.id; 1|spp./mig.yr')
  ) %>% 
  arrange(
    aicc
  ) %>%  
  mutate(
    delta.AICc = round(
      akaike.weights(aicc)$deltaAIC,
      0
    ),
    AICc.weights = round(
      akaike.weights(aicc)$weights,
      2
    ),
    no.param = arrange(
      AIC.rma(
        mod0,
        mod1,
        mod2,
        mod3,
        mod4
      ),
      AIC
    )$df
  ) %>% 
  relocate(
    c(
      delta.AICc, 
      AICc.weights
    ),
    .after = aicc
  ) %>%  
  relocate(
    c(
      Intercept,
      Fixed.Effects,
      Random.Effects,
      no.param
      
    ),
    .before = log.likliehood
  ) %>% 
  # ii---generate output table
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    Intercept = 'Intercept',
    Fixed.Effects = 'Fixed Effects',
    Random.Effects = 'Random Effects',
    no.param = 'No. params.',
    log.likliehood = 'log-liklihood'
  ) %>% 
  mk_par( 
    j = 6,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 7,
    i = 1,
    value = as_paragraph(
      paste0('\u394','AIC'),
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j =8,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c'),
      ' wt.'
    ),
    part = 'header'
  ) %>%
  set_formatter(
    log.likliehood = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f', big.mark = ',')),
    aicc = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f', big.mark = ',')),
    delta.AICc = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f', big.mark = ',')),
    AICc.weights = function(x) ifelse(is.na(x),'', formatC(x,digits = 2, format = 'f', big.mark = ','))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  width(
    j = 2, 
    1.1
  ) |> 
  width(
    j = 3, 
    2.4
  ) |> 
  width(
    j = 4, 
    1.0
  ) %>% 
  print()

# save_as_image(
#   'Output\\Tables\\Table 6.3.png',
#   res = 2000
# )

# hypothesis testing ------------------------------------------------------
#1---multi-model inference
## helper functions necessary to coordinate metafor and MuMin
# eval(metafor:::.MuMIn)

#a---fit full model
## update best supported model from above, with moderators
top.mod <- bsFitStats %>% 
  arrange(
    `AICc:`
  ) %>% 
  head(1) %>% 
  rownames_to_column() %>%  
  dplyr::select(
    rowname
  ) %>%
  as.character()

fullMod <- update(
  eval(
    parse(
      text = top.mod
    )
  ),
  mods = ~ factor(spp.code) + wtt + pitph
)

## Not run:
#b---estimate variance inflation factors
#i---calculate VIF in closed-form and simulate props < vif
# table.3 <- vif.rma(
#   fullMod,
#   sim = TRUE,
#   parallel = 'snow',
#   ncpus = detectCores(),
#   seed = 1234
# ) %>%
#   #ii---manipulated data
#   as.data.frame() %>%
#   mutate(
#     mods = c(
#       'WTT',
#       'PITPH'
#     )
#   ) %>%
#   dplyr::select(
#     mods,
#     vif,
#     prop
#   ) %>%
#   #iii---generate table
#   flextable() %>%
#   align(
#     i = 1,
#     j = 1,
#     align = 'center',
#     part =  'header'
#   ) %>%
#   align(
#     align = 'center',
#     part = 'all'
#   ) %>%
#   set_header_labels(
#     mods = 'Moderator',
#     vif = 'VIF'
#   ) %>%
#   mk_par(
#     j =3,
#     i = 1,
#     value = as_paragraph(
#       'Prop.',
#       as_sub('sims.'),
#       ' < VIF'
#     ),
#     part = 'header'
#   ) %>%
#   set_formatter(
#     mods = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
#     vif = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
#     prop = function(x) ifelse(is.na(x),'', formatC(x,digits = 2, format = 'f'))
#   ) %>%
#   fontsize(
#     size = 12,
#     part = 'all'
#   ) %>%
#   flextable::font(
#     fontname = 'Times New Roman',
#     part = 'all'
#   ) %>%
#   border_remove() %>%
#   hline(
#     i=1,
#     part='header',
#     border = fp_border(
#       color='black',
#       width = 1
#     )
#   ) %>%
#   hline_top(
#     part='header',
#     border = fp_border(
#       color='black',
#       width = 2
#     )
#   ) %>%
#   hline_bottom(
#     part='body',
#     border = fp_border(
#       color='black',
#       width = 2
#     )
#   ) %>%
#   autofit() %>%
#   print()

#iv---output table
# save_as_image(
#   'Output/Tables/sarMetaTbl3.png',
#   res = 2000
# )
# End(**Not run**)

# c---conduct MMI and create and output summary table
table.6.4 <- dredge(
  fullMod,
  trace=2
) %>%
  as.data.frame() %>%
  mutate(
    `factor(spp.code)` = c(
      NA,
      update(
        fullMod,
        ~ factor(spp.code) + wtt + pitph
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      update(
        fullMod,
        ~ factor(spp.code) + wtt
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      NA,
      NA,
      update(
        fullMod,
        ~ factor(spp.code) + pitph
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      update(
        fullMod,
        ~ factor(spp.code)
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      NA
    )
  ) %>% 
  relocate(`factor(spp.code)`, .before = pitph) %>%
  # i---generate output table
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    `(Intercept)` = 'Intercept',
    `factor(spp.code)` = 'Spp.',
    pitph = 'PITPH',
    wtt = 'WTT',
    df = 'No. params.',
    logLik = 'log-liklihood'
  ) %>%
  mk_par( 
    j = 7,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 8,
    i = 1,
    value = as_paragraph(
      paste0('\u394','AIC'),
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 9,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c'),
      ' wt.'
    ),
    part = 'header'
  ) %>%
  set_formatter(
    `factor(spp.code)` = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    pitph = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    wtt = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    df = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f')),
    AICc = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f')),
    delta = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f')),
    weight = function(x) ifelse(is.na(x),'', formatC(x,digits = 2, format = 'f'))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>%
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%
  width(
    j = 4, 
    0.95
  ) %>%
  width(
    j = 5, 
    1.05
  ) %>%
  print()

# save_as_image(
#   'Output\\Tables\\Table 6.4.png',
#   res = 2000
# )

# d---specify best model
# i---update full model based on model selection table
bestMod <- update(
  eval(
    parse(
      text = top.mod
    )
  ),
  mods = ~ wtt + pitph
)

# e---plot model-averaged and best model coefficients and CIs
# i---estimate model-averaged coefficients
avgMod <- dredge(
  fullMod,
  trace = 3
) %>% 
  subset(
    delta <= 4.0 # specify confidence set to include models where dAICc <= 4.0
    # recalc.weights = FALSE
  ) %>% 
  model.avg()

#ii---data manipulation
figure.6.3 <- bind_rows(
  data.frame(
    selType = as.factor(rep('best',3)),
    param.lbs = c(
      'Spp.',
      'WTT',
      'PITPH'
    ),
    param = c(Spp. = NA,bestMod$beta[2:3,]),
    lcl = c(
      NA,
      confint(
        bestMod, 
        fixed = TRUE,
        level = 0.95 # specify 95% confidence limits
      ) %>% 
        as.data.frame() %>% 
        as.data.frame() %>% 
        slice(-1) %>% 
        dplyr::select(
          ci.lb
        ) %>% 
        unlist() %>% 
        unname()
    ),
    ucl = c(
      NA,
      confint(
        bestMod, 
        fixed = TRUE,
        level = 0.95 # specify 90% confidence limits
      ) %>% 
        as.data.frame() %>% 
        as.data.frame() %>% 
        slice(-1) %>% 
        dplyr::select(
          ci.ub
        ) %>% 
        unlist() %>% 
        unname()
    )
  ) %>% 
    rename(
      Estimate = 3
    ) %>% 
    `rownames<-`( NULL ),
  coefTable(
    avgMod, 
    full = TRUE
  ) %>%
    as.data.frame() %>% 
    dplyr::select(
      Estimate
    ) %>% 
    tibble::rownames_to_column('param.lbs') %>% 
    left_join(
      confint(
        avgMod, 
        full = TRUE,
        level = 0.95 # specify 90% confidence limits
      ) %>% 
        as.data.frame() %>% 
        tibble::rownames_to_column('param.lbs'),
      by = 'param.lbs'
    ) %>%  
    as.data.frame() %>%  
    rename(
      lcl = 3,
      ucl = 4
    ) %>%
    mutate(
      param.lbs = c(
        'Intercept',
        'PITPH',
        'WTT',
        'Spp.'
      )
    ) %>% 
    tail(4) %>% 
    mutate(
      selType = as.factor(rep('average',4))
    ) %>% 
    relocate(
      selType,
      .before = param.lbs
    )
) %>% 
  as.data.frame() %>% 
  filter(
    param.lbs != 'Intercept'
  ) %>% 
  #iii---generate plot
  ggplot(aes(x=Estimate, y = param.lbs, color = selType))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 18,vjust = 1,color = 'black',family = 'Calibri'),
        axis.title.x = element_text(face = 'bold', size = 18,vjust = -1,color = 'black',family = 'Calibri'),
        axis.text.x = element_text(face = 'bold',size = 16,angle = 0,hjust = 0.5,vjust = 0.5,color = 'black',family = 'Calibri'),
        axis.text.y = element_text(face = 'bold',size = 16,color = 'black',family = 'Calibri'),
        axis.ticks.length = unit(2,'mm'),
        legend.position = 'right',
        legend.title = element_blank(),
        legend.text = element_text(size = 16,family = 'Calibri',face = 'bold'))+
  scale_colour_manual(name='Transplanted', values = c('best' = '#D55E00','average' = 'black'), breaks = c('best', 'average')) +
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_errorbar(aes(xmax=ucl,xmin=lcl), height = 0, position = position_dodge(width = -0.3), linewidth = 1.0, orientation = 'y') +
  geom_point(size = 3.5, position = position_dodge(width = -0.3)) +
  # geom_point(size = 4,shape = 1, position = position_dodge(width = 0.3)) +
  labs(y = 'Parameter',x = 'Estimate') +
  scale_y_discrete(limits = c('WTT','PITPH','Spp.')) + 
  scale_x_continuous(limits=c(-0.6,0.6),
                     breaks = seq(-0.6,0.6,0.3),
                     labels = function(x) format(x, scientific = FALSE)
  )
# geom_magnify(
#   from = c(xmin = -0.25, xmax = 0.25, ymin = 0.9, ymax = 2.1),
#   to = c(1.2 - 0.5, 1.2 + 0.5,1.5, 2.5)
# )

# iv---print figure for review
print(figure.6.3)

#v--- output figure
# png(filename='Output\\Figures\\Figure 6.3.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.3)
# dev.off()
# print(figure.6.3)

# assumption assessment ---------------------------------------------------
## group-specific interactions
### define base model
baseMod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph
)

### define model with stock x WTT interaction
wtt_grp_inter.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + pitph + wtt:factor(orig_grp)
)

### define model with stock x PITPH interaction
pitph_grp_inter.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph:factor(orig_grp)
)

## species specific interactions
### define model with spp x WTT interaction
wtt_spp_inter.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph + wtt:factor(spp.code)
)

### define model with spp x PITPH interaction
pitph_spp_inter.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph + pitph:factor(spp.code)
)

## species specific interactions
### define model with WTT x PITPH interaction
add_assump.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph + wtt:pitph
)

## species specific interactions
### define model with distance fixed effect
dist_assump.mod <- update(
  bestMod,
  mods = ~ factor(spp.code) + wtt + pitph + dist_rkm
)

## AR assumption
### define model with distance fixed effect
redREmodAR <- rma.mv(
  yi = yi,
  V = vi,
  slab = css.grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(orig_grp)/factor(obs.id),
    ~ mig.yr|factor(spp.code)
  ),
  struct = 'AR',
  mods = ~ factor(spp.code) + wtt + pitph,
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

table.6.2 <- data.frame(
  mod = c(
    'logit(SAR) ~ Spp. + WTT + PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + WTT:Origin grp.',
    'logit(SAR) ~ Spp. + WTT + PITPH + PITPH:Origin grp.',
    'logit(SAR) ~ Spp. + WTT + PITPH + Spp.:WTT',
    'logit(SAR) ~ Spp. + WTT + PITPH + Spp.:PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + WTT:PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + Distance',
    'logit(SAR) ~ Spp. + WTT + PITPH (AR-1 on nested mig. yr. RE)' 
  ),
  assump = c(
    'base model',
    'origin-specific slopes (Origin x WTT)',
    'origin-specific slopes (Origin x PITPH)',
    'species-specific slopes (Spp. x WTT)',
    'species-specific slopes (Spp. x PITPH)',
    'additivity',
    'distance',
    'AR(1) on nested mig. yr. random effect'
  )
) %>% 
  flextable() %>%  
  align(
    i = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  # align(
  #   j = c(2,3),
  #   align = 'center',
  #   part = 'body'
  # ) %>% 
  set_header_labels(
    mod = 'Model Form',
    assump = 'Assumption'
  ) %>% 
  # set_formatter(
  #   # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
  #   bic = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f'))
  # ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  width(
    j = c(
      1
    ), 
    4.6
  ) %>%
  width(
    j = c(
      2
    ), 
    3.0
  ) %>% 
  print()

# save_as_image(
#   'Output\\Tables\\Table 6.2.png',
#   res = 2000
# )

Table.6.5 <- data.frame(
  mod = c(
    'logit(SAR) ~ Spp. + WTT + PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + WTT:Origin grp.',
    'logit(SAR) ~ Spp. + WTT + PITPH + PITPH:Origin grp.',
    'logit(SAR) ~ Spp. + WTT + PITPH + Spp.:WTT',
    'logit(SAR) ~ Spp. + WTT + PITPH + Spp.:PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + WTT:PITPH',
    'logit(SAR) ~ Spp. + WTT + PITPH + Distance',
    'logit(SAR) ~ Spp. + WTT + PITPH (AR-1 on  nested mig. yr. RE)'
  ),
  assump = c(
    'base model',
    'origin-specific slopes (Origin x WTT)',
    'origin-specific slopes (Origin x PITPH)',
    'species-specific slopes (Spp. x WTT)',
    'species-specific slopes (Spp. x PITPH)',
    'additivity',
    'distance',
    'AR(1) on nested mig. yr. random effect'
  ),
  bic = c(
    BIC.rma(baseMod),
    BIC.rma(wtt_grp_inter.mod),
    BIC.rma(pitph_grp_inter.mod),
    BIC.rma(wtt_spp_inter.mod),
    BIC.rma(pitph_spp_inter.mod),
    BIC.rma(add_assump.mod),
    BIC.rma(dist_assump.mod),
    BIC.rma(redREmodAR)
  )
) %>% 
  arrange(
    bic
  )%>% 
  flextable() %>%  
  align(
    i = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    j = c(3),
    align = 'center',
    part = 'body'
  ) %>% 
  set_header_labels(
    mod = 'Model Form',
    assump = 'Assumption',
    bic = 'BIC'
  ) %>% 
  set_formatter(
    # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    bic = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f'))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  width(
    j = c(
      1
    ), 
    4.7
  ) %>%
  width(
    j = c(
      2
    ), 
    2.8
  ) %>% 
  print()

# save_as_image(
#   'Output\\Tables\\Table 6.5.png',
#   res = 2000
# )

# variance components -----------------------------------------------------
#c--- variance component tables (i.e., Table 6.6)
#i--- full model (all moderators)
fullVC.mod <- update(
  redREmodAR,
  mods = ~ wtt + pitph
)

fullVC.mod2 <- update(
  redREmodAR,
  mods = ~ factor(spp.code) + wtt + pitph
)

fullVC.tbl <- data.frame(
  sigNm = c(
    '\\sigma^2_1',
    '\\sigma^2_2',
    '\\tau^2'
  ),
  sigEst = c(paste0(round(fullVC.mod$sigma2,4),'/',round(fullVC.mod2$sigma2,4)),paste0(round(fullVC.mod$tau2,4),'/',round(fullVC.mod2$tau2,4))),
  sQRtsigEst = c(paste0(round(sqrt(fullVC.mod$sigma2),4),'/',round(sqrt(fullVC.mod2$sigma2),4)),paste0(round(sqrt(fullVC.mod$tau2),4),'/',round(sqrt(fullVC.mod2$tau2),4))),
  nLvls = c(fullVC.mod2$s.nlevels.f,paste0(fullVC.mod2$g.nlevels.f[1],'|',fullVC.mod2$g.nlevels.f[2]))
) %>% 
  mutate(
    fac = c(
      'Origin grp.',
      'Origin grp./obs.',
      'mig. yr.|spp. code'
    )
  ) %>% 
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    sigNm = 'Parameter',
    sigEst = 'Estimate',
    sQRtsigEst = 'Sqr. Rt. Estimate',
    nLvls = 'No. Levels',
    fac = 'Factor'
  ) %>% 
  set_formatter(
    # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    # sQRtsigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    nLvls = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f'))
  ) %>% 
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>%
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%  
  mk_par(
    j = 'sigNm',
    value = as_paragraph(as_equation(sigNm))
  ) %>% 
  autofit() %>%
  set_caption(
    caption = as_paragraph(
      as_chunk(
        "Full models from MCS", 
        props = fp_text_default(
          bold = TRUE, 
          font.family = 'Times New Roman', 
          font.size = 14
        )
      )
    ),
    align_with_table = FALSE,
    word_stylename = "Table Caption",
    fp_p = fp_par(text.align = "left", padding = 3)
  ) %>% 
  print()

# save_as_docx(
#   path = 'Output/Tables/Variance Components/fullVC_6.6.docx'
# )

#ii--- reduced model (no moderators)
redVC_nm.mod <- update(
  redREmodAR,
  ~ 1
)

redVC_nm.mod2 <- update(
  redREmodAR,
  ~ factor(spp.code)
)

redVC_nm.tbl <- data.frame(
  sigNm = c(
    '\\sigma^2_1',
    '\\sigma^2_2',
    '\\tau^2'
  ),
  sigEst = c(paste0(round(redVC_nm.mod$sigma2,4),'/',round(redVC_nm.mod2$sigma2,4)),paste0(round(redVC_nm.mod$tau2,4),'/',round(redVC_nm.mod2$tau2,4))),
  sQRtsigEst = c(paste0(round(sqrt(redVC_nm.mod$sigma2),4),'/',round(sqrt(redVC_nm.mod2$sigma2),4)),paste0(round(sqrt(redVC_nm.mod$tau2),4),'/',round(sqrt(redVC_nm.mod2$tau2),4))),
  nLvls = c(redVC_nm.mod2$s.nlevels.f,paste0(redVC_nm.mod2$g.nlevels.f[1],'|',redVC_nm.mod2$g.nlevels.f[2]))
) %>% 
  mutate(
    fac = c(
      'Origin grp.',
      'Origin grp./obs.',
      'mig. yr.|spp. code'
    )
  ) %>% 
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    sigNm = 'Parameter',
    sigEst = 'Estimate',
    sQRtsigEst = 'Sqr. Rt. Estimate',
    nLvls = 'No. Levels',
    fac = 'Factor'
  ) %>% 
  set_formatter(
    # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    # sQRtsigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    nLvls = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f'))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%  
  mk_par(
    j = 'sigNm',
    value = as_paragraph(as_equation(sigNm))
  ) %>%
  autofit() %>% 
  set_caption(
    caption = as_paragraph(
      as_chunk(
        "Reduced model (mods.:MCS1 = none/MCS2 = Spp.)", 
        props = fp_text_default(
          bold = TRUE, 
          font.family = 'Times New Roman', 
          font.size = 14
        )
      )
    ),
    align_with_table = FALSE,
    word_stylename = "Table Caption",
    fp_p = fp_par(text.align = "left", padding = 3)
  ) %>% 
  print()

# save_as_docx(
#   path = 'Output/Tables/Variance Components/redVC_nm_6.6.docx'
# )

#ii--- reduced model (PITPH only)
redVC_pitph.mod <- update(
  redREmodAR,
  ~ pitph
)

redVC_pitph.mod2 <- update(
  redREmodAR,
  ~ factor(spp.code) + pitph
)

redVC_pitph.tbl <- data.frame(
  sigNm = c(
    '\\sigma^2_1',
    '\\sigma^2_2',
    '\\tau^2'
  ),
  sigEst = c(paste0(round(redVC_pitph.mod$sigma2,4),'/',round(redVC_pitph.mod2$sigma2,4)),paste0(round(redVC_pitph.mod$tau2,4),'/',round(redVC_pitph.mod2$tau2,4))),
  sQRtsigEst = c(paste0(round(sqrt(redVC_pitph.mod$sigma2),4),'/',round(sqrt(redVC_pitph.mod2$sigma2),4)),paste0(round(sqrt(redVC_pitph.mod$tau2),4),'/',round(sqrt(redVC_pitph.mod2$tau2),4))),
  nLvls = c(redVC_pitph.mod2$s.nlevels.f,paste0(redVC_pitph.mod2$g.nlevels.f[1],'|',redVC_pitph.mod2$g.nlevels.f[2]))
) %>% 
  mutate(
    fac = c(
      'Origin grp.',
      'Origin grp./obs.',
      'mig. yr.|spp. code'
    )
  ) %>% 
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    sigNm = 'Parameter',
    sigEst = 'Estimate',
    sQRtsigEst = 'Sqr. Rt. Estimate',
    nLvls = 'No. Levels',
    fac = 'Factor'
  ) %>% 
  set_formatter(
    # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    # sQRtsigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    nLvls = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f'))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%  
  mk_par(
    j = 'sigNm',
    value = as_paragraph(as_equation(sigNm))
  ) %>% 
  autofit() %>%
  set_caption(
    caption = as_paragraph(
      as_chunk(
        'Reduced model (mods.:MCS1 = PITPH/MCS2 = Spp. + PITPH)', 
        props = fp_text_default(
          bold = TRUE, 
          font.family = 'Times New Roman', 
          font.size = 14
        )
      )
    ),
    align_with_table = FALSE,
    word_stylename = "Table Caption",
    fp_p = fp_par(text.align = "left", padding = 3)
  ) %>%
  print()

# save_as_docx(
#   path = 'Output/Tables/Variance Components/redVC_pitph_6.6.docx'
# )

#ii--- reduced model (WTT only)
redVC_wtt.mod <- update(
  redREmodAR,
  ~ wtt
)

redVC_wtt.mod2 <- update(
  redREmodAR,
  ~ factor(spp.code) + wtt
)

redVC_wtt.tbl <- data.frame(
  sigNm = c(
    '\\sigma^2_1',
    '\\sigma^2_2',
    '\\tau^2'
  ),
  sigEst = c(paste0(round(redVC_wtt.mod$sigma2,4),'/',round(redVC_wtt.mod2$sigma2,4)),paste0(round(redVC_wtt.mod$tau2,4),'/',round(redVC_wtt.mod2$tau2,4))),
  sQRtsigEst = c(paste0(round(sqrt(redVC_wtt.mod$sigma2),4),'/',round(sqrt(redVC_wtt.mod2$sigma2),4)),paste0(round(sqrt(redVC_wtt.mod$tau2),4),'/',round(sqrt(redVC_wtt.mod2$tau2),4))),
  nLvls = c(redVC_wtt.mod2$s.nlevels.f,paste0(redVC_wtt.mod2$g.nlevels.f[1],'|',redVC_wtt.mod2$g.nlevels.f[2]))
) %>%  
  mutate(
    fac = c(
      'Origin grp.',
      'Origin grp./obs.',
      'mig. yr.|spp. code'
    )
  ) %>% 
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    sigNm = 'Parameter',
    sigEst = 'Estimate',
    sQRtsigEst = 'Sqr. Rt. Estimate',
    nLvls = 'No. Levels',
    fac = 'Factor'
  ) %>% 
  set_formatter(
    # sigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    # sQRtsigEst = function(x) ifelse(is.na(x),'', formatC(x,digits = 4, format = 'f')),
    nLvls = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f'))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>%  
  mk_par(
    j = 'sigNm',
    value = as_paragraph(as_equation(sigNm))
  ) %>% 
  autofit() %>% 
  set_caption(
    caption = as_paragraph(
      as_chunk(
        'Reduced model (mods.:MCS1 = WTT/MCS2 = Spp. + WTT)',
        props = fp_text_default(
          bold = TRUE, 
          font.family = 'Times New Roman', 
          font.size = 14
        )
      )
    ),
    align_with_table = FALSE,
    word_stylename = "Table Caption",
    fp_p = fp_par(text.align = "left", padding = 3)
  ) %>%
  print()

# save_as_docx(
#   path = 'Output/Tables/Variance Components/redVC_wtt_6.6.docx'
# )

# c--- estimate pseudo-R^2
# i--- proportion reduction in origin-level RE
spp_or_r2 <- 1-((redVC_nm.mod$sigma2[1] - redVC_nm.mod2$sigma2[1]) / redVC_nm.mod$sigma2[1])
pitph_or_r2 <- 1-((redVC_nm.mod$sigma2[1] - redVC_pitph.mod2$sigma2[1]) / redVC_nm.mod$sigma2[1])
wtt_or_r2 <- 1-((redVC_nm.mod$sigma2[1] - redVC_wtt.mod2$sigma2)[1] / redVC_nm.mod$sigma2[1])
full_or_r2 <- 1-((redVC_nm.mod$sigma2[1] - fullVC.mod2$sigma2[1]) / redVC_nm.mod$sigma2[1])

# ii--- proportion reduction in observation-level RE
spp_obs_r2 <- 1-((redVC_nm.mod$sigma2[2] - redVC_nm.mod2$sigma2[2]) / redVC_nm.mod$sigma2[2])
pitph_obs_r2 <- 1-((redVC_nm.mod$sigma2[2] - redVC_pitph.mod2$sigma2[2]) / redVC_nm.mod$sigma2[2])
wtt_obs_r2 <- 1-((redVC_nm.mod$sigma2[2] - redVC_wtt.mod2$sigma2)[2] / redVC_nm.mod$sigma2[2])
full_obs_r2 <- 1-((redVC_nm.mod$sigma2[2] - fullVC.mod2$sigma2[2]) / redVC_nm.mod$sigma2[2])

# iii--- proportion reduction in total variance
spp_r2 <- 1-((sum(redVC_nm.mod$sigma2[1]) - sum(redVC_nm.mod2$sigma2)[1]) / sum(redVC_nm.mod$sigma2)[1])
pitph_r2 <- 1-((sum(redVC_nm.mod$sigma2[1]) - sum(redVC_pitph.mod2$sigma2)[1]) / sum(redVC_nm.mod$sigma2)[1])
wtt_r2 <- 1-((sum(redVC_nm.mod$sigma2[1]) - sum(redVC_wtt.mod2$sigma2)[1]) / sum(redVC_nm.mod$sigma2)[1])
full_r2 <- 1-((sum(redVC_nm.mod$sigma2[1]) - sum(fullVC.mod2$sigma2)[1]) / sum(redVC_nm.mod$sigma2)[1])

# d--- plot pseudo-R^2
# i--- set fill parameters
mod.fill <- c(
  "Spp.+WTT+PITPH"="black",
  "Spp."="#56B4E9",
  'Spp.+PITPH' = '#D55E00',
  'Spp.+WTT' = '#009E73'
)

# i--- individual REs
r2_indRE.plot <- data.frame(
  mod = c(
    'Spp.',
    'Spp.+PITPH',
    'Spp.+WTT',
    'Spp.+WTT+PITPH',
    'Spp.',
    'Spp.+PITPH',
    'Spp.+WTT',
    'Spp.+WTT+PITPH'
  ),
  fac = c(
    'origin',
    'origin',
    'origin',
    'origin',
    'obs.',
    'obs.',
    'obs.',
    'obs.'
  ),
  est = c(
    pmin(as.numeric(spp_or_r2),1),
    pmin(as.numeric(pitph_or_r2),1),
    pmin(as.numeric(wtt_or_r2),1),
    pmin(as.numeric(full_or_r2),1),
    pmin(as.numeric(spp_obs_r2),1),
    pmin(as.numeric(pitph_obs_r2),1),
    pmin(as.numeric(wtt_obs_r2),1),
    pmin(as.numeric(full_obs_r2),1)
  )
) %>%  
  mutate(
    index = if_else(
      fac == 'origin',
      1,
      if_else(
        fac == 'obs.',
        2,
        NA
      )
    )
  ) %>% 
  mutate(
    mod = factor(
      mod,
      levels = c(
        "Spp.",
        "Spp.+PITPH",
        "Spp.+WTT",
        'Spp.+WTT+PITPH'
      ),
      ordered = TRUE
    )
  ) %>%  
  mutate(
    mod_index = if_else(
      mod == 'Spp.',
      1,
      if_else(
        mod == 'Spp.+PITPH',
        2,
        if_else(
          mod == 'Spp.+WTT',
          3,
          if_else(
            mod == 'Spp.+WTT+PITPH',
            4,
            NA
          )
        )
      )
    )
  ) %>% 
  ggplot(aes(x = factor(index), y = est, fill = mod)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 26,margin = margin(t = 0, r = 10, b = 0, l = 0), family = "Calibri"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 22,colour = "black", family = "Calibri"),
        axis.text.y = element_text(size = 22,colour = "black", family = "Calibri"),
        legend.text = element_text(size = 16, family = "Calibri"),
        axis.ticks.length=unit(2.0, "mm"),
        # legend.position = 'none',
        legend.position = c(0.88,0.90),
        legend.title = element_text(size = 16, family = "Calibri", face = 'bold'),
        legend.key.width = unit(1.0,"cm"),
        plot.margin = unit(c(0.15,0.20,0.10,0.10), "in"))+
  labs(x = "Component", y = bquote("1-R"["pseudo"]^"2")) +
  theme(plot.title = element_text(hjust = 0.5,size = 28)) +
  geom_bar(stat="identity", color="black", position="dodge")+
  scale_fill_manual(name="Moderators",values = mod.fill)+
  scale_y_continuous(limits=c(0,1.0),breaks = seq(0,1.0,0.1),expand=c(0,0)) +
  scale_x_discrete(labels=c('Origin','Obs.'))

# iv---print figure for review
print(r2_indRE.plot)

#v--- output figure
# png(filename='Output\\Figures\\r2_indREPlot.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(r2_indRE.plot)
# dev.off()
# print(r2_indRE.plot)

# i--- overall variance
r2_totVar.plot <- data.frame(
  mod = c(
    'Spp.',
    'Spp.+PITPH',
    'Spp.+WTT',
    'Spp.+WTT+PITPH'
  ),
  est = c(
    pmin(as.numeric(spp_r2),1),
    pmin(as.numeric(pitph_r2),1),
    pmin(as.numeric(wtt_r2),1),
    pmin(as.numeric(full_r2),1)
  )
) %>%  
  mutate(
    mod = factor(
      mod,
      levels = c(
        "Spp.",
        "Spp.+PITPH",
        "Spp.+WTT",
        'Spp.+WTT+PITPH'
      ),
      ordered = TRUE
    )
  ) %>%  
  mutate(
    mod_index = if_else(
      mod == 'Spp.',
      1,
      if_else(
        mod == 'Spp.+PITPH',
        2,
        if_else(
          mod == 'Spp.+WTT',
          3,
          if_else(
            mod == 'Spp.+WTT+PITPH',
            4,
            NA
          )
        )
      )
    )
  ) %>% 
  ggplot(aes(x = factor(mod_index), y = est, fill = mod)) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "black"),
        axis.title.y = element_text(size = 26,margin = margin(t = 0, r = 10, b = 0, l = 0), family = "Calibri"),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 22,colour = "black", family = "Calibri"),
        axis.text.y = element_text(size = 22,colour = "black", family = "Calibri"),
        legend.text = element_text(size = 16, family = "Calibri"),
        axis.ticks.length=unit(2.0, "mm"),
        legend.position = 'none',
        # legend.position = c(0.82,0.88),
        legend.title = element_text(size = 16, family = "Calibri"),
        legend.key.width = unit(1.5,"cm"),
        plot.margin = unit(c(0.15,0.10,0.10,0.10), "in"))+
  labs(x = "Model", y = bquote("1-R"["pseudo"]^"2")) +
  theme(plot.title = element_text(hjust = 0.5,size = 28)) +
  geom_bar(stat="identity", color="black", position="dodge")+
  scale_fill_manual(name="Moderators",values = mod.fill)+
  scale_y_continuous(limits=c(0,1.0),breaks = seq(0,1.0,0.1),expand=c(0,0)) +
  scale_x_discrete(labels=c("Spp.",'Spp.+PITPH','Spp.+WTT','Spp.+WTT+PITPH'))

print(r2_totVar.plot)

# iv---print figure for review
print(r2_totVar.plot)

#v--- output figure
# png(filename='Output\\Figures\\r2_totVarPlot.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(r2_totVar.plot)
# dev.off()
# print(r2_totVar.plot)

# final model simplification ---------------------------------------------------
## given no little group effect in variance components, examine model with
## reduced error structure

# a--- original model with nested obs.id
origREmod <- update(
  redREmodAR,
  mods = ~ factor(spp.code) + wtt + pitph
)

# a--- reduced model with obs.id
redREmod <- rma.mv(
  yi = yi,
  V = vi,
  slab = css.grp,
  data = sar.meta.es,
  random = list(
    ~ 1 | factor(obs.id),
    ~ mig.yr|factor(spp.code)
  ),
  struct = 'AR',
  mods = ~ factor(spp.code) + wtt + pitph,
  test = 't',
  dfs = 'contain',
  method = 'ML'
)

bsFitStats <- fitstats(
  origREmod,
  redREmod
) %>% 
  t() %>% 
  as.data.frame()

table.6.7 <- bsFitStats %>% 
  dplyr::select('logLik:',
                'BIC:',
                'AICc:'
  ) %>%  
  rename(
    log.likliehood = 1,
    bic = 2,
    aicc = 3
  ) %>% 
  mutate(
    'Intercept' = c('+', '+'),
    'Fixed.Effects' = c('Spp.; WTT; PITPH','Spp.; WTT; PITPH'),
    'Random.Effects' = c('1|Origin Grp./Obs.ID; Mig.Yr.|Spp.','1|Obs.ID; Mig.Yr.|Spp.')
  ) %>% 
  arrange(
    aicc
  ) %>%  
  mutate(
    delta.AICc = round(
      akaike.weights(aicc)$deltaAIC,
      0
    ),
    AICc.weights = round(
      akaike.weights(aicc)$weights,
      2
    ),
    no.param = arrange(
      AIC.rma(
        origREmod,
        redREmod
      ),
      AIC
    )$df
  ) %>% 
  relocate(
    c(
      delta.AICc, 
      AICc.weights
    ),
    .after = aicc
  ) %>%  
  relocate(
    c(
      Intercept,
      Fixed.Effects,
      Random.Effects,
      no.param
      
    ),
    .before = log.likliehood
  ) %>% 
  # ii---generate output table
  flextable() %>%  
  align(
    i = 1,
    j = 1,
    align = 'center',
    part =  'header'
  ) %>% 
  align(
    align = 'center',
    part = 'all'
  ) %>% 
  set_header_labels(
    Intercept = 'Intercept',
    Fixed.Effects = 'Fixed Effects',
    Random.Effects = 'Random Effects',
    no.param = 'No. params.',
    log.likliehood = 'log-liklihood',
    bic = 'BIC'
  ) %>% 
  mk_par( 
    j = 7,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 8,
    i = 1,
    value = as_paragraph(
      paste0('\u394','AIC'),
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j =9,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c'),
      ' wt.'
    ),
    part = 'header'
  ) %>%
  set_formatter(
    log.likliehood = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f', big.mark = ',')),
    bic = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f', big.mark = ',')),
    aicc = function(x) ifelse(is.na(x),'', formatC(x,digits = 1, format = 'f', big.mark = ',')),
    delta.AICc = function(x) ifelse(is.na(x),'', formatC(x,digits = 0, format = 'f', big.mark = ',')),
    AICc.weights = function(x) ifelse(is.na(x),'', formatC(x,digits = 2, format = 'f', big.mark = ','))
  ) %>%  
  fontsize(
    size = 12,
    part = 'all'
  ) %>% 
  flextable::font(
    fontname = 'Times New Roman',
    part = 'all'
  ) %>% 
  border_remove() %>% 
  hline(
    i=1,
    part='header',
    border = fp_border(
      color='black',
      width = 1
    )
  ) %>% 
  hline_top(
    part='header',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  hline_bottom(
    part='body',
    border = fp_border(
      color='black',
      width = 2
    )
  ) %>% 
  width(
    j = 2, 
    1.5
  ) |> 
  width(
    j = 3, 
    2.6
  ) |> 
  width(
    j = 4, 
    1.0
  ) %>%
  width(
    j = 5, 
    1.1
  ) %>% 
  print()

# save_as_docx(
#   path = 'Output\\Tables\\Table 6.7.png'
# )

# model verification (average predictions) --------------------------------
# a---specify component models
# i---first component model
predMod1 <- update(
  redREmod,
  ~ wtt + pitph
)

# ii---second component model
predMod2 <- update(
  redREmod,
  ~ factor(spp.code) + pitph + wtt
)

# b---estimate BLUPs (best linear unbiased predictions for random effects)
# i---first component model
randEffEsts_mod1 <- ranef(
  predMod1
)

# ii---second component model
randEffEsts_mod2 <- ranef(
  predMod2
)

# b---generate weighted predictions
# i---first component model
# *---generate predictions data set
sarMeta.pred_mod1 <- predict(
  predMod1,
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = sar.meta.es %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = sar.meta.es %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = sar.meta.es %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = sar.meta.es %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = sar.meta.es %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    nestGrpObs = paste0(
      sar.meta.es %>%
        dplyr::select(
          obs.id,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          obs.id
        ) %>% 
        unlist() %>% 
        unname()
    ),
    nestSppMY = paste0(
      sar.meta.es %>%
        dplyr::select(
          mig.yr,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          mig.yr
        ) %>% 
        unlist() %>% 
        unname(),
      ' | ',
      sar.meta.es %>%
        dplyr::select(
          spp.code,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          spp.code
        ) %>% 
        unlist() %>% 
        unname()
    )
  ) %>% 
  left_join(
    tibble::rownames_to_column(randEffEsts_mod1$`~mig.yr | factor(spp.code`[1], 'nestSppMY') %>%  
      as.data.frame(),
    by = 'nestSppMY'
  ) %>% 
  rename(
    fnestSppMY = intrcpt
  ) %>%
  left_join(
    tibble::rownames_to_column(randEffEsts_mod1$`factor(obs.id)`[1], 'nestGrpObs') %>%  
      as.data.frame(),
    by = 'nestGrpObs'
  ) %>%  
  rename(
    fnestGrpObs = intrcpt
  ) %>% 
  # left_join(
  #   tibble::rownames_to_column(randEffEsts_mod1$orig_grp[1], 'orig_grp') %>%  
  #     as.data.frame(),
  #   by = 'orig_grp'
  # ) %>% 
  # rename(
  #   forig_grp = intrcpt
  # ) %>% 
  # left_join(
  #   tibble::rownames_to_column(randEffEsts_mod1$'factor(spp.code)'[1], 'spp.code') %>%  
  #     as.data.frame(),
  #   by = 'spp.code'
  # ) %>% 
  # rename(
  #   fspp.code = intrcpt
  # ) %>% 
  rowwise() %>% 
  mutate(
    flogitPred = rowSums(
      pick(
        pred, 
        fnestSppMY, 
        fnestGrpObs,
        # fspp.code
      ), 
      na.rm = FALSE
    )
  ) %>%
  mutate(
    arithPred = inv.logit(pred),
    farithPred = inv.logit(flogitPred)
  ) %>% 
  mutate(
    zone=replace(zone, zone=='SNAK', 'Snake R.')
  ) %>% 
  mutate(
    zone=replace(zone, zone=='MCOL', 'Middle Columbia R.')
  ) %>%
  mutate(
    zone=replace(zone, zone=='UCOL', 'Upper Columbia R.')
  ) %>%
  mutate(
    across(
      zone,
      factor
    )
  ) %>% 
  as.data.frame()

# *---generate model verification data set
modVerif_mod1 <- sarMeta.pred_mod1 %>%  
  mutate(
    arithObs = na.omit(sar.meta.es$sar.est)/100,
    spp = sar.meta.es %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(na.omit(sar.meta.es$sar.est)/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      avgMod$msTable[5] %>% 
        head(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_farithPred = farithPred * weight,
    wt_logitPred = pred * weight,
    wt_flogitPred = flogitPred * weight
  )

# ii---second component model
# *---generate predictions data set
sarMeta.pred_mod2 <- predict(
  predMod2,
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = sar.meta.es %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = sar.meta.es %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = sar.meta.es %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = sar.meta.es %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = sar.meta.es %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    nestGrpObs = paste0(
      sar.meta.es %>%
        dplyr::select(
          obs.id,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          obs.id
        ) %>% 
        unlist() %>% 
        unname()
    ),
    nestSppMY = paste0(
      sar.meta.es %>%
        dplyr::select(
          mig.yr,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          mig.yr
        ) %>% 
        unlist() %>% 
        unname(),
      ' | ',
      sar.meta.es %>%
        dplyr::select(
          spp.code,
          cases
        ) %>% 
        na.omit() %>% 
        dplyr::select(
          spp.code
        ) %>% 
        unlist() %>% 
        unname()
    )
  ) %>% 
  left_join(
    tibble::rownames_to_column(randEffEsts_mod2$`~mig.yr | factor(spp.code`[1], 'nestSppMY') %>%  
      as.data.frame(),
    by = 'nestSppMY'
  ) %>% 
  rename(
    fnestSppMY = intrcpt
  ) %>%
  left_join(
    tibble::rownames_to_column(randEffEsts_mod2$`factor(obs.id)`[1], 'nestGrpObs') %>%  
      as.data.frame(),
    by = 'nestGrpObs'
  ) %>%  
  rename(
    fnestGrpObs = intrcpt
  ) %>% 
  # left_join(
  #   tibble::rownames_to_column(randEffEsts_mod2$orig_grp[1], 'orig_grp') %>%  
  #     as.data.frame(),
  #   by = 'orig_grp'
  # ) %>% 
  # rename(
  #   forig_grp = intrcpt
  # ) %>% 
  # left_join(
  #   tibble::rownames_to_column(randEffEsts_mod2$'factor(spp.code)'[1], 'spp.code') %>%  
  #     as.data.frame(),
  #   by = 'spp.code'
  # ) %>% 
  # rename(
  #   fspp.code = intrcpt
  # ) %>% 
  rowwise() %>% 
  mutate(
    flogitPred = rowSums(
      pick(
        pred, 
        fnestSppMY, 
        fnestGrpObs
        # fspp.code
      ), 
      na.rm = FALSE
    )
  ) %>%
  mutate(
    arithPred = inv.logit(pred),
    farithPred = inv.logit(flogitPred)
  ) %>% 
  mutate(
    zone=replace(zone, zone=='SNAK', 'Snake R.')
  ) %>% 
  mutate(
    zone=replace(zone, zone=='MCOL', 'Middle Columbia R.')
  ) %>%
  mutate(
    zone=replace(zone, zone=='UCOL', 'Upper Columbia R.')
  ) %>%
  mutate(
    across(
      zone,
      factor
    )
  ) %>% 
  as.data.frame()

# *---generate model verification data set
modVerif_mod2 <- sarMeta.pred_mod2 %>%  
  mutate(
    arithObs = na.omit(sar.meta.es$sar.est)/100,
    spp = sar.meta.es %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(na.omit(sar.meta.es$sar.est)/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      avgMod$msTable[5] %>% 
        head(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_farithPred = farithPred * weight,
    wt_logitPred = pred * weight,
    wt_flogitPred = flogitPred * weight
  )

# iii---combine predictions from component models
modVerif_avg <- bind_rows(
  modVerif_mod1 %>% 
    dplyr::select(
      wt_arithPred,
      wt_farithPred,
      wt_logitPred,
      wt_flogitPred
    ) %>%  
    rownames_to_column(), 
  modVerif_mod2 %>%
    dplyr::select(
      wt_arithPred,
      wt_farithPred,
      wt_logitPred,
      wt_flogitPred
    ) %>% 
    rownames_to_column()
) %>% 
  group_by(rowname) %>%
  summarise_all(sum) %>% 
  arrange(
    as.numeric(rowname)
  ) %>% 
  dplyr::select(
    -rowname
  ) %>% 
  mutate(
    zone = modVerif_mod1$zone,
    spp = modVerif_mod1$spp,
    arithObs = modVerif_mod1$arithObs,
    logitObs = modVerif_mod1$logitObs
  ) %>% 
  rename(
    arithPred = wt_arithPred,
    farithPred = wt_farithPred,
    logitPred = wt_logitPred,
    flogitPred = wt_flogitPred
  )

# c---model performance
# i---estimate index of agreement
# *---marginal predictions; arithmetic scale
sarMeta.IOAarithFix <- modStats(
  modVerif_avg,
  mod = 'arithPred',
  obs = 'arithObs',
  statistic = c('IOA')
) %>%
  as.data.frame() %>% 
  dplyr::select(
    IOA
  ) %>% 
  round(2) %>% 
  as.numeric()

# *---marginal predictions; logit scale
sarMeta.IOAlogitFix <- modStats(
  modVerif_avg,
  mod = 'logitPred',
  obs = 'logitObs',
  statistic = c('IOA')
) %>%
  as.data.frame() %>% 
  dplyr::select(
    IOA
  ) %>% 
  round(2) %>% 
  as.numeric()

# *---conditional predictions; arithmetic scale
sarMeta.IOAarithRand <- modStats(
  modVerif_avg,
  mod = 'farithPred',
  obs = 'arithObs',
  statistic = c('IOA')
) %>%
  as.data.frame() %>% 
  dplyr::select(
    IOA
  ) %>% 
  round(2) %>% 
  as.numeric()

# *---conditional predictions; logit scale
sarMeta.IOAlogitRand <- modStats(
  modVerif_avg,
  mod = 'flogitPred',
  obs = 'logitObs',
  statistic = c('IOA')
) %>%
  as.data.frame() %>% 
  dplyr::select(
    IOA
  ) %>% 
  round(2) %>% 
  as.numeric()

# ii---generate plots
# *---marginal predictions; arithmetic scale
figure.6.4_a <- ggplot(
  data = modVerif_avg %>% 
    mutate(
      sppCol = paste0(spp,', ',zone)
    ),
  aes(
    x = arithObs, y = arithPred, fill = as.factor(sppCol), color = as.factor(sppCol), alpha = as.factor(sppCol)
  )
) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 26,vjust = 1,margin = margin(t = 0, r = 40, b = 0, l = 0),family = 'Calibri'),
        axis.title.x = element_text(face = 'bold', size = 26,vjust = -1,margin = margin(t = 4, r = 0, b = 0, l = 0),family = 'Calibri'),
        axis.text.x = element_text(face = 'bold',size = 22,color='black', vjust=0.5,family = 'Calibri'),
        axis.text.y = element_text(face = 'bold',size = 22,color='black',family = 'Calibri'),
        legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'),
        legend.text=element_text(face = 'bold',size = 15,color='black',family = 'Calibri'),
        legend.position = 'inside',
        legend.key.spacing.y = unit(0.2, "cm"),
        legend.position.inside = c(0.05,0.7),
        legend.key = element_rect(fill = NA, color = NA),
        legend.background = element_rect(fill = "transparent"),
        axis.ticks.length = unit(0.15, 'cm'))+
  labs(title ='Marginal Predictions', y = expression(SAR['predicted']), x = expression(SAR['estimated'])) +
  theme(plot.title = element_text(hjust = 0.5,size = 24,face = 'bold',family = 'Calibri')) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(shape = 21,size = 4.5,stroke=0.5) +
  annotate("text", x = 0.008, y = 0.12, label = list(bquote(d[r]~ '='~.(sarMeta.IOAarithFix))), size = 8, family = 'Calibri') +
  scale_color_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 'black',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = 'black',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = 'black',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_alpha_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 1.0,
      'steelhead, Snake R.' = 0.7,
      'Chinook, Middle Columbia R.' = 1.0,
      'steelhead, Middle Columbia R.' = 0.7,
      'Chinook, Upper Columbia R.' = 1.0,
      'steelhead, Upper Columbia R.' = 0.7
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_fill_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = '#999999',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = '#E69F00',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = '#56B4E9',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_y_continuous(limits=c(0,0.12),breaks = seq(0,0.12,0.02),labels = scales::percent)+
  scale_x_continuous(limits=c(0,0.12),breaks = seq(0,0.12,0.02),labels = scales::percent)

# *---print figure for review
print(figure.6.4_a)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_a.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(r2_totVar.plot)
# dev.off()
# print(r2_totVar.plot)

# *---marginal predictions; logit scale
figure.6.4_b <- ggplot(
  data = modVerif_avg %>% 
    mutate(
      sppCol = paste0(spp,', ',zone)
    ) %>% 
    filter(
      logitObs != '-Inf'
    ),
  aes(
    x = logitObs, y = logitPred, fill = as.factor(sppCol), color = as.factor(sppCol), alpha = as.factor(sppCol)
  )
) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 26,vjust = 1,margin = margin(t = 0, r = 40, b = 0, l = 0),family = 'Calibri'),
        axis.title.x = element_text(face = 'bold', size = 26,vjust = -1,margin = margin(t = 4, r = 0, b = 0, l = 0),family = 'Calibri'),
        axis.text.x = element_text(face = 'bold',size = 22,color='black', vjust = 0.5,family = 'Calibri'),
        axis.text.y = element_text(face = 'bold',size = 22,color='black',family = 'Calibri'),
        legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'),
        legend.text=element_text(face = 'bold',size = 18,color='black',family = 'Calibri'),
        legend.position = 'none',
        legend.position.inside = c(0.8,0.4),
        axis.ticks.length = unit(0.15, 'cm'))+
  labs(title ='', y = expression(logit(SAR['predicted'])), x = expression(logit(SAR['estimated']))) +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = 'bold',family = 'Calibri')) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(shape = 21,size = 4.5,stroke=0.5) +
  annotate("text", x = -7.6, y = -1.0, label = list(bquote(d[r]~ '='~.(sarMeta.IOAlogitFix))), size = 8, family = 'Calibri') +
  scale_color_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 'black',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = 'black',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = 'black',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_alpha_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 1.0,
      'steelhead, Snake R.' = 0.7,
      'Chinook, Middle Columbia R.' = 1.0,
      'steelhead, Middle Columbia R.' = 0.7,
      'Chinook, Upper Columbia R.' = 1.0,
      'steelhead, Upper Columbia R.' = 0.7
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_fill_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = '#999999',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = '#E69F00',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = '#56B4E9',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_y_continuous(limits=c(-8,-1),breaks = seq(-8,-1,1.0))+
  scale_x_continuous(limits=c(-8,-1),breaks = seq(-8,-1,1.0))


# *---print figure for review
print(figure.6.4_b)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_b.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(r2_totVar.plot)
# dev.off()
# print(r2_totVar.plot)

figure.6.4_c <- ggplot(
  data = modVerif_avg %>% 
    mutate(
      sppCol = paste0(spp,', ',zone)
    ),
  aes(
    x = arithObs, y = farithPred, fill = as.factor(sppCol), color = as.factor(sppCol), alpha = as.factor(sppCol)
  )
) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 26,vjust = 1,margin = margin(t = 0, r = 40, b = 0, l = 0),family = 'Calibri'),
        axis.title.x = element_text(face = 'bold', size = 26,vjust = -1,margin = margin(t = 4, r = 0, b = 0, l = 0),family = 'Calibri'),
        axis.text.x = element_text(face = 'bold',size = 22,color='black', vjust=0.5,family = 'Calibri'),
        axis.text.y = element_text(face = 'bold',size = 22,color='black',family = 'Calibri'),
        legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, 'cm'),
        legend.text=element_text(face = 'bold',size = 18,color='black',family = 'Calibri'),
        legend.position = 'none',
        legend.position.inside = c(0.8,0.4),
        axis.ticks.length = unit(0.15, 'cm'))+
  labs(title ='Conditional Predictions', y = '', x = expression(SAR['estimated'])) +
  theme(plot.title = element_text(hjust = 0.5,size = 24,face = 'bold',family = 'Calibri')) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(shape = 21,size = 4.5,stroke=0.5) +
  annotate("text", x = 0.008, y = 0.12, label = list(bquote(d[r]~ '='~.(sarMeta.IOAarithRand))), size = 8, family = 'Calibri') +
  scale_color_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 'black',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = 'black',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = 'black',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_alpha_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 1.0,
      'steelhead, Snake R.' = 0.7,
      'Chinook, Middle Columbia R.' = 1.0,
      'steelhead, Middle Columbia R.' = 0.7,
      'Chinook, Upper Columbia R.' = 1.0,
      'steelhead, Upper Columbia R.' = 0.7
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_fill_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = '#999999',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = '#E69F00',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = '#56B4E9',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_y_continuous(limits=c(0,0.12),breaks = seq(0,0.12,0.02),labels = scales::percent)+
  scale_x_continuous(limits=c(0,0.12),breaks = seq(0,0.12,0.02),labels = scales::percent)


# *---print figure for review
print(figure.6.4_c)

# *---conditional predictions; logit scale
figure.6.4_d <- ggplot(
  data = modVerif_avg %>% 
    mutate(
      sppCol = paste0(spp,', ',zone)
    ) %>% 
    filter(
      logitObs != '-Inf'
    ),
  aes(
    x = logitObs, y = flogitPred, fill = as.factor(sppCol), color = as.factor(sppCol), alpha = as.factor(sppCol)
  )
) +
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = 'black'),
        axis.title.y = element_text(face = 'bold', size = 26,vjust = 1,margin = margin(t = 0, r = 40, b = 0, l = 0),family = 'Calibri'),
        axis.title.x = element_text(face = 'bold', size = 26,vjust = -1,margin = margin(t = 4, r = 0, b = 0, l = 0),family = 'Calibri'),
        axis.text.x = element_text(face = 'bold',size = 22,color='black', vjust=0.5,family = 'Calibri'),
        axis.text.y = element_text(face = 'bold',size = 22,color='black',family = 'Calibri'),
        legend.title = element_blank(),
        plot.margin = margin(0.1, 0.1, 0.15, 0.1, 'cm'),
        legend.text=element_text(face = 'bold',size = 18,color='black',family = 'Calibri'),
        legend.position = 'none',
        legend.position.inside = c(0.8,0.4),
        axis.ticks.length = unit(0.15, 'cm'))+
  labs(title ='', y = '', x = expression(logit(SAR['estimated']))) +
  theme(plot.title = element_text(hjust = 0.5,size = 16,face = 'bold',family = 'Calibri')) +
  geom_abline(intercept = 0, slope = 1)+
  geom_point(shape = 21,size = 4.5,stroke=0.5) +
  annotate("text", x = -7.6, y = -1.0, label = list(bquote(d[r]~ '='~.(sarMeta.IOAlogitRand))), size = 8, family = 'Calibri') +
  scale_color_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 'black',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = 'black',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = 'black',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_alpha_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = 1.0,
      'steelhead, Snake R.' = 0.7,
      'Chinook, Middle Columbia R.' = 1.0,
      'steelhead, Middle Columbia R.' = 0.7,
      'Chinook, Upper Columbia R.' = 1.0,
      'steelhead, Upper Columbia R.' = 0.7
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_fill_manual(
    name = '',
    values = c(
      'Chinook, Snake R.' = '#999999',
      'steelhead, Snake R.' = '#999999',
      'Chinook, Middle Columbia R.' = '#E69F00',
      'steelhead, Middle Columbia R.' = '#E69F00',
      'Chinook, Upper Columbia R.' = '#56B4E9',
      'steelhead, Upper Columbia R.' = '#56B4E9'
    ),
    breaks = c(
      'steelhead, Middle Columbia R.',
      'Chinook, Middle Columbia R.',
      'steelhead, Upper Columbia R.',
      'Chinook, Upper Columbia R.',
      'steelhead, Snake R.',
      'Chinook, Snake R.'
    )
  ) +
  scale_y_continuous(limits=c(-8,-1),breaks = seq(-8,-1,1.0))+
  scale_x_continuous(limits=c(-8,-1),breaks = seq(-8,-1,1.0))

# *---print figure for review
print(figure.6.4_d)

# *---combined plot
figure.6.4 <- ggarrange(
  figure.6.4_a,
  figure.6.4_c,
  figure.6.4_b,
  figure.6.4_d,
  ncol = 2,
  nrow = 2
)

# *---print figure for review
print(figure.6.4)

# simulations (these were presented in Chapter 9 of the 2025 report)
# 1---specify simulation parameters
# a---species
# *---0 = Chinook
# *---1 = steelhead
sppNum <- 0

# b---number of replications
reps <- 1000

# 2---estimate distribution (normal) parameters for random effects
# a---first component model
# *---observation-level random effect
paramObs_mod1 <- enorm(
  method = 'mle',
  randEffEsts_mod1$`factor(obs.id)` %>%
    pull(intrcpt)
)

# *---steelhead migration year random effect
paramSTLmy_mod1 <- enorm(
  method = 'mle',
  randEffEsts_mod1$`~mig.yr | factor(spp.code)` %>%
    tibble::rownames_to_column() %>%
    filter(
      grepl('ST', rowname)
    ) %>%
    pull(intrcpt)
)

# *---Chinook migration year random effect
paramCHNmy_mod1 <- enorm(
  method = 'mle',
  randEffEsts_mod1$`~mig.yr | factor(spp.code)` %>%
    tibble::rownames_to_column() %>%
    filter(
      grepl('CH', rowname)
    ) %>%
    pull(intrcpt)
)

# *---CSS group random effect
# paramCSSgrp_mod1 <- enorm(
#   method = 'mle',
#   randEffEsts_mod1$orig_grp %>%
#     pull(intrcpt)
# )

# b---second component model
# *---observation-level random effect
paramObs_mod2 <- enorm(
  method = 'mle',
  randEffEsts_mod2$`factor(obs.id)` %>%
    pull(intrcpt) 
)

# *---steelhead migration year random effect
paramSTLmy_mod2 <- enorm(
  method = 'mle',
  randEffEsts_mod2$`~mig.yr | factor(spp.code)` %>% 
    tibble::rownames_to_column() %>% 
    filter(
      grepl('ST', rowname)
    ) %>%
    pull(intrcpt) 
)

# *---Chinook migration year random effect 
paramCHNmy_mod2 <- enorm(
  method = 'mle',
  randEffEsts_mod2$`~mig.yr | factor(spp.code)` %>% 
    tibble::rownames_to_column() %>% 
    filter(
      grepl('CH', rowname)
    ) %>%
    pull(intrcpt) 
)

# *---CSS group random effect 
# paramCSSgrp_mod2 <- enorm(
#   method = 'mle',
#   randEffEsts_mod2$orig_grp %>%
#     pull(intrcpt) 
# )

# 3---simulation routines
# a---first component model
# i---generate predictions from fixed effects
fixPredFun_mod1 <- function(sppNum,i){
  sarMeta.pred <- predict(
    predMod1,
    newmods = cbind(seq(2,20,1),c(rep(i,19))),
    # transf = transf.ilogit,
    addx = TRUE
  ) %>% 
    as.data.frame()
  # print(sarMeta.pred)
}

# ii---create output data frame to store replications
fp.out.preds_mod1 <- lapply(
  seq(0.5,3,0.5),
  function(i){
    cbind(c(seq(2,20,1)),c(rep(i,19)))
  }
)

fp.out.preds_mod1 <- do.call(
  rbind,
  fp.out.preds_mod1
) %>% 
  as.data.frame()

# iii---apply function across a range of PITPH values
# *---generate predictions based on fixed effects
outFun_mod1 <- function(){
  fp.out_mod1 <- lapply(
    seq(0.5,3,0.5), 
    fixPredFun_mod1,
    sppNum = sppNum
  ) %>% 
    bind_rows(.id = 'column_label')
  
  # *---add-in draws from simulated distributions (random-normal) of random effects
  fp.out2_mod1 <- fp.out_mod1 %>% 
    mutate(
      fObs.id_mod1 = rnorm(
        nrow(fp.out_mod1),
        mean = paramObs_mod1$parameters[1],
        sd = paramObs_mod1$parameters[2]
      ),
      fMig.yr_mod1 = rnorm(
        nrow(fp.out_mod1), 
        mean = if_else(sppNum == 1, paramSTLmy_mod1$parameters[1], paramCHNmy_mod1$parameters[1]), 
        sd = if_else(sppNum == 1, paramSTLmy_mod1$parameters[2], paramCHNmy_mod1$parameters[2])
      )
      # fspp.code_mod1 = if_else(
      #   sppNum == 1, 
      #   as.numeric(
      #     randEffEsts_mod1$`factor(spp.code)` %>%
      #       pull(intrcpt) %>% 
      #       tail(1)
      #   ),
      # as.numeric(
      #   randEffEsts_mod1$`factor(spp.code)` %>%
      #     pull(intrcpt) %>% 
      #     head(1)
      # )
      # )
      # forig_grp_mod1 = rnorm(
      #   nrow(fp.out_mod1), 
      #   mean = paramCSSgrp_mod1$parameters[1], 
      #   sd = paramCSSgrp_mod1$parameters[2]
      # )
    ) %>% 
    rowwise() %>%
    mutate(
      logitPred_mod1 = rowSums(
        pick(
          pred, 
          fObs.id_mod1, 
          fMig.yr_mod1
          # fspp.code_mod1
          # forig_grp_mod1
        ), 
        na.rm = FALSE
      )
    ) %>%
    # *---back-transform predictions
    mutate(
      arithPred_mod1 = inv.logit(logitPred_mod1)
    ) %>% 
    as.data.frame()
  
  # *---combine data into final (mean) predictions data frame  
  fp.out.preds_mod1 <<- cbind(
    fp.out.preds_mod1,
    fp.out2_mod1[,ncol(fp.out2_mod1)]
  )
}

# iv---generate replicates
replicate(reps,outFun_mod1())

# v---create data set of weighted (AICc) predictions
mod1.out <- fp.out.preds_mod1

colnames(mod1.out) <- make.names(
  colnames(mod1.out),
  unique = TRUE
)

mod1.out <- mod1.out %>% 
  mutate(
    weight = as.numeric(
      avgMod$msTable[5] %>% 
        head(1)
    )
  ) %>% 
  relocate(
    weight,
    .after = V2
  ) %>% 
  mutate(across(everything(), ~ .x * weight, .names = 'wt_{.col}')) %>% 
  dplyr::select(
    -c(
      wt_V1,
      wt_V2,
      wt_weight
    )
  ) %>% 
  dplyr::select(
    starts_with('V') | starts_with('wt')
  )

# b---second component model
# i---generate predictions from fixed effects
fixPredFun_mod2 <- function(sppNum,i){
  sarMeta.pred <- predict(
    predMod2,
    newmods = cbind(rep(sppNum,19),seq(2,20,1),c(rep(i,19))),
    # transf = transf.ilogit,
    addx = TRUE
  ) %>% 
    as.data.frame()
  # print(sarMeta.pred)
}

# ii---create output data frame to store replications
fp.out.preds_mod2 <- lapply(
  seq(0.5,3,0.5),
  function(i){
    cbind(c(seq(2,20,1)),c(rep(i,19)))
  }
) 

fp.out.preds_mod2 <- do.call(
  rbind,
  fp.out.preds_mod2
) %>% 
  as.data.frame()

# iii---apply function across a range of PITPH values
# *---generate predictions based on fixed effects
outFun_mod2 <- function(){
  fp.out_mod2 <- lapply(
    seq(0.5,3,0.5), 
    fixPredFun_mod2,
    sppNum = sppNum
  ) %>% 
    bind_rows(.id = 'column_label')
  
  # *---add-in draws from simulated distributions (random-normal) of random effects
  fp.out2_mod2 <- fp.out_mod2 %>% 
    mutate(
      fObs.id_mod2 = rnorm(
        nrow(fp.out_mod2),
        mean = paramObs_mod2$parameters[1],
        sd = paramObs_mod2$parameters[2]
      ),
      fMig.yr_mod2 = rnorm(
        nrow(fp.out_mod2), 
        mean = if_else(sppNum == 1, paramSTLmy_mod2$parameters[1], paramCHNmy_mod2$parameters[1]), 
        sd = if_else(sppNum == 1, paramSTLmy_mod2$parameters[2], paramCHNmy_mod2$parameters[2])
      )
      # fspp.code_mod2 = if_else(
      #   sppNum == 1, 
      #   as.numeric(
      #     randEffEsts_mod2$`factor(spp.code)` %>%
      #       pull(intrcpt) %>% 
      #       tail(1)
      #   ),
      #   as.numeric(
      #     randEffEsts_mod2$`factor(spp.code)` %>%
      #       pull(intrcpt) %>% 
      #       head(1)
      #   )
      # )
      # forig_grp_mod2 = rnorm(
      #   nrow(fp.out_mod2), 
      #   mean = paramCSSgrp_mod2$parameters[1], 
      #   sd = paramCSSgrp_mod2$parameters[2]
      # )
    ) %>% 
    rowwise() %>%
    mutate(
      logitPred_mod2 = rowSums(
        pick(
          pred, 
          fObs.id_mod2, 
          fMig.yr_mod2
          # forig_grp_mod2
        ), 
        na.rm = FALSE
      )
    ) %>%
    # *---back-transform predictions
    mutate(
      arithPred_mod2 = inv.logit(logitPred_mod2)
    ) %>% 
    as.data.frame()
  
  # *---combine data into final (mean) predictions data frame  
  fp.out.preds_mod2 <<- cbind(
    fp.out.preds_mod2,
    fp.out2_mod2[,ncol(fp.out2_mod2)]
  )
}

# iv---generate replicates
replicate(reps,outFun_mod2())

# v---create data set of weighted (AICc) predictions
mod2.out <- fp.out.preds_mod2

colnames(mod2.out) <- make.names(
  colnames(mod2.out),
  unique = TRUE
)

mod2.out <- mod2.out %>% 
  mutate(
    weight = as.numeric(
      avgMod$msTable[5] %>% 
        tail(1)
    )
  ) %>% 
  relocate(
    weight,
    .after = V2
  ) %>% 
  mutate(across(everything(), ~ .x * weight, .names = 'wt_{.col}')) %>% 
  dplyr::select(
    -c(
      wt_V1,
      wt_V2,
      wt_weight
    )
  ) %>% 
  dplyr::select(
    starts_with('V') | starts_with('wt')
  )

# 3---combine output from component models into ensamble predictions
comb.mod <- bind_cols(
  mod1.out$V1,
  mod1.out$V2,
  mod1.out[,3:ncol(mod1.out)] + mod2.out[,3:ncol(mod2.out)]
)
# a---generate lower hdi for each wtt x pitph combination (rows)
lhdi.out <- lapply(
  seq(1.0,nrow(comb.mod),1.0),
  function(i){
    HDInterval::hdi(as.numeric(comb.mod[i,3:ncol(comb.mod)]))[1]
  }
) %>% 
  bind_rows() %>% 
  as.data.frame() %>% 
  rename(
    lhdi = lower
  )

# b---generate upper hdi for each wtt x pitph combination (rows)
uhdi.out <- lapply(
  seq(1.0,nrow(comb.mod),1.0),
  function(i){
    HDInterval::hdi(as.numeric(comb.mod[i,3:ncol(comb.mod)]))[2]
  }
) %>% 
  bind_rows() %>%  
  as.data.frame() %>%  
  rename(
    uhdi = upper
  )

# c---create summary output data frame
out.dat <- data.frame(
  wtt = comb.mod[,1],
  pitph = comb.mod[,2],
  mean = rowMeans(comb.mod[,3:ncol(comb.mod)]),
  median = rowMedians(as.matrix(comb.mod[,3:ncol(comb.mod)])),
  lhdi = pull(lhdi.out),
  uhdi = pull(uhdi.out)
)

###--- save (as .rds) summary output frame
saveRDS(
  out.dat, 
  paste0(
    "Output\\Datasets\\simsOutwAR_",
    if_else(sppNum == 1,"stl","chn"),
    reps,
    ".rds"
  )
)

# d---generate heat map
# i---read-in data
sarMeta_simsPlot.dat <- readRDS(
  paste0(
    "Output\\Datasets\\simsOutwAR_",
    if_else(sppNum == 1,"stl","chn"),
    reps,
    ".rds"
  )
)

# ii---generate plot
sarMeta_HMp.plot <- ggplot(sarMeta_simsPlot.dat, aes(x=pitph, y=wtt, fill=mean)) + 
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.y = element_text(face = 'bold', size = 14, margin = margin(t = 0, r = 5, b = 0, l = 0),family = 'Calibri'),
    axis.title.x = element_text(face = 'bold', size = 14, margin = margin(t = 0, r = 0, b = 0, l = 0),family = 'Calibri'),
    axis.text.x = element_text(face = 'bold',size = 12,color='black', hjust=0.5,family = 'Calibri'),
    axis.text.y = element_text(face = 'bold',size = 12,color='black',family = 'Calibri'),
    legend.title = element_text(size = 12, hjust = 0.5, face = 'bold'),
    plot.title = element_text(hjust = 0,size = 14,face = 'bold',family = 'Calibri')) +
  labs(title ='', y = 'WTT', x = 'PITPH') +
  geom_tile() +
  geom_text(aes(label=paste0(format(round(mean*100, 1), trim = TRUE), '%')), fontface = 2, size = 4.0, family = 'Calibri', color = 'black') +
  scale_fill_gradientn(colors = rainbow(if_else(sppNum == 1,12.0,3.0)), name = 'Mean SAR', breaks = c(0.01,0.02,0.03,0.04,0.05), labels = c('1.0%','2.0%','3.0%','4.0%','5.0%')) +
  scale_y_reverse(breaks = seq(2.0,20.0,1), expand = c(0,0)) +
  scale_x_continuous(breaks = seq(0.5,3.0,0.5), position = 'top', expand = c(0,0))

# iii--- print figure for review
print(sarMeta_HMp.plot)