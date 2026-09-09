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

# x---summarize marine covariates
## NPGO
npgo.dat <- read_lines("https://o3d.org/npgo/data/NPGO.txt") %>% 
  as.data.frame() %>% 
  tail(-26) %>% 
  rename(
    'data' = 1
  ) %>% 
  cSplit('data', sep=" ", type.convert=FALSE) %>% 
  rename(
    'year' = 1,
    'month' = 2,
    'npgo' = 3
  ) %>% 
  mutate(
    year = as.numeric(year),
    month = as.numeric(month),
    npgo = as.numeric(npgo)
  ) %>% 
  filter(
    month == 5
  ) %>% 
  rename(
    'npgo.may' = 'npgo'
  ) %>% 
  dplyr::select(
    year,
    npgo.may
  )

## Upwelling
upwel.dat <- read.csv('https://psl.noaa.gov/data/timeseries/month/data/bakun.45125.csv') %>% 
  rename(
    'date' = 1,
    'upwel' = 2
  ) %>% 
  mutate(
    year = year(date),
    month =  month(date),
    day = day(date)
  ) %>% 
  dplyr::select(
    -c(
      date,
      day
    )
  ) %>% 
  relocate(
    upwel,
    .after = month
  ) %>% 
  filter(
    month == 4
  ) %>% 
  rename(
    'upwel.apr' = 'upwel'
  )%>% 
  dplyr::select(
    year,
    upwel.apr
  )

## SST
sst.dat <- read.csv('https://www.psl.noaa.gov/thredds/ncss/grid/Datasets/icoads/1degree/global/enh/sst.mean.nc?var=sst&latitude=46.5&longitude=%20-124.5&time_start=1960-01-01T00:00:00Z&time_end=2026-04-01T00:00:00Z&&&accept=csv') %>%
  rename(
    'lat' = 3,
    'long' = 4,
    'sst' = 5
  ) %>% 
  mutate(
    year = year(time),
    month = month(time)
  ) %>% 
  relocate(
    sst,
    .after = month
  ) %>% 
  dplyr::select(
    -c(
      time,
      station
    )
  ) %>% 
  filter(
    month >= 5 & month <=9
  ) %>% 
  group_by(year) %>%  
  dplyr::summarize(
    mean_sst = mean(sst)
  ) %>% 
  rename(
    'sst.maysep' = 'mean_sst'
  )%>% 
  dplyr::select(
    year,
    sst.maysep
  )

## combined covs
combMarCovs <- upwel.dat %>% 
  full_join(
    npgo.dat,
    by = 'year'
  ) %>% 
  full_join(
    sst.dat,
    by = 'year'
  ) %>% 
  rename(
    'mig.yr' = "year"
  )

# x1---append marine covariates to summary effects table
sar.meta.es <- sar.meta.es %>% 
  left_join(
    combMarCovs,
    by = 'mig.yr'
)

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

# c---conduct MMI and create and output summary table
table.6.4 <- dredge(
  fullMod,
  trace=2
) %>%
  as.data.frame() %>%
  mutate(
    `factor(spp.code)` = c(
      update(
        fullMod,
        ~ factor(spp.code) + wtt + pitph
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      NA,
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

# d---specify best model
# i---update full model based on model selection table
bestMod <- update(
  eval(
    parse(
      text = top.mod
    )
  ),
  mods = ~ factor(spp.code) + wtt + pitph
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
    param = c(bestMod$beta[2:4,]),
    lcl = c(
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
        'Spp.',
        'PITPH',
        'WTT'
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
  scale_x_continuous(limits=c(-2.8,2.8),
                     breaks = seq(-2.8,2.8,0.7),
                     labels = function(x) format(x, scientific = FALSE)
  ) +
geom_magnify(
  from = c(xmin = -0.25, xmax = 0.25, ymin = 0.9, ymax = 2.1),
  to = c(-1.2 - 0.5, -1.2 + 0.5,0.5, 1.5),
  shadow = TRUE
)

# iv---print figure for review
print(figure.6.3)

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

# variance components -----------------------------------------------------
#c--- variance component tables (i.e., Table 6.6)
#i--- full model (all moderators)
fullVC.mod <- update(
  redREmodAR,
  mods = ~ factor(spp.code) + wtt + pitph
)

fullVC.mod2 <- update(
  redREmodAR,
  mods = ~ wtt + pitph
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
  ~ factor(spp.code)
)

redVC_nm.mod2 <- update(
  redREmodAR,
  ~ 1
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
        "Reduced model (mods.:MCS1 = Spp./MCS2 = none)", 
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
  ~ factor(spp.code) + pitph
)

redVC_pitph.mod2 <- update(
  redREmodAR,
  ~ pitph
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
        'Reduced model (mods.:MCS1 = Spp. + PITPH/MCS2 = PITPH)', 
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
  ~ factor(spp.code) + wtt
)

redVC_wtt.mod2 <- update(
  redREmodAR,
  ~ wtt
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
        'Reduced model (mods.:MCS1 = Spp. + WTT/MCS2 = WTT)',
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

# save_as_image(
#   path = 'Output\\Tables\\Table 6.7.png',
#   res = 2000
# )

# model verification (average predictions) --------------------------------
# a---specify component models
# i---first component model
predMod1 <- update(
  redREmod,
  ~ factor(spp.code) + wtt + pitph
)

# ii---second component model
predMod2 <- update(
  redREmod,
  ~ pitph + wtt
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
        tail(1)
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
# print(figure.6.4_a)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_a.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.4_a)
# dev.off()
# print(figure.6.4_a)

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
# print(figure.6.4_b)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_b.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.4_b)
# dev.off()
# print(figure.6.4_b)

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
# print(figure.6.4_c)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_c.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.4_c)
# dev.off()
# print(figure.6.4_c)

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
# print(figure.6.4_d)

#v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_d.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.4_d)
# dev.off()
# print(figure.6.4_d)

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
# print(figure.6.4)

# *---export combined figure
# png(filename='Output/Figures/Verification/figure.6.4.png',
#     type='cairo',
#     units='in',
#     width=16,
#     height=12,
#     res=300)
# print(figure.6.4)
# dev.off()
# print(figure.6.4)





















# NEW WORK FOR 2026 - ASSESS MARINE COVARIATES
marCov.mod <- update(
  redREmod,
  mods = ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may + sst.maysep
)

# c---conduct MMI and create and output summary table
table.6.8_marCov <- dredge(
  marCov.mod,
  trace=2
) %>%
  subset(
    delta<4.0,
    recalc.weights = TRUE
    ) %>% 
  as.data.frame() %>%
  mutate(
    `factor(spp.code)` = c(
      NA,
      update(
        marCov.mod,
        ~ factor(spp.code) + wtt + pitph + npgo.may
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      NA,
      update(
        marCov.mod,
        ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      update(
        marCov.mod,
        ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may + sst.maysep
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      update(
        marCov.mod,
        ~ factor(spp.code) + wtt + pitph + npgo.may + sst.maysep
      )[1] |>  
        as.data.frame() |>  
        head(2) |>  
        tail(1) |>  
        as.numeric(),
      NA,
      NA
    )
  ) %>% 
  relocate(`factor(spp.code)`, .before = pitph) %>%
  relocate(pitph, .after = `factor(spp.code)`) %>% 
  relocate(wtt, .after = pitph) %>% 
  relocate(npgo.may, .after = wtt) %>% 
  relocate(upwel.apr, .after = npgo.may) %>% 
  relocate(sst.maysep, .after = upwel.apr) %>% 
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
    npgo.may = 'NPGO (May)',
    upwel.apr = 'Upwel. (Apr.)',
    sst.maysep = 'SST (May-Sep.)',
    df = 'No. params.',
    logLik = 'log-liklihood'
  ) %>%
  mk_par( 
    j = 10,
    i = 1,
    value = as_paragraph(
      'AIC',
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 11,
    i = 1,
    value = as_paragraph(
      paste0('\u394','AIC'),
      as_sub('c')
    ),
    part = 'header'
  ) %>% 
  mk_par( 
    j = 12,
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
    npgo.may = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    upwel.apr = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    sst.maysep = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
    logLik = function(x) ifelse(is.na(x),'', formatC(x,digits = 3, format = 'f')),
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
#   path = 'Output\\Tables\\Table 6.8.png',
#   res = 2000
# )

# e---plot model-averaged and best model coefficients and CIs
# i---estimate model-averaged coefficients
avgMod_marCov <- dredge(
  marCov.mod,
  trace = 3
) %>% 
  subset(
    delta <= 4.0, # specify confidence set to include models where dAICc <= 4.0
    recalc.weights = TRUE
  ) %>% 
  model.avg()

#ii---data manipulation
figure.6.5_marCov <- bind_rows(
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
        'Spp.',
        'PITPH',
        'WTT'
      )
    ) %>% 
    tail(4) %>% 
    mutate(
      selType = as.factor(rep('w/o Mar. Covs.',4))
    ) %>% 
    relocate(
      selType,
      .before = param.lbs
    ),
  coefTable(
    avgMod_marCov, 
    full = TRUE
  ) %>%
    as.data.frame() %>% 
    dplyr::select(
      Estimate
    ) %>% 
    tibble::rownames_to_column('param.lbs') %>% 
    left_join(
      confint(
        avgMod_marCov, 
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
        'NPGO (May)',
        'PITPH',
        'WTT',
        'Spp.',
        'Upwel. (Apr.)',
        'SST (May-Sep.)' 
      )
    ) %>% 
    tail(7) %>% 
    mutate(
      selType = as.factor(rep('w/ Mar. Covs.',7))
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
  scale_colour_manual(name='Transplanted', values = c('w/o Mar. Covs.' = '#D55E00','w/ Mar. Covs.' = 'black'), breaks = c('w/o Mar. Covs.', 'w/ Mar. Covs.')) +
  geom_vline(xintercept = 0,linetype = 'dashed') +
  geom_errorbar(aes(xmax=ucl,xmin=lcl), height = 0, position = position_dodge(width = -0.3), linewidth = 1.0, orientation = 'y') +
  geom_point(size = 3.5, position = position_dodge(width = -0.3)) +
  # geom_point(size = 4,shape = 1, position = position_dodge(width = 0.3)) +
  labs(y = 'Parameter',x = 'Estimate') +
  scale_y_discrete(limits = c('SST (May-Sep.)','Upwel. (Apr.)','NPGO (May)','WTT','PITPH','Spp.')) + 
  scale_x_continuous(limits=c(-0.7,0.7),
                     breaks = seq(-0.7,0.7,0.2),
                     labels = function(x) format(x, scientific = FALSE)
  ) 

# iv---print figure for review
print(figure.6.5_marCov)

# v--- output figure
png(filename='Output\\Figures\\Figure 6.5.png',
    type='cairo',
    units='in',
    width=8,
    height=6,
    res=300)
print(figure.6.5_marCov)
dev.off()
print(figure.6.5_marCov)










# model verification (average predictions) --------------------------------
## develop training set
valid.train <- sar.meta.es %>%
  group_by(css.grp) %>%
  slice_head(prop = 0.60) %>%
  ungroup()

## develop testing set
valid.test <- anti_join(
  sar.meta.es,
  valid.train,
  by = c('css.grp','mig.yr')
)

# a---specify component models
## fit full model to training set
marCov.mod <- update(
  redREmod,
  mods = ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may + sst.maysep,
  data = valid.train
)

## estimate model average coefficients
marCov.avg <- dredge(
  marCov.mod,
  trace = 3
) %>% 
  subset(
    delta <= 4.0, # specify confidence set to include models where dAICc <= 4.0
    recalc.weights = TRUE
  ) %>% 
  {print(.)} %>% 
  model.avg()

## fit component models to training set
predMod1_marCov <- update(
  marCov.mod,
  ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may,
  data = valid.train
)

predMod2_marCov <- update(
  marCov.mod,
  ~ wtt + pitph + upwel.apr + npgo.may,
  data = valid.train
)

predMod3_marCov <- update(
  marCov.mod,
  ~ factor(spp.code) + wtt + pitph + upwel.apr + npgo.may + sst.maysep,
  data = valid.train
)

predMod4_marCov <- update(
  marCov.mod,
  ~ factor(spp.code) + wtt + pitph + npgo.may,
  data = valid.train
)

predMod5_marCov <- update(
  marCov.mod,
  ~ wtt + pitph + npgo.may,
  data = valid.train
)

predMod6_marCov <- update(
  marCov.mod,
  ~ wtt + pitph + upwel.apr + npgo.may + sst.maysep,
  data = valid.train
)

predMod7_marCov <- update(
  marCov.mod,
  ~ factor(spp.code) + wtt + pitph + npgo.may + sst.maysep,
  data = valid.train
)

predMod8_marCov <- update(
  marCov.mod,
  ~ wtt + pitph + npgo.may + sst.maysep,
  data = valid.train
)

# b---generate weighted predictions
# i---first component model
# *---generate predictions data set
cMod1.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    spp.code,
    wtt,
    pitph,
    upwel.apr,
    npgo.may,
  ) %>%
  mutate(
    spp.code = case_match(
      spp.code,
      'CH' ~ '0',
      'ST' ~ '1',
      .default = spp.code
    )
  ) %>%
  mutate(
    spp.code = as.numeric(spp.code)
  )

## generate predictions and manipulate data
Mod1_marCov.pred <- predict(
  predMod1_marCov,
  newmods = as.matrix(cMod1.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod1_marCov.verif <- Mod1_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )










cMod2.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    wtt,
    pitph,
    upwel.apr,
    npgo.may,
  )
  # mutate(
  #   spp.code = case_match(
  #     spp.code,
  #     'CH' ~ '0',
  #     'ST' ~ '1',
  #     .default = spp.code
  #   )
  # ) %>%
  # mutate(
  #   spp.code = as.numeric(spp.code)
  # )

## generate predictions and manipulate data
Mod2_marCov.pred <- predict(
  predMod2_marCov,
  newmods = as.matrix(cMod2.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod2_marCov.verif <- Mod2_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(2) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---third component model
# *---generate predictions data set
cMod3.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    spp.code,
    wtt,
    pitph,
    upwel.apr,
    npgo.may,
    sst.maysep
  ) %>%
  mutate(
    spp.code = case_match(
      spp.code,
      'CH' ~ '0',
      'ST' ~ '1',
      .default = spp.code
    )
  ) %>%
  mutate(
    spp.code = as.numeric(spp.code)
  )

## generate predictions and manipulate data
Mod3_marCov.pred <- predict(
  predMod3_marCov,
  newmods = as.matrix(cMod3.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod3_marCov.verif <- Mod3_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(3) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---fourth component model
# *---generate predictions data set
cMod4.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    spp.code,
    wtt,
    pitph,
    npgo.may
  ) %>%
  mutate(
    spp.code = case_match(
      spp.code,
      'CH' ~ '0',
      'ST' ~ '1',
      .default = spp.code
    )
  ) %>%
  mutate(
    spp.code = as.numeric(spp.code)
  )

## generate predictions and manipulate data
Mod4_marCov.pred <- predict(
  predMod4_marCov,
  newmods = as.matrix(cMod4.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod4_marCov.verif <- Mod4_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(4) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---fifth component model
# *---generate predictions data set
cMod5.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    wtt,
    pitph,
    npgo.may
  )

## generate predictions and manipulate data
Mod5_marCov.pred <- predict(
  predMod5_marCov,
  newmods = as.matrix(cMod5.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod5_marCov.verif <- Mod5_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(5) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---sixth component model
# *---generate predictions data set
cMod6.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    wtt,
    pitph,
    upwel.apr,
    npgo.may,
    sst.maysep
  )

## generate predictions and manipulate data
Mod6_marCov.pred <- predict(
  predMod6_marCov,
  newmods = as.matrix(cMod6.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod6_marCov.verif <- Mod6_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(6) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---seventh component model
# *---generate predictions data set
cMod7.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    spp.code,
    wtt,
    pitph,
    npgo.may,
    sst.maysep
  ) %>%
  mutate(
    spp.code = case_match(
      spp.code,
      'CH' ~ '0',
      'ST' ~ '1',
      .default = spp.code
    )
  ) %>%
  mutate(
    spp.code = as.numeric(spp.code)
  )

## generate predictions and manipulate data
Mod7_marCov.pred <- predict(
  predMod7_marCov,
  newmods = as.matrix(cMod7.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod7_marCov.verif <- Mod7_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(7) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )




# i---eigth component model
# *---generate predictions data set
cMod8.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    wtt,
    pitph,
    npgo.may,
    sst.maysep
  )

## generate predictions and manipulate data
Mod8_marCov.pred <- predict(
  predMod8_marCov,
  newmods = as.matrix(cMod8.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod8_marCov.verif <- Mod8_marCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
  ) %>% 
  mutate(
    spp=replace(spp, spp=='CH', 'Chinook')
  ) %>% 
  mutate(
    spp=replace(spp, spp=='ST', 'steelhead')
  ) %>% 
  mutate(
    weight = as.numeric(
      marCov.avg$msTable[5] %>% 
        head(8) %>% 
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )



# iii---combine predictions from component models
marCov.verif_avg <- bind_rows(
  Mod1_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(), 
  Mod2_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod3_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod4_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod5_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod6_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod7_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
  Mod8_marCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(),
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
    zone = Mod1_marCov.verif$zone,
    spp = Mod1_marCov.verif$spp,
    arithObs = Mod1_marCov.verif$arithObs,
    logitObs = Mod1_marCov.verif$logitObs
  ) %>% 
  rename(
    arithPred = wt_arithPred,
    logitPred = wt_logitPred
  ) %>% 
  as.data.frame()


marCov.IOAarithFix <- modStats(
  marCov.verif_avg,
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

# ii---generate plots
# *---marginal predictions; arithmetic scale
figure.6.4_a_marCov <- ggplot(
  data = marCov.verif_avg %>% 
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
  annotate("text", x = 0.008, y = 0.12, label = list(bquote(d[r]~ '='~.(marCov.IOAarithFix))), size = 8, family = 'Calibri') +
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
  scale_y_continuous(limits=c(0,0.12),breaks = seq(0,0.06,0.12),labels = scales::percent)+
  scale_x_continuous(limits=c(0,0.12),breaks = seq(0,0.12,0.12),labels = scales::percent)

# *---print figure for review
print(figure.6.4_a_marCov)

 #v--- output figure
# png(filename='Output\\Figures\\Verification\\Figure.6.4_a_marCov.png',
#     type='cairo',
#     units='in',
#     width=8,
#     height=6,
#     res=300)
# print(figure.6.4_a_marCov)
# dev.off()
# print(figure.6.4_a_marCov)











# a---specify component models
# i---first component model
predMod1_nMarCov <- update(
  redREmod,
  ~ factor(spp.code) + wtt + pitph,
  data = valid.train 
)

# ii---second component model
predMod2_nMarCov <- update(
  redREmod,
  ~ pitph + wtt,
  data = valid.train
)

# b---generate weighted predictions
# i---first component model
# *---generate predictions data set
ncMod1.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    spp.code,
    wtt,
    pitph
  ) %>%
  mutate(
    spp.code = case_match(
      spp.code,
      'CH' ~ '0',
      'ST' ~ '1',
      .default = spp.code
    )
  ) %>%
  mutate(
    spp.code = as.numeric(spp.code)
  )

## generate predictions and manipulate data
Mod1_nMarCov.pred <- predict(
  predMod1_nMarCov,
  newmods = as.matrix(ncMod1.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod1_nMarCov.verif <- Mod1_nMarCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
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
    wt_logitPred = pred * weight
  )




# i---second component model
# *---generate predictions data set
ncMod2.dat <- as.data.frame(valid.test) %>% 
  dplyr::select(
    wtt,
    pitph
  )

## generate predictions and manipulate data
Mod2_nMarCov.pred <- predict(
  predMod2_nMarCov,
  newmods = as.matrix(ncMod2.dat),
  addx = TRUE
) %>% 
  as.data.frame() %>%   
  mutate(
    zone = valid.test %>%
      dplyr::select(
        zone,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        zone
      ) %>% 
      unlist() %>% 
      unname(),
    orig_grp = valid.test %>%
      dplyr::select(
        orig_grp,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        orig_grp
      ) %>% 
      unlist() %>% 
      unname(),
    mig.yr = valid.test %>%
      dplyr::select(
        mig.yr,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        mig.yr
      ) %>% 
      unlist() %>% 
      unname(),
    obs.id = valid.test %>%
      dplyr::select(
        obs.id,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        obs.id
      ) %>% 
      unlist() %>% 
      unname(),
    spp.code = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname()
  ) %>%
  mutate(
    arithPred = inv.logit(pred)
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
Mod2_nMarCov.verif <- Mod2_nMarCov.pred %>%  
  mutate(
    arithObs = valid.test$sar.est/100,
    spp = valid.test %>%
      dplyr::select(
        spp.code,
        cases
      ) %>% 
      # na.omit() %>% 
      dplyr::select(
        spp.code
      ) %>% 
      unlist() %>% 
      unname(),
    logitObs = logit(valid.test$sar.est/100)
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
        tail(1)
    )
  ) %>% 
  mutate(
    wt_arithPred = arithPred * weight,
    wt_logitPred = pred * weight
  )









# iii---combine predictions from component models
nMarCov.verif_avg <- bind_rows(
  Mod1_nMarCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
    ) %>%  
    rownames_to_column(), 
  Mod2_nMarCov.verif %>% 
    dplyr::select(
      wt_arithPred,
      wt_logitPred
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
    zone = Mod1_nMarCov.verif$zone,
    spp = Mod1_nMarCov.verif$spp,
    arithObs = Mod1_nMarCov.verif$arithObs,
    logitObs = Mod1_nMarCov.verif$logitObs
  ) %>% 
  rename(
    arithPred = wt_arithPred,
    logitPred = wt_logitPred
  ) %>% 
  as.data.frame()


nMarCov.IOAarithFix <- modStats(
  nMarCov.verif_avg,
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