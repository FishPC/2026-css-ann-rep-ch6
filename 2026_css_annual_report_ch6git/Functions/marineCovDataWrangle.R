library(readr)
library(dplyr)
library(tidyr)
library(splitstackshape)
library(janitor)
library(lubridate)

# Read unstructured lines
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

pdo.dat <- read_lines("https://www.ncei.noaa.gov/pub/data/cmb/ersst/v5/v6/index/ersst.v6.pdo.dat") %>% 
  as.data.frame() %>% 
  tail(-1) %>% 
  rename(
    'data' = 1
  ) %>% 
  cSplit('data', sep=" ", type.convert=FALSE) %>% 
  row_to_names(row_number = 1) %>% 
  rename_with(tolower)

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
  summarize(
    mean_sst = mean(sst)
  ) %>% 
  rename(
    'sst.maysep' = 'mean_sst'
  )%>% 
  dplyr::select(
    year,
    sst.maysep
  )

combMarCovs <- upwel.dat %>% 
  full_join(
    npgo.dat,
    by = 'year'
  ) %>% 
  full_join(
    sst.dat,
    by = 'year'
  )
