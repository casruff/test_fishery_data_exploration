

if(!require("readr")) {
  install.packages("readr")
  library("readr")
}
if(!require("here")) {
  install.packages("here")
  library("here")
}
## set data dir
datadir <- here("data")

## first & last years of fish data
yr_frst <- 2013
yr_last <- 2020


## data file names
## 1. file with escapement data
fn_trs <- "sthd_trs.csv"

## 2. file with age comp data
fn_tst_fsh <- "skgt_sthd_tst_fsh.csv"


## read in data
dat_trs <- read_csv(file.path(datadir,fn_trs))
dat_tst_fsh <- read_csv(file.path(datadir,fn_tst_fsh))

## get flow data (currently using mainstem at mount vernon)
## flow gage ID
flow_site <- 12200500  
## get URL for flow data from USGS
flow_url <- paste0("https://waterdata.usgs.gov/nwis/dv",
                   "?cb_00060=on",
                   "&format=rdb",
                   "&site_no=",flow_site,
                   "&begin_date=",yr_frst,"-01-01",
                   "&end_date=",yr_last,"-12-31")

## raw flow data from USGS
flow_raw <- read_lines(flow_url)
## lines with metadata
hdr_flow <- which(lapply(flow_raw, grep, pattern = "\\#")==1, arr.ind = TRUE)
## print flow metadata
print(flow_raw[hdr_flow], quote = FALSE)

## flow data for years of interest
dat_flow <-  read_tsv(flow_url,
                      col_names = FALSE,
                      col_types = "ciDdc",
                      skip = max(hdr_flow)+2)
colnames(dat_flow) <- unlist(strsplit(tolower(flow_raw[max(hdr_flow)+1]),
                                      split = "\\s+"))
head(dat_flow)

## keep only relevant columns
dat_flow <- dat_flow[c("datetime", grep("[0-9]$", colnames(dat_flow), value = TRUE))]
## nicer column names
colnames(dat_flow) <- c("date","flow")
## convert cubic feet to cubic meters
dat_flow$flow <- dat_flow$flow / 35.3147
## flow by year & month
dat_flow$year <- as.integer(format(dat_flow$date,"%Y"))
dat_flow$day <- as.integer(format(dat_flow$date,"%j"))
dat_flow$month <- as.integer(format(dat_flow$date,"%m"))
dat_flow <- dat_flow[,c("year","month","flow")]

dat_tst_fsh


