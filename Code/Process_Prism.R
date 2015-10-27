library(devtools)
install_github(repo = "prism", username = "ropensci")
library(prism)

options(prism.path = "prismtmp")
get_prism_monthlys(type="tmean", year = 1980:2014, mon = 7, keepZip=F)
ls_prism_data()

# use regular expressions to grep through the list and get data only from one month

loc <- round(c(df[which(df$featureid == featureids[1]), ]$NodeLon, df[which(df$featureid == featureids[1]), ]$NodeLat), digits = 4)

to_slice <- grep("_[0-9]{4}[0][7]",ls_prism_data()[,1],value=T)
to_slice = grep("tmean",to_slice, value = T)
p <- prism_slice(loc,to_slice)
p + stat_smooth(method="lm",se=F) + theme_bw() + ggtitle("Average July temperature in Boulder, CO 1980-2014")


foo <- extract_prism(loc, to_slice)

extract_prism <- function (location, prismfile) 
{
  if (!is.null(dim(prismfile))) {
    stop("You must enter a vector of file names, not a data frame, try  ls_prism_data()[,1]")
  }
  meta_d <- unlist(prism_md(prismfile, returnDate = T))
  meta_names <- unlist(prism_md(prismfile))[1]
  param_name <- strsplit(meta_names, "-")[[1]][3]
  pstack <- prism_stack(prismfile)
  data <- unlist(raster::extract(pstack, matrix(location, nrow = 1), 
                         buffer = 10))
  data <- as.data.frame(data)
  data$date <- as.Date(meta_d)
  data <- data[order(data$date), ]
  data <- data %>%
    group_by(date) %>%
    dplyr::rename(value = data) %>%
    dplyr::summarise(value = mean(data))
  return(data)
}
