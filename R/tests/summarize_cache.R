library(xts)
library(lubridate)

out_lines <- c()
files <- list.files('data/cache', pattern='\\.rds$', full.names=TRUE)
if (length(files) == 0) {
  message('No cache files found in data/cache. Skipping cache summarization test.')
  return(invisible(NULL))
}

for (f in files) {
  s <- tryCatch(readRDS(f), error = function(e) e)
  if (inherits(s, 'error')) {
    out_lines <- c(out_lines, paste('---', f, '---', 'READ ERROR:', s$message))
    next
  }

  # support old plain xts or list(data=xts, meta=list)
  if (is.list(s) && !is.null(s$data)) {
    obj <- s$data
    meta <- s$meta
  } else {
    obj <- s
    meta <- NULL
  }

  cls <- class(obj)[1]
  if (inherits(obj, 'xts') || inherits(obj, 'zoo')) {
    idx <- index(obj)
    start_date <- tryCatch(as.character(min(as.Date(idx))), error = function(e) NA)
    end_date <- tryCatch(as.character(max(as.Date(idx))), error = function(e) NA)
    n <- NROW(obj)
    na_count <- sum(is.na(obj))
    na_pct <- round(100 * na_count / max(n,1), 4)
    periodicity <- tryCatch(periodicity(obj)$scale, error = function(e) NA)
    # check monthly conventions
    is_all_first <- all(format(as.Date(idx), '%d') == '01')
    is_all_monthend <- all(as.Date(idx) == (ceiling_date(as.Date(idx), 'month') - days(1)))

    values <- as.numeric(coredata(obj))
    stats <- c(mean = round(mean(values, na.rm=TRUE), 6), sd = round(sd(values, na.rm=TRUE), 6), median = round(median(values, na.rm=TRUE),6), min = round(min(values, na.rm=TRUE),6), max = round(max(values, na.rm=TRUE),6))

    out_lines <- c(out_lines, paste('---', f, '---'))
    out_lines <- c(out_lines, paste('Class:', cls))
    if (!is.null(meta)) out_lines <- c(out_lines, paste('Meta.source:', meta$source, 'meta.from:', meta$from, 'meta.to:', meta$to))
    out_lines <- c(out_lines, paste('Start:', start_date, 'End:', end_date, 'N:', n, 'NA_count:', na_count, 'NA%:', na_pct, 'Per:', periodicity))
    out_lines <- c(out_lines, paste('MonthlyIndex:first?', is_all_first, 'month_end?', is_all_monthend))
    out_lines <- c(out_lines, paste('Stats: mean=', stats['mean'], ' sd=', stats['sd'], ' median=', stats['median'], ' min=', stats['min'], ' max=', stats['max']))

  } else if (is.data.frame(obj) || is.matrix(obj)) {
    n <- nrow(obj)
    na_count <- sum(is.na(obj))
    na_pct <- round(100 * na_count / max(n,1), 4)
    out_lines <- c(out_lines, paste('---', f, '---'))
    out_lines <- c(out_lines, paste('Class:', class(obj)[1], 'Rows:', n, 'NA_count:', na_count, 'NA%:', na_pct))
    if (!is.null(meta)) out_lines <- c(out_lines, paste('Meta.source:', meta$source, 'meta.from:', meta$from, 'meta.to:', meta$to))
  } else {
    out_lines <- c(out_lines, paste('---', f, '---'))
    out_lines <- c(out_lines, paste('Class:', class(obj)[1], 'Unable to summarize automatically'))
  }
}

cat(paste(out_lines, collapse='\n'))
