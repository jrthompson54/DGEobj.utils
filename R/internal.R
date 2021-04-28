.tsmsg <- function(...) {
    # works like message() but prepends a timestamp
    message(date(), ": ", ...)
}
