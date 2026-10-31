# pull Caldas Novas out of the CHIRPS workbook into a small CSV
suppressMessages(library(readxl))

ct <- c("guess","numeric","skip","skip","skip","numeric","numeric",
        "skip","skip","skip","skip","numeric","skip","skip","skip","text","text")

d <- read_excel("data/rainfall.xlsx", sheet = "bra-rainfall-subnat-5ytd", col_types = ct)
cat("full data sheet rows:", nrow(d), "\n")

ca <- d[!is.na(d$municipality) & tolower(d$municipality) == "caldas novas", ]
cat("Caldas Novas rows:", nrow(ca), "\n")
cat("date range:", format(min(ca$date)), "->", format(max(ca$date)), "\n")
cat("states seen:", paste(unique(ca$state), collapse=", "), "\n")

write.csv(ca, "data/rainfall_caldas.csv", row.names = FALSE)
cat("Wrote rainfall_caldas.csv\n")
print(utils::head(as.data.frame(ca)))
print(utils::tail(as.data.frame(ca)))
