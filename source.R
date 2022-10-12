## ---- echo=FALSE, warning=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(tidyverse, quietly = TRUE))
suppressPackageStartupMessages(library(magrittr, quietly = TRUE))
suppressPackageStartupMessages(library(lubridate, quietly = TRUE))
suppressPackageStartupMessages(library(broom, quietly = TRUE))
suppressPackageStartupMessages(library(modelr, quietly = TRUE))

assert <- function(expr, expected) {
	if (!expr) stop(paste0("ERROR: ", expr))
}

knit_hooks$set(plot = function(x, options) {
	base = sub("\\s+$", "", hook_plot_md(x, options))
    paste0(base, "{#fig:", options$label, "}")
})


options(width = 50)
Sys.setenv(LANG = "en")


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
my_data <- read_csv(file = "data/data.csv")


## -----------------------------------------------------------------------------
spec(my_data)


## -----------------------------------------------------------------------------
my_data <- read_csv(file = "data/data.csv",
                    show_col_types = FALSE)


## -----------------------------------------------------------------------------
read_csv(
    "A, B, C,    D
     1, a, a,  1.2
     2, b, b,  2.1
     3, c, c, 13.0
", show_col_types = FALSE)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_names = FALSE,
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data-no-header.csv",
    col_names = FALSE,
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data-no-header.csv",
    col_names = c("X", "Y", "Z", "W"),
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_names = c("X", "Y", "Z", "W"),
    skip = 1,
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    "A, B, C,    D  # this is a comment
    1, a, a,  1.2   # another comment
    2, b, b,  2.1
    3, c, c, 13.0",
    comment = "#",
    show_col_types = FALSE)


## ---- warning=TRUE------------------------------------------------------------
read_csv(
    "A, B, C,    D  # this is a comment
# whole line comment
    1, a, a,  1.2   # another comment
    2, b, b,  2.1
    3, c, c, 13.0",
    comment = "#",
    show_col_types = FALSE)


## -----------------------------------------------------------------------------
read_csv(
    "A, B, C,    D  # this is a comment
    # the indentation is a potential problem; missing columns?
    1, a, a,  1.2   # another comment
    2, b, b,  2.1
    3, c, c, 13.0",
    comment = "#",
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = "????"
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = "dccd"
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = "iccd"
)


## ---- warning=TRUE------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = "icci"
)


## -----------------------------------------------------------------------------
read_csv(
    'A, B, C,    D,      E
    $1,a,a,1.2%,"1,100,200"
    $2,b,b,2.1%,"140,000"
    $3,c,c,13.0%,"2,005,000"',
    col_types = "nccnn"
    )


## -----------------------------------------------------------------------------
read_csv(
    'A, B, C,    D,      E
    $1,a,a," 1,2%","1.100.200"
    $2,b,b," 2,1%","  140.000"
    $3,c,c,"13,0%","2.005.000"',
    locale=locale(decimal_mark = ",", grouping_mark = "."),
    col_types = "nccnn")


## -----------------------------------------------------------------------------
read_csv(
    'A, B, C,    D
    TRUE, a, a,  1.2
    false, b, b,  2.1
    true, c, c, 13',
    show_col_types = FALSE
)


## -----------------------------------------------------------------------------
read_csv(
    'A, B, C,    D
    1, a, a,  1.2
    0, b, b,  2.1
    1, c, c, 13',
    col_types = "lccn"
)


## -----------------------------------------------------------------------------
read_csv(
    'D, T, t
    "2018-08-23", "2018-08-23T14:30", 14:30',
    col_types = "DTt"
)


## -----------------------------------------------------------------------------
read_csv(
    'D, t
    "23 Oct 2018", 2pm',
    col_types = "Dt",
    locale = locale(
        date_format = "%d %b %Y",
        time_format = "%I%p"
    )
)


## -----------------------------------------------------------------------------
read_csv(
    'A, B, C,    D
    1, a, a,  1.2
    0, b, b,  2.1
    1, c, c, 13',
    col_types = "lcfn")


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_type = "_cc-"
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = cols(A = "c")
)
read_csv(
    file = "data/data.csv",
    col_types = cols(A = "c", D = "c")
)


## -----------------------------------------------------------------------------
read_csv(
    file = "data/data.csv",
    col_types = cols(A = col_integer())
)
read_csv(
    file = "data/data.csv",
    col_types = cols(D = col_character())
)


## -----------------------------------------------------------------------------
my_data <- read_csv(
    file = "data/data.csv",
    col_types = cols(C = col_factor())
)
my_data$C


## -----------------------------------------------------------------------------
my_data <- read_csv(
    "A, B, C
     Foo, 12.4, Medium
     Bar, 5.2,  Small
     Baz, 42.0, Large
    ",
    col_types = cols(
        C = col_factor(levels = c("Small", "Medium", "Large"))
    )
)
my_data$C


## -----------------------------------------------------------------------------
my_data <- read_csv(
    file = "data/data.csv",
    col_types = cols(
        B = col_factor(ordered = TRUE),
        C = col_factor(levels = c("c", "b", "a"))
    )
)
my_data$B
my_data$C


## -----------------------------------------------------------------------------
parse_time(c("18:00"), format = "%R")
parse_time(c("18:00:30"), format = "%T")


## -----------------------------------------------------------------------------
parse_time(c("6 pm"),  format = "%I %p")


## -----------------------------------------------------------------------------
parse_time(c("18"), format = "%R")


## -----------------------------------------------------------------------------
parse_date(c("1975-02-05"), format = "%Y-%m-%d")


## -----------------------------------------------------------------------------
parse_date(c("75-02"), format = "%y-%m")


## -----------------------------------------------------------------------------
parse_date(c("15/02/75"), format = "%d/%m/%y")
parse_date(c("02/15/75"), format = "%m/%d/%y")


## -----------------------------------------------------------------------------
parse_date(c("Feb 15 1975"), format = "%b %d %Y", locale = locale("en"))
parse_date(c("15. feb. 1975"), format = "%d. %b %Y", locale = locale("da"))
parse_date(c("February 15 1975"), format = "%B %d %Y", locale = locale("en"))
parse_date(c("15. feb. 1975"), format = "%d. %b %Y", locale = locale("da"))
parse_date(c("Oct 15 1975"), format = "%b %d %Y", locale = locale("en"))
parse_date(c("15. okt. 1975"), format = "%d. %b %Y", locale = locale("da"))
parse_date(c("October 15 1975"), format = "%B %d %Y", locale = locale("en"))
parse_date(c("15. oktober 1975"), format = "%d. %B %Y", locale = locale("da"))


## -----------------------------------------------------------------------------
parse_datetime(c("Feb 15 1975 18:00 US/Pacific"), format = "%b %d %Y %R %Z")
parse_datetime(c("Feb 15 1975 18:00 -0800"), format = "%b %d %Y %R %z")


## -----------------------------------------------------------------------------
parse_datetime(c("Jun 15 1975 18:00 US/Pacific"), format = "%b %d %Y %R %Z")
parse_datetime(c("Jun 15 1975 18:00 -0700"), format = "%b %d %Y %R %z")


## -----------------------------------------------------------------------------
parse_datetime(
  c("Aug 15 1975 18:00"),
  format = "%b %d %Y %R",
  locale = locale(tz = "US/Pacific")
)
parse_datetime(
  c("Aug 15 1975 18:00 US/Pacific"),
  format = "%b %d %Y %R %Z"
)


## -----------------------------------------------------------------------------
x <- parse_datetime(
  c("Aug 15 1975 18:00"),
  format = "%b %d %Y %R",
  locale = locale(tz = "US/Pacific")
)
y <- parse_datetime(
  c("Aug 15 1975 18:00 US/Pacific"),
  format = "%b %d %Y %R %Z"
)
x == y


## -----------------------------------------------------------------------------
read_table(
    "A  B  C  D
     1  2  3  4
    15 16 17 18"
)


## -----------------------------------------------------------------------------
read_table(
    "A  B  C  D
     1  2  3  4
    15 16 17 18",
    col_names = FALSE
)


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
x <- read_csv(file = "data/data.csv", show_col_types = FALSE)
x


## -----------------------------------------------------------------------------
class(x)


## -----------------------------------------------------------------------------
y <- as.data.frame(x)
y
class(y)


## -----------------------------------------------------------------------------
z <- as_tibble(y)
z


## -----------------------------------------------------------------------------
x <- tibble(
    x = 1:100,
    y = x^2,
    z = y^2
)
x


## -----------------------------------------------------------------------------
print(x, n = 2)


## -----------------------------------------------------------------------------
print(x, n = 2, width = 15)


## -----------------------------------------------------------------------------
tribble(
    ~x, ~y,  ~z,
     1, 10, 100,
     2, 20, 200,
     3, 30, 300
)


## -----------------------------------------------------------------------------
x <- read_csv(file = "data/data.csv", show_col_types = FALSE)
y <- as.data.frame(x)
x["A"]
y["A"]
x[1]
y[1]


## -----------------------------------------------------------------------------
x[["A"]]
y[["A"]]


## -----------------------------------------------------------------------------
x$A
y$A


## -----------------------------------------------------------------------------
x[c("A", "C")]
y[c("A", "C")]
x[1:2]
y[1:2]


## -----------------------------------------------------------------------------
x[1:2,1:2]
y[1:2,1:2]


## -----------------------------------------------------------------------------
x[1:2,]
y[1:2,]


## -----------------------------------------------------------------------------
x[1:2,2]
y[1:2,2]


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
tbl <- tribble(
  ~sample, ~min_size, ~max_size, ~min_weight, ~max_weight,
  "foo", 13, 16, 45.2, 67.2,
  "bar", 12, 17, 83.1, 102.5
)
tbl |> select(sample, min_size, min_weight)


## -----------------------------------------------------------------------------
tbl |> select(min_size : max_weight)


## -----------------------------------------------------------------------------
tbl |> select(sample : min_size, min_weight : max_weight)


## -----------------------------------------------------------------------------
tbl |> select(!sample)


## -----------------------------------------------------------------------------
tbl |> select(!(sample:min_size))


## -----------------------------------------------------------------------------
tbl |> select(sample, !(sample:min_size))


## -----------------------------------------------------------------------------
tbl |> select( 
  sample:min_weight   # sample, min_size, max_size, and min_weight
  &                   # intersect with
  max_size:max_weight # max_size, min_weight, max_weight
)


## -----------------------------------------------------------------------------
tbl |> select( 
  sample:min_weight   # sample, min_size, max_size, and min_weight
  |                   # union with
  max_size:max_weight # max_size, min_weight, max_weight
)


## -----------------------------------------------------------------------------
tbl |> select(starts_with("min"))


## -----------------------------------------------------------------------------
tbl |> select(ends_with("weight"))


## -----------------------------------------------------------------------------
tbl |> select(contains("_"))


## -----------------------------------------------------------------------------
tbl |> select(matches(".*_.*"))


## -----------------------------------------------------------------------------
tbl |> select(everything())


## -----------------------------------------------------------------------------
tbl |> select(last_col())


## -----------------------------------------------------------------------------
tbl |> select(last_col(0))
tbl |> select(last_col(1))
tbl |> select(last_col(2))
tbl |> select(last_col(3))


## -----------------------------------------------------------------------------
tbl |> select(last_col(3):last_col())


## -----------------------------------------------------------------------------
vars <- c("min_size", "min_weight")


## -----------------------------------------------------------------------------
tbl |> select(all_of(vars))
tbl |> select(any_of(vars))


## ---- error=TRUE--------------------------------------------------------------
vars <- c(vars, "foo")
tbl |> select(all_of(vars))
tbl |> select(any_of(vars))


## -----------------------------------------------------------------------------
tbl |> select(where(is.numeric))


## -----------------------------------------------------------------------------
tbl |> select(where(\(x) is.numeric(x) && max(x) > 100.0))


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
mean_income <- tribble(
    ~country,  ~`2001`, ~`2002`, ~`2003`, ~`2004`, ~`2005`,
    "Numenor",  123456,  132654,  321646,  324156,  325416,
    "Westeros", 314256,  432165,  546123,  465321,  561423,
    "Narnia",   432156,  342165,  564123,  543216,  465321,
    "Gondor",   531426,  321465,  235461,  463521,  561423,
    "Laputa",    14235,   34125,   45123,   51234,   54321
)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    cols = c(`2001`, `2001`, `2002`,
             `2003`, `2004`, `2005`)
)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    cols = `2001`:`2005`
)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    cols = !country
)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    starts_with("2")
)
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    !starts_with("c")
)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    !country,
    names_transform = as.integer
)


## -----------------------------------------------------------------------------
tidy_income <- mean_income |> pivot_longer(
    names_to = "year",
    values_to = "income",
    !country
)

tidy_income |> pivot_wider(
    names_from = "year",
    values_from = "income"
)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~date,
    "11/5",
    "4/7",
    "21/12"
)
separate(tbl, date, into = c("day", "month"))


## -----------------------------------------------------------------------------
separate(tbl, date, into = c("day", "month"), remove = FALSE)


## -----------------------------------------------------------------------------
separate(tbl, date, into = c("day", "month"), convert = TRUE)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~date,
    "11th of month 5",
    "4th of month 7",
    "21st of month 12"
)
separate(tbl, date, into = c("day", "month"), 
         sep = "[[:alpha:][:space:]]+")


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~date,
    "11th of May",
    "4th of July",
    "21st of December"
)


## -----------------------------------------------------------------------------
tbl |> extract(
    col = date,
    into = c("day", "month"),
    regex = paste0(
        "([[:digit:]]+)",        # First group, the day
        "[[:alpha:][:space:]]+", # Stuff we ignore; will match "th of" and such
        "[[:space:]]",           # The final space before the month name
        "([[:alpha:]]+)"         # Second group, the month
    )
)


## -----------------------------------------------------------------------------
tbl |> extract(
    col = date,
    into = c("day", "fluf", "month"),
    regex = paste0(
        "([[:digit:]]+)",        # First group, the day
        "([[:alpha:]]+)",        # st, nd, rd, th, that kind of stuff goes in fluf
        "[[:alpha:][:space:]]+", # stuff between day and month; not a group so we throw it away
        "([[:alpha:]]+)"         # Third group, the month
    )
)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~date,
    "11/5",
    "4/7",
    "21/12"
)
tbl2 <- tbl |> separate(col = date, into = c("day", "month"))


## -----------------------------------------------------------------------------
tbl2 |> unite(col = "date", day, month)


## -----------------------------------------------------------------------------
tbl2 |> unite(col = "date", day, month, sep = "/")


## -----------------------------------------------------------------------------
tbl2 |> unite(
   col = "date", day, month, 
   sep = "/", remove = FALSE
)


## -----------------------------------------------------------------------------
military_casualties <- tribble(
    ~war, ~groups, ~deaths,
    
    'WW1',
    "Allied Powers/Central Powers",
    "5.7,4.0",
    
    'WW2', 
    "Germany/Japan/USSR/British Empire/USA", 
    "5.3,2.1,10.7,0.6,0.4"
)


## -----------------------------------------------------------------------------
military_casualties |> separate_rows(
    groups, deaths,
    sep = "/|,"
)


## -----------------------------------------------------------------------------
military_casualties |> separate_rows(
    groups, deaths,
    sep = "/|,",
    convert = TRUE
)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~A, ~B, ~C,
     1, 11, 21,
     2, 11, 22,
     4, 13, 32
)
tbl |> expand(A, B)


## -----------------------------------------------------------------------------
tbl |> expand(A = 1:4, B)


## -----------------------------------------------------------------------------
crossing(A = 1:3, B = 11:13)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~A, ~B, ~C,
    1, 11, 21,
    2, 11, 22,
    2, 11, 12,
    4, 13, 42,
    4, 13, 32
)
expand(tbl, nesting(A, B))


## -----------------------------------------------------------------------------
complete(tbl, A = 1:4)
complete(tbl, B = 11:13)
complete(tbl, A = 1:4, B = 11:13)


## -----------------------------------------------------------------------------
mean_income <- tribble(
    ~country,  ~`2002`, ~`2003`, ~`2004`, ~`2005`,
    "Numenor",  123456,  132654,      NA,  324156,
    "Westeros", 314256,  NA,          NA,  465321,
    "Narnia",   432156,  NA,          NA,      NA,
    "Gondor",   531426,  321465,  235461,  463521,
    "Laputa",    14235,   34125,   45123,   51234,
)
drop_na(mean_income)


## -----------------------------------------------------------------------------
mean_income |> pivot_longer(
    names_to = "year",
    values_to = "mean_income",
    -country
) |> drop_na()


## -----------------------------------------------------------------------------
replace_na(mean_income, list(
    `2003` = mean(mean_income$`2003`, na.rm = TRUE)
))

replace_na(mean_income, list(
    `2003` = mean(mean_income$`2003`, na.rm = TRUE),
    `2004` = mean(mean_income$`2004`, na.rm = TRUE),
    `2005` = mean(mean_income$`2005`, na.rm = TRUE)
))


## -----------------------------------------------------------------------------
tbl <- read_csv(
    "year, quarter, income
    2011, Q1, 13
    , Q2, 12
    , Q3, 14
    , Q4, 11
    2012, Q1, 12
    , Q2, 14
    , Q3, 15
    , Q4, 17",
    show_col_types = FALSE
)
spec(tbl)


## -----------------------------------------------------------------------------
fill(tbl, year)


## -----------------------------------------------------------------------------
tbl <- tribble(
    ~A, ~B, ~C,
    1, 11, 21,
    1, 11, 12,
    2, 42, 22,
    2, 15, 22,
    2, 15, 32,
    4, 13, 32
)


## -----------------------------------------------------------------------------
nested_tbl <- tbl |> nest(BC = c(B, C))
nested_tbl
nested_tbl[[1,2]]


## -----------------------------------------------------------------------------
nested_tbl |> unnest(cols = c(BC))


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
write_csv(
    tribble(
        ~country,  ~`2002`, ~`2003`, ~`2004`, ~`2005`,
        "Numenor",  123456,  132654,      NA,  324156,
        "Westeros", 314256,  NA,          NA,  465321,
        "Narnia",   432156,  NA,          NA,      NA,
        "Gondor",   531426,  321465,  235461,  463521,
        "Laputa",    14235,   34125,   45123,   51234,
    ),
    "data/income.csv"
)


## -----------------------------------------------------------------------------
mydata <- read_csv("data/income.csv",
                   col_types = "cdddd")
mydata <- pivot_longer(
    data = mydata,
    names_to = "year",
    values_to = "mean_income",
    !country
)
mydata <- drop_na(mydata)


## -----------------------------------------------------------------------------
mydata <- drop_na(
    pivot_longer(
        data = read_csv("data/income.csv",
                        col_types = "cdddd"),
        names_to = "year",
        values_to = "mean_income",
        !country
    )   
)


## ---- echo = FALSE------------------------------------------------------------
f <- invisible
x <- 5


## -----------------------------------------------------------------------------
x %>% f()


## -----------------------------------------------------------------------------
x |> f()


## -----------------------------------------------------------------------------
f(x)


## ---- echo = FALSE------------------------------------------------------------
f <- g <- identity
h <- invisible


## -----------------------------------------------------------------------------
x %>% f() %>% g() %>% h()


## -----------------------------------------------------------------------------
x |> f() |> g() |> h()


## -----------------------------------------------------------------------------
h(g(f(x)))


## -----------------------------------------------------------------------------
x %>% f %>% g %>% h


## -----------------------------------------------------------------------------
"hello," %>% (\(y){ paste(y, "world!")})


## -----------------------------------------------------------------------------
"hello," |> (\(y){ paste(y, "world!")})()


## -----------------------------------------------------------------------------
mydata <- read_csv("data/income.csv", col_types = "cdddd") %>%
    pivot_longer(names_to = "year", values_to = "mean_income", !country) %>%
    drop_na()


## -----------------------------------------------------------------------------
mydata <- read_csv("data/income.csv",
                   col_types = "cdddd") |>
    pivot_longer(names_to = "year", values_to = "mean_income", !country) |>
    drop_na()

mydata <- read_csv("data/income.csv",
                   col_types = "cdddd") %>%
    pivot_longer(names_to = "year", values_to = "mean_income", !country) %>%
    drop_na()

mydata <- read_csv("data/income.csv", 
                   col_types = "cdddd") %>%
    pivot_longer(data = ., names_to = "year", values_to = "mean_income", !country) %>%
    drop_na()

mydata <- read_csv("data/income.csv", 
                   col_types = "cdddd") %>%
    pivot_longer(names_to = "year", values_to = "mean_income", !country, data = .) %>%
    drop_na()


## -----------------------------------------------------------------------------
rnorm(5) %>% tibble(x = ., y = .)


## -----------------------------------------------------------------------------
rnorm(5) %>% tibble(x = ., y = abs(.))


## -----------------------------------------------------------------------------
rnorm(5) %>% tibble(x = sin(.), y = abs(.))


## -----------------------------------------------------------------------------
rnorm(5) %>% { tibble(x = sin(.), y = abs(.)) }


## -----------------------------------------------------------------------------
h <- . %>% f() %>% g()


## -----------------------------------------------------------------------------
pipeline <- . %>% 
    pivot_longer(names_to = "year", values_to = "mean_income", !country) %>%
    drop_na()
pipeline


## ---- include=FALSE-----------------------------------------------------------
suppressPackageStartupMessages(library(magrittr, quietly = TRUE))


## -----------------------------------------------------------------------------
mydata <- read_csv("data/income.csv", col_types = "cdddd")
mydata %<>% pipeline()


## -----------------------------------------------------------------------------
mydata <- read_csv("data/income.csv", col_types = "cdddd")
mydata <- mydata %>% pipeline()


## -----------------------------------------------------------------------------
mydata <- tibble(x = rnorm(5), y = rnorm(5))
mydata %>% { .$x - .$y }


## -----------------------------------------------------------------------------
mydata %$% { x - y }


## ---- eval=FALSE--------------------------------------------------------------
## tidy_income <- read_csv("data/income.csv", col_types = "cdddd") |>
##     pivot_longer(names_to = "year", values_to = "mean_income", !country)
## 
## # Summarize and visualize data before we continue
## summary(tidy_income)
## ggplot(tidy_income, aes(x = year, y = mean_income)) + geom_point()
## 
## # Continue processing
## tidy_income |>
##     drop_na() |>
##     write_csv("data/tidy-income.csv")


## ---- eval=FALSE--------------------------------------------------------------
## mydata <- read_csv("data/income.csv", col_types = "cdddd") %>%
##     pivot_longer(names_to = "year", values_to = "mean_income", !country) %T>%
##     # Summarize and then continue
##     { print(summary(.)) } %T>%
##     # Plot and then continue
##     { print(ggplot(., aes(x = year, y = mean_income)) + geom_point()) } %>%
##     drop_na() %>% write_csv("data/tidy-income.csv")


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
is_even <- function(x) x %% 2 == 0
1:6 |> keep(is_even)
1:6 |> discard(is_even)


## -----------------------------------------------------------------------------
1:6 |> keep(negate(is_even))
1:6 |> discard(negate(is_even))


## -----------------------------------------------------------------------------
y <- list(NULL, 1:3, NULL)
y |> compact()


## -----------------------------------------------------------------------------
x <- y <- 1:3
names(y) <- c("one", "two", "three")


## -----------------------------------------------------------------------------
names(x)
names(y)


## -----------------------------------------------------------------------------
z <- list(x = x, y = y)
z |> compact(names)


## -----------------------------------------------------------------------------
is_even <- function(x) x %% 2 == 0
1:4 |> map(is_even)


## -----------------------------------------------------------------------------
1:4 |> map_lgl(is_even)


## -----------------------------------------------------------------------------
1:4 %% 2 == 0


## -----------------------------------------------------------------------------
n <- seq.int(100, 1000, 300) # the different n value we want
n

# Are we sampling for each of our n values?
rnorm(n) # no, but length(n) is used for the number of samples in rnorm()


## -----------------------------------------------------------------------------
sem <- function(n) sd(rnorm(n = n)) / sqrt(n)
n |> map_dbl(sem)


## -----------------------------------------------------------------------------
1:3 %>% map_dbl(identity) %T>% print() %>% class()
1:3 %>% map_chr(identity) %T>% print() %>% class()
1:3 %>% map_int(identity) %T>% print() %>% class()


## -----------------------------------------------------------------------------
x <- tibble(a = 1:2, b = 3:4)
list(a = x, b = x) |> map_dfr(identity)
list(a = x, b = x) |> map_dfc(identity)


## -----------------------------------------------------------------------------
list(a = x, b = x) |> map_df(identity)


## -----------------------------------------------------------------------------
x <- list(1:3, 4:6)
x |> map_dbl(1)
x |> map_dbl(3)


## -----------------------------------------------------------------------------
x <- list(
    c(a = 42, b = 13),
    c(a = 24, b = 31)
)
x |> map_dbl("a")
x |> map_dbl("b")


## -----------------------------------------------------------------------------
a <- tibble(foo = 1:3, bar = 11:13)
b <- tibble(foo = 4:6, bar = 14:16)
ab <- list(a = a, b = b)
ab |> map("foo")


## -----------------------------------------------------------------------------
ab |> map_depth(0, length)
ab |> length()


## -----------------------------------------------------------------------------
ab |> map_depth(1, sum) |> unlist()
ab |> map_depth(2, sum) |> unlist()
ab |> map_depth(3, sum) |> unlist()


## -----------------------------------------------------------------------------
is_even <- function(x) x %% 2 == 0
add_one <- function(x) x + 1
1:6 |> map_if(is_even, add_one) |> as.numeric()


## -----------------------------------------------------------------------------
1:6 |> keep(is_even) |> map_dbl(add_one)


## -----------------------------------------------------------------------------
add_two <- function(x) x + 2
1:6 |>
    map_if(is_even, add_one, .else = add_two) |>
    as.numeric()


## -----------------------------------------------------------------------------
1:6 |> map_at(2:5, add_one) |> as.numeric()


## -----------------------------------------------------------------------------
list(a = 1:3, b = 4:6) |> map(print) |> invisible()
list(a = 1:3, b = 4:6) |> lmap(print) |> invisible()


## -----------------------------------------------------------------------------
f <- function(x) list("foo")
1:2 |> lmap(f)
f <- function(x) list("foo", "bar")
1:2 |> lmap(f)


## -----------------------------------------------------------------------------
list(a = 1:3, b = 4:8) |> map(length) |> unlist()


## -----------------------------------------------------------------------------
wrapped_length <- function(x) {
    x |> length() |>         # get the length of x (result will be numeric)
         list()   |>         # go back to a list
         set_names(names(x)) # and give it the original name
}
list(a = 1:3, b = 4:8) |> lmap(wrapped_length) |> unlist()


## -----------------------------------------------------------------------------
wrapped_length <- function(x) {
    x |> pluck(1) |>         # pluck(x, 1) from purrr does that same as x[[1]]
         length() |>         # now we have the underlying vector and get the length
         list()   |>         # go back to a list
         set_names(names(x)) # and give it the original name
}
list(a = 1:3, b = 4:8) %>% lmap(wrapped_length) %>% unlist()


## -----------------------------------------------------------------------------
1:3 |> map(print) |> invisible()
1:3 |> walk(print)


## -----------------------------------------------------------------------------
x <- 1:3
y <- 3:1
map2_dbl(x, y, `+`)


## -----------------------------------------------------------------------------
list(x, y) |> pmap_dbl(`+`)


## -----------------------------------------------------------------------------
z <- 4:6
f <- function(x, y, z) x + y - z
list(x, y, z) |> pmap_dbl(f)


## -----------------------------------------------------------------------------
x <- c("foo", "bar", "baz")
f <- function(x, i) paste0(i, ": ", x)
x |> imap_chr(f)


## -----------------------------------------------------------------------------
modify2(1:3, 3:1, `+`)

x <- c("foo", "bar", "baz")
f <- function(x, i) paste0(i, ": ", x)
x |> imodify(f)


## -----------------------------------------------------------------------------
pair <- function(first, second) {
    structure(list(first = first, second = second),
              class = "pair")
}
toString.pair <- function(x, ...) {
    first <- toString(x$first, ...)
    rest <- toString(x$second, ...)
    paste('[', first, ', ', rest, ']', sep = '')
}
print.pair <- function(x, ...) {
    x |> toString() |> cat() |> invisible()
}


## -----------------------------------------------------------------------------
1:4 |> reduce(pair)


## -----------------------------------------------------------------------------
1:4 |> rev() |> reduce(pair)


## -----------------------------------------------------------------------------
1:4 |> reduce(pair, .dir = "backward")


## -----------------------------------------------------------------------------
1:3 |> reduce(pair, .init = 0)
1:3 |> rev() |> reduce(pair, .init = 4)
1:3 |> reduce(pair, .init = 4, .dir = "backward")


## -----------------------------------------------------------------------------
# additional arguments
loud_pair <- function(acc, next_val, volume) {
    # Build a pair
    ret <- pair(acc, next_val)
    # Announce that pair to the world
    ret |> toString() |>
        paste(volume, '\n', sep = '') |>
        cat()
    # Then return the new pair
    ret
}


## -----------------------------------------------------------------------------
1:3 |>
    reduce(loud_pair, volume = '!') |>
    invisible()
1:3 |>
    reduce(loud_pair, volume = '!!') |>
    invisible()


## -----------------------------------------------------------------------------
volumes <- c('!', '!!')
1:3 |> reduce2(volumes, loud_pair) |> invisible()
1:3 |>
    reduce2(c('!', '!!', '!!!'), .init = 0, loud_pair) |>
    invisible()


## -----------------------------------------------------------------------------
res <- 1:3 |> accumulate(pair)
print(res[[1]])
print(res[[2]])
print(res[[3]])

res <- 1:3 |> accumulate(pair, .init = 0)
print(res[[1]])
print(res[[4]])

res <- 1:3 |> accumulate(
    pair, .init = 0,
    .dir = "backward"
)
print(res[[1]])
print(res[[4]])


## -----------------------------------------------------------------------------
greater_than_three <- function(x) 3 < x
less_than_three <- function(x) x < 3

1:6 |> keep(greater_than_three)
1:6 |> keep(less_than_three)


## -----------------------------------------------------------------------------
1:6 |> keep(partial(`<`, 3))


## -----------------------------------------------------------------------------
`<`


## -----------------------------------------------------------------------------
1:6 |> keep(partial(`<`, e2 = 3))


## -----------------------------------------------------------------------------
1:6 |> map_dbl(partial(`+`, 2))
1:6 |> map_dbl(partial(`-`, 1))
1:3 |> map_dbl(partial(`-`, e1 = 4))
1:3 |> map_dbl(partial(`-`, e2 = 4))


## -----------------------------------------------------------------------------
1:3 |>
    map_dbl(partial(`+`, 2)) |>
    map_dbl(partial(`*`, 3))


## -----------------------------------------------------------------------------
1:3 |> map_dbl(
    compose(partial(`*`, 3), partial(`+`, 2))
)


## -----------------------------------------------------------------------------
is_even <- function(x) x %% 2 == 0
1:6 |> keep(is_even)


## -----------------------------------------------------------------------------
1:6 |> keep(function(x) x %% 2 == 0)


## -----------------------------------------------------------------------------
1:6 |> keep(\(x) x %% 2 == 0)


## -----------------------------------------------------------------------------
is_even_lambda <- ~ .x %% 2 == 0
1:6 |> keep(is_even_lambda)


## -----------------------------------------------------------------------------
1:6 |> keep(~ .x %% 2 == 0)


## ---- error=TRUE--------------------------------------------------------------
f <- function(x) 2 * x
f <- 5
f(2)


## -----------------------------------------------------------------------------
f <- function(x) 2 * x
g <- function() {
    f <- 5 # not a function
    f(2)   # will look for a function
}
g()


## -----------------------------------------------------------------------------
1:4 |> map_dbl(~ .x / 2)
1:3 |> map_dbl(~ 2 + .x)
1:3 |> map_dbl(~ 4 - .x)
1:3 |> map_dbl(~ .x - 4)


## -----------------------------------------------------------------------------
1:3 |> map_dbl(~ 3 * (.x + 2))


## -----------------------------------------------------------------------------
map2_dbl(1:3, 1:3, ~ .x + .y)


## -----------------------------------------------------------------------------
list(1:3, 1:3, 1:3) %>% pmap_dbl(~ .1 + .2 + .3)


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
iris_df <- as_tibble(iris)
print(iris_df, n = 3)
head(iris_df$Species)


## -----------------------------------------------------------------------------
iris_df |>
    select(Sepal.Length, Species) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(!Species) |>
    print(n = 3)

iris_df |>
    select(!c(Species, Sepal.Length)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(-Species) |>
    print(n = 3)

iris_df |>
    select(-c(Species, Sepal.Length)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(!Species, !Sepal.Length) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(-Species, -Sepal.Length) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(1) |>
    print(n = 3)

iris_df |>
    select(1:3) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(Petal.Length:Species) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(starts_with("Petal")) |>
    print(n = 3)

iris_df |>
    select(-starts_with("Petal")) |>
    print(n = 3)

iris_df |>
    select(starts_with("Petal"), Species) |>
    print(n = 3)

iris_df |>
    select(starts_with("PETAL", ignore.case = TRUE)) |>
    print(n = 3)

iris_df |>
    select(starts_with("S")) |>
    print(n = 3)

iris_df |>
    select(ends_with("Length")) |>
    print(n = 3)

iris_df |>
    select(contains("ng")) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(sepal_length = Sepal.Length,
           sepal_width = Sepal.Width) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    rename(sepal_length = Sepal.Length,
           sepal_width = Sepal.Width) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |> distinct(Species)


## -----------------------------------------------------------------------------
iris_df |>
    filter(Species == "setosa") |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter(Species == "setosa") |>
    select(ends_with("Length"), Species) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter(Species != "setosa") |>
    distinct(Species) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter(Sepal.Length > 5, Petal.Width < 0.4) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter(between(Sepal.Width, 2, 2.5)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter(str_starts(Species, "v")) |>
    print(n = 3)

iris_df |>
    filter(str_ends(Species, "r")) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    select(-Species) |>
    filter_all(any_vars(. > 5)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter_at(vars(-Species), any_vars(. > 5)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter_at(c("Petal.Length", "Sepal.Length"),
              any_vars(. > 0)) |>
    print(n = 3)

iris_df |>
    filter_at(vars(Petal.Length, Sepal.Length),
              any_vars(. > 0)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
iris_df |>
    filter_if(is.numeric, all_vars(. < 5)) |>
    print(n = 3)


## -----------------------------------------------------------------------------
df <- tribble(
    ~A, ~B, ~C,
     1,  2,  3,
     4,  5, NA,
    11, 12, 13,
    22, 22,  1
)
df |> filter_all(all_vars(. > 3))


## -----------------------------------------------------------------------------
df |> filter_if(~ all(!is.na(.)), all_vars(. > 3))


## -----------------------------------------------------------------------------
df |> filter_all(all_vars(is.na(.) | . > 3))


## -----------------------------------------------------------------------------
iris_df |>
    arrange(Petal.Length) |>
    print(n = 5)

iris_df |>
    arrange(Sepal.Length, Petal.Length) |>
    print(n = 5)


## -----------------------------------------------------------------------------
iris_df |>
    arrange(desc(Petal.Length)) |>
    print(n = 5)


## -----------------------------------------------------------------------------
df <- tribble(
    ~height, ~width,
         10,     12,
         42,     24,
         14,     12
)


## -----------------------------------------------------------------------------
df |> mutate(area = height * width)


## -----------------------------------------------------------------------------
cm_per_inch <- 2.54
df |> mutate(
    height_cm = cm_per_inch * height,
    width_cm  = cm_per_inch * width,
    area_cm   = height_cm * width_cm
)


## ---- error=TRUE--------------------------------------------------------------
df %>% mutate(
    area_cm   = height_cm * width_cm,
    height_cm = cm_per_inch * height,
    width_cm  = cm_per_inch * width
)


## -----------------------------------------------------------------------------
df |> mutate(cm_per_inch * height)


## ---- echo=FALSE, warning=FALSE-----------------------------------------------
suppressPackageStartupMessages(library(units, quietly = TRUE))


## -----------------------------------------------------------------------------
df |> mutate(
    height_in = units::as_units(height, "in"),
    width_in  = units::as_units(width, "in"),
    area_in   = height_in * width_in,
    height_cm = units::set_units(height_in, "cm"),
    width_cm  = units::set_units(width_in, "cm"),
    area_cm   = units::set_units(area_in, "cm^2")
)


## -----------------------------------------------------------------------------
df |> transmute(
    height_in = units::as_units(height, "in"),
    width_in  = units::as_units(width, "in"),
    area_in   = height_in * width_in,
    height_cm = units::set_units(height_in, "cm"),
    width_cm  = units::set_units(width_in, "cm"),
    area_cm   = units::set_units(area_in, "cm^2")
)


## ---- error=TRUE--------------------------------------------------------------
df <- tibble(
    x = rnorm(3, mean = 12, sd = 5),
)
df |> mutate(~ if (.x < 0) -.x else .x)


## -----------------------------------------------------------------------------
df |> mutate(abs(x))


## -----------------------------------------------------------------------------
df |> mutate(ifelse(x < 0, -x, x))


## -----------------------------------------------------------------------------
my_abs <- Vectorize(\(x) if (x < 0) -x else x)
df |> mutate(my_abs(x))


## -----------------------------------------------------------------------------
df <- tibble(x = rnorm(100))
df |>
    mutate(
        x_category = case_when(
            x - mean(x) < -2 * sd(x) ~ "small",
            x - mean(x) >  2 * sd(x) ~ "large",
            TRUE                     ~ "medium"
        )
    ) |>
    print(n = 3)


## -----------------------------------------------------------------------------
df <- tibble(x = rnorm(100), y = rnorm(100))
df |> summarise(mean_x = mean(x), mean_y = mean(y))


## -----------------------------------------------------------------------------
classify <- function(x) {
    case_when(
        x - mean(x) < -2 * sd(x) ~ "small",
        x - mean(x) >  2 * sd(x) ~ "large",
        TRUE           ~ "medium"
    )
}
df |> mutate(x_category = classify(x)) |>
      group_by(x_category) |>
      print(n = 3)


## -----------------------------------------------------------------------------
df |> mutate(x_category = classify(x)) |>
      group_by(x_category) |>
      summarise(mean_x = mean(x), no_x = n())


## -----------------------------------------------------------------------------
df |> mutate(x_category = classify(x)) |>
      group_by(x_category) |>
      group_vars()


## -----------------------------------------------------------------------------
df <- tibble(x = rnorm(100), y = rnorm(100)) |>
      mutate(
        x_category = classify(x),
        y_category = classify(y)
      )

df |> group_by(x_category, y_category) |>
      group_vars()


## -----------------------------------------------------------------------------
df |> group_by(x_category, y_category) |>
      summarise(mean_x = mean(x), mean_y = mean(y))


## -----------------------------------------------------------------------------
# Adding two categories at the same time
df |> group_by(x_category, y_category) |> group_vars()
# Adding categories in two steps
df |> group_by(x_category) |> group_by(y_category, .add = TRUE) |> group_vars()
# Replacing one grouping with another
df |> group_by(x_category) |> group_by(y_category) |> group_vars()


## -----------------------------------------------------------------------------
# Drop the last group we made
df |> group_by(x_category, y_category) |>
      summarise(mean_x = mean(x), mean_y = mean(y), .groups = "drop_last") |>
      group_vars()

# Drop all groups
df |> group_by(x_category, y_category) |>
      summarise(mean_x = mean(x), mean_y = mean(y), .groups = "drop") |>
      group_vars()

# Keep all groups
df |> group_by(x_category, y_category) |>
      summarise(mean_x = mean(x), mean_y = mean(y), .groups = "keep") |>
      group_vars()


## -----------------------------------------------------------------------------
df |> group_by(x_category, y_category) |>
      ungroup() |> # Remove the grouping again
      group_vars()


## -----------------------------------------------------------------------------
df |> group_by(x_category) |>
      mutate(mean_x = mean(x), mean_y = mean(y)) |>
      ungroup() |>
      print(n = 5)


## -----------------------------------------------------------------------------
df |> mutate(mean_x = mean(x), mean_y = mean(y)) |> # Calculate mean x and y for entire data
      distinct(mean_x, mean_y)                      #  Get unique values for printing


## -----------------------------------------------------------------------------
df |> group_by(x_category) |>                       # Group by x categories only
      mutate(mean_x = mean(x), mean_y = mean(y)) |> # Calculate mean x and y for each x group
      distinct(mean_x, mean_y)                      # Get unique values for printing


## -----------------------------------------------------------------------------
df |> mutate(mean_y = mean(y)) |> # Get the mean for y over all data points
      group_by(x_category)     |> # Then group so the following works per group
      mutate(mean_x = mean(x)) |> # Get the mean x for each group (not globally)
      distinct(mean_x, mean_y)    # Get the unique values so we can print the result


## -----------------------------------------------------------------------------
df |> group_by(x_category)     |>
      mutate(mean_x = mean(x)) |>
      group_by(y_category)     |>
      mutate(mean_y = mean(y)) |>
      distinct(
        x_category, mean_x,
        y_category, mean_y
      )


## -----------------------------------------------------------------------------
df1 <- tibble(
    A = paste0("a", 1:2),
    B = paste0("b", 1:2)
)
df2 <- tibble(
    A = paste0("a", 3:4),
    B = paste0("b", 3:4)
)
df3 <- tibble(
    C = paste0("c", 1:2),
    D = paste0("d", 1:2)
)
bind_rows(df1, df2)
bind_cols(df1, df3)


## -----------------------------------------------------------------------------
grades_maths <- tribble(
    ~name,            ~grade,
    "Marko Polo",     "D",
    "Isaac Newton",   "A+",
    "Charles Darwin", "B"
)
grades_biology <- tribble(
    ~name,            ~grade,
    "Marko Polo",     "F",
    "Isaac Newton",   "D",
    "Charles Darwin", "A+"
)

inner_join(grades_maths, grades_biology, by = "name")


## -----------------------------------------------------------------------------
grades_maths2 <- tribble(
    ~name,            ~grade,
    "Marko Polo",     "D",
    "Isaac Newton",   "A+", # so good at physics
    "Isaac Newton",   "A+", # that he got an A+ twice
    "Charles Darwin", "B"
)
grades_biology2 <- tribble(
    ~name,            ~grade,
    "Marko Polo",     "F",
    "Isaac Newton",   "D",
    "Charles Darwin", "A+", # so good at biology that we
    "Charles Darwin", "A+"  # listed him twice
)
inner_join(grades_maths2, grades_biology2, by = "name")
inner_join(grades_maths2, grades_biology2, by = "grade")


## -----------------------------------------------------------------------------
inner_join(grades_maths2, grades_biology2, by = "grade") |> distinct()


## -----------------------------------------------------------------------------
inner_join(
    grades_maths, grades_biology,
    by = "name", suffix = c(".maths", ".biology")
)


## -----------------------------------------------------------------------------
grades_geography <- tribble(
    ~name,            ~grade,
    "Marko Polo",     "A",
    "Charles Darwin", "A",
    "Immanuel Kant",  "A+"
)
grades_physics <- tribble(
    ~name,            ~grade,
    "Isaac Newton",    "A+",
    "Albert Einstein", "A+",
    "Charles Darwin",  "C"
)

inner_join(
    grades_geography, grades_physics,
    by = "name", suffix = c(".geography", ".physics")
)


## -----------------------------------------------------------------------------
full_join(
    grades_geography, grades_physics,
    by = "name", suffix = c(".geography", ".physics")
)


## -----------------------------------------------------------------------------
left_join(
    grades_geography, grades_physics,
    by = "name", suffix = c(".geography", ".physics")
)
right_join(
    grades_maths, grades_physics,
    by = "name", suffix = c(".maths", ".physics")
)


## -----------------------------------------------------------------------------
semi_join(
    grades_maths2, grades_biology2,
    by = "name", suffix = c(".geography", ".physics")
)


## -----------------------------------------------------------------------------
inner_join(
    grades_maths2, grades_biology2,
    by = "name", suffix = c(".geography", ".physics")
) |> select(1:2)


## -----------------------------------------------------------------------------
anti_join(
    grades_maths2, grades_physics,
    by = "name", suffix = c(".geography", ".physics")
)


## -----------------------------------------------------------------------------
grades <- list(
    grades_maths, grades_biology,
    grades_geography, grades_physics
)
grades |>
    reduce(full_join, by = "name") |>
    rename_at(2:5, ~ c("maths", "biology", "geography", "physics"))


## -----------------------------------------------------------------------------
mean_income <- tribble(
    ~country,  ~`2002`, ~`2003`, ~`2004`, ~`2005`,
    "Numenor",  123456,  132654,      NA,  324156,
    "Westeros", 314256,  NA,          NA,  465321,
    "Narnia",   432156,  NA,          NA,      NA,
    "Gondor",   531426,  321465,  235461,  463521,
    "Laputa",    14235,   34125,   45123,   51234,
)


## -----------------------------------------------------------------------------
mean_income                                                      |>
    pivot_longer(
        names_to = "year",
        values_to = "mean_income",
        !country
    )                                                            |> 
    group_by(
        country
    )                                                            |> 
    mutate(
        mean_per_country = mean(mean_income, na.rm = TRUE),
        mean_income = ifelse(
            is.na(mean_income),
            mean_per_country,
            mean_income
        )
    )                                                            |> 
    pivot_wider(
        names_from = "year", 
        values_from = "mean_income"
    )


## -----------------------------------------------------------------------------
mean_income |>
    pivot_longer(
        names_to = "year",
        values_to = "mean_income",
        !country
    )


## -----------------------------------------------------------------------------
mean_income |>
    pivot_longer(
        names_to = "year",
        values_to = "mean_income",
        !country
    )                                                            |> 
    group_by(
        country
    )                                                            |> 
    summarise(
        per_country_mean = mean(mean_income, na.rm = TRUE)
    )


## -----------------------------------------------------------------------------
mean_income                                                      |>
    pivot_longer(
        names_to = "year",
        values_to = "mean_income",
        !country
    )                                                            |> 
    group_by(
        country
    )                                                            |> 
    mutate(
        mean_per_country = mean(mean_income, na.rm = TRUE),
        mean_income = ifelse(
            is.na(mean_income),
            mean_per_country,
            mean_income
        )
    )


## -----------------------------------------------------------------------------
mean_income                                                      |>
    pivot_longer(
        names_to = "year",
        values_to = "mean_income",
        !country
    )                                                            |> 
    group_by(
        country
    )                                                            |> 
    # Imagine we are doing something here that removes the group...
    ungroup()                                                    |>
    # The mutate that follows is not grouped so the mean is global...
    mutate(
        mean_per_country = mean(mean_income, na.rm = TRUE),
        mean_income = ifelse(
            is.na(mean_income),
            mean_per_country,
            mean_income
        )
    )


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
strings <- c(
    "Give me an ice cream",
    "Get yourself an ice cream",
    "We are all out of ice creams",
    "I scream, you scream, everybody loves ice cream.",
    "one ice cream,
    two ice creams,
    three ice creams",
    "I want an ice cream. Do you want an ice cream?"
)
str_count(strings)


## -----------------------------------------------------------------------------
str_count(strings, boundary("character"))


## -----------------------------------------------------------------------------
str_count(strings, boundary("word"))


## -----------------------------------------------------------------------------
str_count(strings, boundary("line_break"))
str_count(strings, boundary("sentence"))


## -----------------------------------------------------------------------------
str_count(strings, "ice cream")
str_count(strings, "cream") # gets the screams as well


## -----------------------------------------------------------------------------
str_count(strings, ".")


## -----------------------------------------------------------------------------
str_count(strings, fixed("."))


## -----------------------------------------------------------------------------
str_count(strings, "[[:punct:]]")


## -----------------------------------------------------------------------------
str_count(strings, "ice creams?$")


## -----------------------------------------------------------------------------
str_count(strings, "ice creams?[[:punct:]]?$")


## -----------------------------------------------------------------------------
strings <- c(
    "one",
    "two",
    "one two",
    "one   two",
    "one. two."
)
str_split(strings, " ")


## -----------------------------------------------------------------------------
str_split(strings, "[[:space:]]+")


## -----------------------------------------------------------------------------
str_split(strings, boundary("word"))


## -----------------------------------------------------------------------------
macdonald <- "Old MACDONALD had a farm."
str_to_lower(macdonald)


## -----------------------------------------------------------------------------
str_to_upper(macdonald)


## -----------------------------------------------------------------------------
str_to_sentence(macdonald)


## -----------------------------------------------------------------------------
str_to_title(macdonald)


## -----------------------------------------------------------------------------
strings <- c(
    "Give me an ice cream",
    "Get yourself an ice cream",
    "We are all out of ice creams",
    "I scream, you scream, everybody loves ice cream.",
    "one ice cream,
    two ice creams,
    three ice creams",
    "I want an ice cream. Do you want an ice cream?"
)
str_wrap(strings)


## -----------------------------------------------------------------------------
str_wrap(strings, width = 10)


## -----------------------------------------------------------------------------
str_wrap(strings, width = 10, indent = 2)


## -----------------------------------------------------------------------------
str_pad(strings, width = 50)


## -----------------------------------------------------------------------------
str_pad(strings, width = 50, side = "right")


## -----------------------------------------------------------------------------
str_pad(strings, width = 50, side = "both")


## -----------------------------------------------------------------------------
strings |> str_trunc(25) |> str_pad(width = 25, side = "left")


## -----------------------------------------------------------------------------
str_trim(c(
    " one small coke",
    "two large cokes  ",
    " three medium cokes "
))


## -----------------------------------------------------------------------------
str_trim(c(
    " one small  coke",
    "two   large cokes  ",
    " three  medium cokes "
))


## -----------------------------------------------------------------------------
str_squish(c(
    " one small  coke",
    "two   large cokes  ",
    " three  medium cokes "
))


## -----------------------------------------------------------------------------
str_detect(strings, "me")
str_detect(strings, "I")
str_detect(strings, "cream")


## -----------------------------------------------------------------------------
str_detect(strings, "ice cream.")


## -----------------------------------------------------------------------------
str_detect(strings, fixed("ice cream."))


## -----------------------------------------------------------------------------
str_detect(strings, "ice cream\\.")


## -----------------------------------------------------------------------------
str_detect(strings, fixed("ice cream."), negate = TRUE)


## -----------------------------------------------------------------------------
str_starts(strings, "I")
str_ends(strings, fixed("."))


## -----------------------------------------------------------------------------
str_locate(strings, "ice cream")


## -----------------------------------------------------------------------------
strings[6]


## -----------------------------------------------------------------------------
str_locate(strings[6], "ice cream")


## -----------------------------------------------------------------------------
str_locate_all(strings[6], "ice cream")


## -----------------------------------------------------------------------------
ice_cream_locations <- str_locate_all(strings[6], "ice cream")
ice_cream_locations
invert_match(ice_cream_locations[[1]])


## -----------------------------------------------------------------------------
str_extract(strings, "(s|ice )cream\\w*")


## -----------------------------------------------------------------------------
strings[4]
str_extract(strings[4], "(s|ice )cream\\w*")
str_extract_all(strings[4], "(s|ice )cream\\w*")


## -----------------------------------------------------------------------------
lego_str <- str_replace(strings, "ice cream[s]?", "LEGO")
lego_str
lego_str <- str_replace(lego_str, "an LEGO", "a LEGO")
lego_str


## -----------------------------------------------------------------------------
strings %>%
    str_replace("ice cream[s]?", "LEGO") %>%
    str_replace("an LEGO", "a LEGO")


## -----------------------------------------------------------------------------
strings %>%
    str_replace_all("ice cream[s]?", "LEGO") %>%
    str_replace_all("an LEGO", "a LEGO")


## -----------------------------------------------------------------------------
us_dates <- c(
    valentines = "2/14",
    my_birthday = "2/15",
    # no one knows but let's just go with this
    jesus_birthday = "12/24"
)
# US date format to a more sane format
str_replace(us_dates, "(.*)/(.*)", "\\2/\\1")


## -----------------------------------------------------------------------------
str_c(
    "NA",
    str_dup("-NA", times = 7),
    " BATMAN!"
)


## -----------------------------------------------------------------------------
# -- concatenation -------------------------------------------
c("foo", "bar", "baz")
str_c("foo", "bar", "baz")


## -----------------------------------------------------------------------------
my_string <- "this is my string"
my_location <- str_locate(my_string, "my")
my_location
s <- my_location[,"start"]
e <- my_location[,"end"]
str_sub(my_string, s, e)

my_string_location <- str_locate(my_string, "string")
s <- my_string_location[,"start"]
e <- my_string_location[,"end"]
str_sub(my_string, s, e)

your_string <- my_string
s <- my_location[,"start"]
e <- my_location[,"end"]
str_sub(your_string, s, e) <- "your"
your_string

your_banana <- your_string
your_string_location <- str_locate(your_string, "string")
s <- your_string_location[,"start"]
e <- your_string_location[,"end"]
str_sub(your_banana, s, e) <- "banana"
your_banana


## -----------------------------------------------------------------------------
my_string
your_string
your_banana


## -----------------------------------------------------------------------------
macdonald <- "Old MacDonald"
eieio <- "E-I-E-I-O"
str_glue("{macdonald} had a farm. {eieio}")


## -----------------------------------------------------------------------------
str_glue(
    "{macdonald} had a farm. {eieio}",
    macdonald = "Thomas",
    eieio = "He certainly did not!"
)


## -----------------------------------------------------------------------------
str_glue(
    "{str_dup(\"NA-\", times = 7)}NA BATMAN!"
)


## -----------------------------------------------------------------------------
x <- seq(1:10)
str_glue(
  "Holy {mean(x)} BATMAN!"
  )



## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## ---- echo=FALSE--------------------------------------------------------------
set.seed(5)


## -----------------------------------------------------------------------------
factor(c("A", "C", "B"))


## -----------------------------------------------------------------------------
factor(c("A", "C", "B"), levels = c("A", "B", "C", "D"))
factor(c("A", "C", "B"), levels = c("D", "C", "B", "A"))


## -----------------------------------------------------------------------------
factor(c("A", "C", "B"), levels = c("A", "B"))


## -----------------------------------------------------------------------------
f1 <- as.factor(c("C", "B", "B", "A", "D"))
f1 # levels ordered alphabetically


## -----------------------------------------------------------------------------
f2 <- as_factor(c("C", "B", "B", "A", "D"))
f2 # levels in the order they appear in the input


## -----------------------------------------------------------------------------
f1
f2


## -----------------------------------------------------------------------------
fct_c(f1, f2)


## -----------------------------------------------------------------------------
fct_c(f2, f1)


## -----------------------------------------------------------------------------
f3 <- as_factor(c("X", "Y", "A"))


## -----------------------------------------------------------------------------
fct_c(f1, f3)


## -----------------------------------------------------------------------------
fct_c(f3, f1)


## -----------------------------------------------------------------------------
fct_c(f1, f2, f3)
fct_c(f2, f3, f1)
fct_c(f3, f1, f2)


## -----------------------------------------------------------------------------
fct_unify(fs = list(f1, f2, f3))


## -----------------------------------------------------------------------------
fct_collapse(
    fct_c(f3, f1),
    a = c("A", "X"),
    b = c("B", "Y"),
    c = c("C", "D")
)


## -----------------------------------------------------------------------------
fct_collapse(
    fct_c(f3, f1),
    a = c("A", "X"),
    b = c("B", "Y")
)


## -----------------------------------------------------------------------------
f1
fct_recode(f1, a = "A", b = "B")
fct_recode(f1, X = "A", X = "B", Y = "C", Y = "D")


## -----------------------------------------------------------------------------
f <- factor(sample(1:10, 20, replace = TRUE), levels = 1:10)
table(f)
f |> fct_lump(n = 5, other_level = "X") |> table()
f |> fct_lump(n = 2, other_level = "X") |> table()


## -----------------------------------------------------------------------------
f |> fct_lump(n = -2, other_level = "X") |> table()


## -----------------------------------------------------------------------------
f |> fct_lump(prop = 0.1, other_level = "X") |> table()


## -----------------------------------------------------------------------------
f |> fct_lump(prop = -0.1, other_level = "X") |> table()


## -----------------------------------------------------------------------------
f |> fct_other(keep = 1:5, other_level = "X")
f |> fct_other(drop = 1:5, other_level = "X")


## -----------------------------------------------------------------------------
f1
levels(f1) <- LETTERS[1:10]
f1
fct_drop(f1)


## -----------------------------------------------------------------------------
f1
f3
fct_expand(f1, levels(f3))


## -----------------------------------------------------------------------------
f2
fna <- f2
fna[2] <- NA
fna[3] <- NA
fna
fna <- fct_explicit_na(fna, na_level = "Missing")
fna


## -----------------------------------------------------------------------------
f <- factor(sample(LETTERS[1:5], 10, replace = TRUE))
f


## -----------------------------------------------------------------------------
fct_inorder(f)


## -----------------------------------------------------------------------------
fct_inorder(f, ordered = TRUE)


## -----------------------------------------------------------------------------
factor(f, levels = sort(levels(f)))
factor(f, levels = sort(levels(f)), ordered = TRUE)


## -----------------------------------------------------------------------------
table(f)
fct_infreq(f)
f |> fct_infreq() |> table()
f |> fct_infreq(ordered = TRUE)
f |> fct_infreq() |> fct_rev() |> table()


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## ---- echo=FALSE--------------------------------------------------------------
suppressPackageStartupMessages(library(lubridate, quietly = TRUE))


## -----------------------------------------------------------------------------
ymd("1975 Feb 15")
ymd("19750215")
ymd("1975/2/15")
ymd("1975-02-15")


## -----------------------------------------------------------------------------
dmy("150275")
mdy("February 15th 1975")


## -----------------------------------------------------------------------------
dmy_h("15/2/1975 2pm")
dmy_hm("15/2/1975 14:30")
dmy_hms("15/2/1975 14:30:10")


## -----------------------------------------------------------------------------
x <- dmy_hms("15/2/1975 14:30:10")


## -----------------------------------------------------------------------------
c(day(x), month(x), year(x))
c(hour(x), minute(x), second(x))
c(week(x), # The week in the year
  wday(x), # The day in the week
  yday(x)) # The day in the year


## -----------------------------------------------------------------------------
minute(x) <- 15
wday(x) <- 42
x


## -----------------------------------------------------------------------------
dmy_hm(
  "15/2/1975 14:00", 
  tz = "Europe/Copenhagen"
)


## -----------------------------------------------------------------------------
force_tz(
  # This is a date/time in CET
  dmy_hm("15/2/1975 14:00", tz = "Europe/Copenhagen"),
  # It will be moved to GMT
  tz = "Europe/London"
) 


## -----------------------------------------------------------------------------
with_tz(
  # This is a date/time in CET
  dmy_hm("15/2/1975 14:00", tz = "Europe/Copenhagen"),
  # This moves it to the same time but in GMT
  tz = "Europe/London"
)


## -----------------------------------------------------------------------------
start <- dmy("02 11 1949")
end <- dmy("15 02 1975")
interval(start, end)


## -----------------------------------------------------------------------------
start %--% end


## -----------------------------------------------------------------------------
int <- interval(start, end)
int
int_start(int)
int_end(int)


## -----------------------------------------------------------------------------
end %--% start
int_start(start %--% end)
int_start(end %--% start)


## -----------------------------------------------------------------------------
int_flip(end %--% start)


## -----------------------------------------------------------------------------
int_standardize(start %--% end)
int_standardize(end %--% start)


## -----------------------------------------------------------------------------
x <- now()
int <- interval(x, x + minutes(1)) # from now and one minute forward
int_length(int) # the length is one minute, so 60 seconds

int <- interval(x, x + minutes(20)) # now and 20 minutes forward
int_length(int) / 60 # Dividing by 60 to get the length in minutes


## -----------------------------------------------------------------------------
ymd("1867 05 02") %within% int
ymd("1959 04 23") %within% int
x %within% int # start point is inside the interval
(x + minutes(20)) %within% int # end point is inside the interval


## -----------------------------------------------------------------------------
int_start(int) <- dmy("19 Aug 1950")
int
int_end(int) <- dmy("19 Sep 1950")
int


## -----------------------------------------------------------------------------
int_shift(int, months(1))


## -----------------------------------------------------------------------------
int1 <- interval(dmy("19 Aug 1950"), dmy("19 Sep 1950"))
int2 <- interval(dmy("19 oct 1950"), dmy("25 nov 1951"))
int3 <- interval(dmy("19 oct 1948"), dmy("25 aug 1951"))
int4 <- interval(dmy("19 oct 1981"), dmy("25 aug 2051"))

# int1 ends before int2
int_overlaps(int1, int2)
# int3 starts before int1 but they overlap
int_overlaps(int1, int3)
# no overlap, int4 is far in the future compared to int1
int_overlaps(int1, int4)


## -----------------------------------------------------------------------------
c(
  int_aligns(int1, int2),
  int_aligns(int1, int3),
  int_aligns(int1, int4)
)


## -----------------------------------------------------------------------------
int5 <- interval(int_start(int1), int_end(int1) + years(3))
int6 <- int_shift(int5, -years(3))
int7 <- int_shift(int6, -years(3))

c(
  int_aligns(int1, int5), # share start
  int_aligns(int1, int6), # share end
  int_aligns(int1, int7)  # overlaps but does not share endpoints
)


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)
library(broom)
library(modelr)


## -----------------------------------------------------------------------------
model <- lm(disp ~ hp + wt, data = mtcars)
summary(model)


## -----------------------------------------------------------------------------
tidy(model) # transform to tidy tibble
glance(model) # model summaries
augment(model) # add model info to data


## -----------------------------------------------------------------------------
# Build a model where variable x can help us predict response y
dat <- tibble(
  x = runif(50),
  y = 15 * x^2 * x + 42 + rnorm(5)
)
# Fit a linear model to the data (even though y is quadratic in x)
model <- lm(y ~ x, data = dat)
tidy(model)


## -----------------------------------------------------------------------------
add_predictions(dat, model)
add_residuals(dat, model)


## -----------------------------------------------------------------------------
new_dat <- tibble(x = seq(0, 1, length.out = 5))
add_predictions(new_dat, model)


## -----------------------------------------------------------------------------
seq(0, 1, length.out = 5)
seq_range(dat$x, n = 5) # over the range of observations


## -----------------------------------------------------------------------------
# comparing line to a (better) model y ~ x^2 + x + 1
model2 <- lm(y ~ I(x^2) + x, data = dat)
gather_predictions(new_dat, model, model2)
spread_predictions(new_dat, model, model2)


## -----------------------------------------------------------------------------
gather_predictions(dat, model, model2)
spread_predictions(dat, model, model2)


## -----------------------------------------------------------------------------
gather_residuals(dat, model, model2)
spread_residuals(dat, model, model2)


## -----------------------------------------------------------------------------
bootstrap(dat, n = 3)


## -----------------------------------------------------------------------------
crossv_mc(dat, n = 3)


## -----------------------------------------------------------------------------
crossv_kfold(dat, k = 3)
crossv_loo(dat)


## -----------------------------------------------------------------------------
samples <- bootstrap(dat, 3)
fitted_models <- samples |> 
  mutate(
    # Map over all the bootstrap samples and fit each of them
    fits = strap |> map(\(dat) lm(y ~ x, data = dat))
  )
fitted_models
fitted_models$fits[[1]]


## -----------------------------------------------------------------------------
fitted_models$fits |>
  map(glance) |>
  bind_rows()


## -----------------------------------------------------------------------------
get_x <- function(m) {
  tidy(m) |> filter(term == "x") |>
             select(estimate) %>% as.double()
}
fitted_models$fits |> map_dbl(get_x)


## -----------------------------------------------------------------------------
models <- formulae(~y, linear = ~x, quadratic = ~I(x^2) + x)


## -----------------------------------------------------------------------------
fits <- fit_with(dat, lm, models)
fits |> map(glance) |> bind_rows()


## -----------------------------------------------------------------------------
fits |> map_dbl(rmse, data = dat)


## -----------------------------------------------------------------------------
fits |> map_dbl(mae, data = dat)


## -----------------------------------------------------------------------------
fits |> map_dbl(~ glance(.x)$AIC)


## -----------------------------------------------------------------------------
samples <- dat |> crossv_loo()
training_fits <- samples$train |> map(~lm(y ~ x, data = .))
training_fits |> map2_dbl(samples$test, rmse) |> head(10)


## ---- include=FALSE-----------------------------------------------------------
knitr::opts_knit$set(root.dir = "~/Library/Mobile Documents/com~apple~CloudDocs/Books/R Data Science Quick Reference")
library(tidyverse)


## -----------------------------------------------------------------------------
p <- ggplot()


## -----------------------------------------------------------------------------
summary(p)


## -----------------------------------------------------------------------------
dat <- tibble(
  foo = runif(100),
  bar = 20 * foo + 5 + rnorm(100, sd = 10),
  baz = rep(1:2, each = 50)
)
p <- ggplot(data = dat)


## -----------------------------------------------------------------------------
p <- ggplot(data = dat, aes(x = foo, y = bar, color = baz))


## -----------------------------------------------------------------------------
p <- ggplot(data = dat, aes(x = foo, y = bar, color = baz)) +
       geom_point()


## ----simplest_geom_point, fig.cap="Point geometry plot"-----------------------
print(p)


## ----discrete-colour, fig.cap="Discrete colour aesthetics"--------------------
ggplot(data = dat, aes(x = foo, y = bar, color = factor(baz))) +
  geom_point()


## ----line-plot, fig.cap="A line plot"-----------------------------------------
ggplot(data = arrange(dat, foo), 
       aes(x = foo, y = bar, color = factor(baz))) + 
  geom_line()


## ---- warnings=FALSE----------------------------------------------------------
p <- ggplot(data = dat, 
            aes(x = foo, y = bar, color = baz)) + 
  geom_point() + 
  geom_smooth(formula = y ~ x, method = "loess")


## ----combined-point-smooth, fig.cap="Plot with two geometries"----------------
print(p)


## ---- combined-point-smooth-discrete, fig.cap="Plot with two geometries and a discrete colour"----
ggplot(data = dat, aes(x = foo, y = bar, color = factor(baz))) + 
  geom_point() + 
  geom_smooth(formula = y ~ x, method = "loess")


## ----hist-plot, fig.cap="Histogram plot"--------------------------------------
ggplot(data = dat, aes(x = foo)) + 
  geom_histogram(binwidth = 0.05)


## ---- facet-grid, fig.cap="Faceting the plot"---------------------------------
ggplot(data = dat, aes(x = foo, y = bar)) + 
  geom_point() + facet_grid(~ factor(baz))


## ----two-dim-facet, fig.cap="Facet grid for two variables."-------------------
dat2 <- tibble(
  foo = rep(1:5, each = 20),
  bar = rep(1:2, each = 50),
  x = foo * bar + rnorm(100),
  y = -foo
)
ggplot(data = dat2, aes(x = x, y = y)) + 
  geom_point() + facet_grid(factor(foo) ~ factor(bar))


## ----multi-variables-facets, fig.cap="Facet with four variables"--------------
dat3 <- tibble(
  foo = factor(rep(1:5, each = 20)),
  bar = factor(rep(1:2, each = 50)),
  baz = factor(rep(1:5, times = 20)),
  qux = factor(rep(1:2, times = 50)),
  x = rnorm(100),
  y = rnorm(100)
)
ggplot(data = dat3, aes(x = x, y = y)) + 
  geom_point() + facet_grid(foo + bar ~ baz + qux)


## ----coord-flip, fig.cap="Plot with the x-coordinate flipped."----------------
ggplot(data = dat, aes(x = foo, y = bar)) + 
  geom_point() + coord_flip()


## ----polar, fig.cap="Plot in polar coordinates."------------------------------
ggplot(data = dat, aes(x = foo, y = bar)) + 
  geom_point() + coord_polar()


## ----scaled_y_color, fig.cap="Plot with rescaled y-coordinate and colors."----
ggplot(data = dat, aes(x = bar, y = foo, color = baz)) + 
  geom_point() + 
  scale_y_log10() + 
  scale_color_binned()

