% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_6800.R
\name{read_6800_equations}
\alias{read_6800_equations}
\title{Read gas-exchange equations directly from 6800 xlsx files.}
\usage{
read_6800_equations(filename)
}
\arguments{
\item{filename}{an xlsx file containing 6800 gas-exchange data}
}
\value{
A list of gas-equations.

In principle, this can be made to work for the 6400 as well, but since
that instrument uses a variation of the older xls format, it is hard to
get working in practice.
}
\description{
It is recommended to use \code{\link[=create_equations]{create_equations()}} for recalculating gas-exchange
data, but under some conditions, it may be useful to apply exactly the same
equations as used in the xlsx data file.
}
\details{
Currently, this only works for xlsx files stored by the 6800. this function
extracts xlsx formulas from the file and stores them in a list for use by the
\code{\link[=recalculate]{recalculate()}} function. Note there is no guarantee that the extracted
equations work on any other data files. Since newer versions of the 6800
firmware allows defining custom equations, it is not guaranteed that all
equations can be extracted successfully.
}
\examples{
example <- system.file("extdata", "lowo2.xlsx", package = "gasanalyzer")

# get equations stored in the xlsx file
Eqs <- read_6800_equations(example)

#Inpect how stomatal conductance is calculated:
Eqs$GasEx.gsw
}
\seealso{
\code{\link[=create_equations]{create_equations()}}
}
