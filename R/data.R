#' Data from the Box Lunch Study
#'
#' This dataset was originally made public in the R package visTree (Venkatasubramaniam and Wolfson, 2018). The covariates of interest are
#' laboratory- and questionnaire-based measures of eating habits and attitudes towards food. The response is the average number of calories consumed
#' in a 24 hour period. This data comes from the baseline measurements of the Box Lunch Study (French et al., 2014)
#'
#' @docType data
#' @format
#' \describe{
#'   \item{trt}{Treatment}
#'   \item{sex}{Sex}
#'   \item{bmi0}{BMI}
#'   \item{skcal}{Snacking kilo calories.}
#'   \item{srvgfv0}{Serving size of fruits and vegetables.}
#'   \item{srvgssb}{Serving size of beverages}
#'   \item{kcal24h0}{Response variable; kcals consumed in 24 hours.}
#'   \item{disinhibition}{}
#'   \item{wanting}{}
#'      \item{hunger}{}
#'         \item{liking}{}
#'            \item{resteating}{}
#'   \item{edeq01}{}
#'   \item{edeq02}{}
#'   \item{edeq13}{}
#'   \item{edeq14}{}
#'   \item{edeq15}{}
#'   \item{edeq22}{}
#'   \item{edeq23}{}
#'   \item{edeq25}{}
#'   \item{edeq26}{}
#'   \item{cdrsbody0}{Body image}
#'   \item{weighfreq0}{Weighing frequency}
#'   \item{freqff}{Fast food frequency}
#'   \item{age}{Age}
#'   \item{rrvfood}{Relative reinforcement of food}
#'}
#' @keywords datasets
#' @format A data frame with 226 rows and 26 variables
#' @references Ashwini Venkatasubramaniam and Julian Wolfson (2018). visTree: Visualization of Subgroups for Decision Trees. R package version 0.8.1. https://CRAN.R-project.org/package=visTree
#' @references French SA, Mitchell NR, Wolfson J, Harnack LJ, Jeffery RW, Gerlach AF, Blundell JE, Pentel PR. Portion size effects on weight gain in a free living setting. Obesity. 2014;22(6):1400–5.
#' @examples
#' data(blsdata)
"blsdata"
