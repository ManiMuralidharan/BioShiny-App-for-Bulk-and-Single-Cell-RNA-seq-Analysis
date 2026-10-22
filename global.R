# global.R ---------------------------------------------------------------
# Loaded once when the app starts, before ui and server.
# Light packages are attached here; heavy Bioconductor packages are called
# with pkg::fun() inside the modules, so the app starts fast and only needs
# a package when that analysis is actually run.

library(shiny)
library(shinydashboard)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)

options(shiny.maxRequestSize = 500 * 1024^2)  # 500 MB uploads

# Sourced at the top of app.R. Files in R/ are sourced automatically by
# shiny::runApp() (Shiny >= 1.5), so modules do not need source() calls.
