library(shiny)
library(knitr)
ui <- shinyUI(
  fluidPage(
    uiOutput('markdown')
  )
)
server <- function(input, output) {
  output$markdown <- renderUI({
    HTML(markdown::markdownToHTML(knit('index.Rmd', quiet = TRUE)))
  })
}
shinyApp(ui, server)