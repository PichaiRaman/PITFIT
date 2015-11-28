library(shiny)

shinyApp(
  ui = fluidPage(
    textInput("text", "Text", "")
  ),
  server = function(input, output, session) {
    observe({
      query <- parseQueryString(session$clientData$url_search)
      if (!is.null(query[['text']])) {
        updateTextInput(session, "text", value = query[['text']])
      }
    })
  }
)

