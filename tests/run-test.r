library(knitr)
knit("test.rmd","test.md") ## 15 minutes


library(markdown)
options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
options(markdown.HTML.stylesheet=file.path(getwd(), "style.css"))
options(markdown.HTML.header=file.path(getwd(), "collapsible.html"))
markdownToHTML("test.md","test.html")
