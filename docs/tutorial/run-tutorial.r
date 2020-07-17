library(knitr)
knit("tutorial.rmd","tutorial.md") 

library(markdown)
options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
options(markdown.HTML.stylesheet=file.path(getwd(), "style.css"))
markdownToHTML("tutorial.md","tutorial.html")

