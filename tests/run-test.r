library(knitr)
knit("test.rmd","test.md") ## 10-15 minutes

library(markdown)
options(markdown.HTML.options=union('toc', getOption("markdown.HTML.options")))
options(markdown.HTML.stylesheet=file.path(getwd(), "style.css"))
markdownToHTML("test.md","test.html")

dir.create("test-output")
R.utils::copyFile("test.md", "test-output")
R.utils::copyFile("test.html", "test-output")
R.utils::copyDirectory("figure", "test-output/figure")
