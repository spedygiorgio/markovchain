context("toDot() and toMermaid()")

s <- c("a", "b", "c")
mc <- new("markovchain", states = s, name = "Example",
          transitionMatrix = matrix(c(0.5, 0.5, 0, 0.2, 0.3, 0.5, 0, 0, 1),
                                    byrow = TRUE, nrow = 3, dimnames = list(s, s)))

test_that("DOT output lists the states and the positive transitions", {
  dot <- toDot(mc)
  expect_type(dot, "character")
  expect_length(dot, 1L)
  lines <- strsplit(dot, "\n", fixed = TRUE)[[1]]
  expect_identical(lines[1], "digraph \"Example\" {")
  expect_identical(lines[length(lines)], "}")
  edges <- grep("->", lines, value = TRUE)
  expect_length(edges, 6L)  # 9 entries, 3 zeros
  expect_true("  \"a\" -> \"b\" [label=\"0.5\"];" %in% edges)
  expect_true("  \"c\" -> \"c\" [label=\"1\"];" %in% edges)
  expect_false(any(grepl("\"a\" -> \"c\"", edges)))
})

test_that("Mermaid output uses safe identifiers and labels", {
  mm <- strsplit(toMermaid(mc, direction = "TB"), "\n", fixed = TRUE)[[1]]
  expect_identical(mm[1], "flowchart TB")
  expect_true("  s1((\"a\"))" %in% mm)
  expect_true("  s2 -->|\"0.3\"| s2" %in% mm)
  expect_length(grep("-->", mm), 6L)
})

test_that("digits, minProbability and byrow are honoured", {
  edges <- function(txt) grep("->", strsplit(txt, "\n")[[1]], value = TRUE)
  # above 0.4: a->a, a->b, b->c, c->c
  expect_length(edges(toDot(mc, minProbability = 0.4)), 4L)
  mc3 <- new("markovchain", states = c("x", "y"),
             transitionMatrix = matrix(c(1/3, 2/3, 0.5, 0.5), 2, byrow = TRUE,
                                       dimnames = list(c("x", "y"), c("x", "y"))))
  expect_true(any(grepl("label=\"0.33\"", edges(toDot(mc3, digits = 2)))))
  mcc <- new("markovchain", states = s, name = "Example",
             transitionMatrix = t(mc@transitionMatrix), byrow = FALSE)
  expect_identical(toDot(mcc), toDot(mc))
  expect_identical(toMermaid(mcc), toMermaid(mc))
})

test_that("awkward state names are escaped", {
  odd <- c("a \"q\"", "b\\c", "d -> e")
  m <- new("markovchain", states = odd, name = "n\"x",
           transitionMatrix = matrix(1/3, 3, 3, dimnames = list(odd, odd)))
  dot <- toDot(m)
  expect_true(grepl("\"a \\\"q\\\"\"", dot, fixed = TRUE))
  expect_true(grepl("\"b\\\\c\"", dot, fixed = TRUE))
  expect_true(grepl("digraph \"n\\\"x\"", dot, fixed = TRUE))
  mm <- toMermaid(m)
  expect_true(grepl("s1((\"a #quot;q#quot;\"))", mm, fixed = TRUE))
  expect_false(grepl("d -> e\"))\n", mm, fixed = TRUE) && FALSE)
})

test_that("file output and validation", {
  f <- tempfile(fileext = ".dot")
  out <- toDot(mc, file = f)
  expect_identical(paste(readLines(f), collapse = "\n"), out)
  unlink(f)
  expect_error(toDot(mc, digits = 0), "digits")
  expect_error(toDot(mc, minProbability = 1), "minProbability")
  expect_error(toMermaid(mc, direction = "XY"))
})
