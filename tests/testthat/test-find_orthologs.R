test_that("Test Human to fly conversions", {
  res = find_orthologs(genes = c(1232, 7316))
  res %>%
    filter(symbol == "CCR3") %>%
    pull(to_symbol) %>%
    expect_equal("No Ortholog Data")

  res %>%
    filter(symbol == "UBC" & confidence == "high") %>%
    pull(to_symbol) %>%
    expect_equal("Ubi-p63E")

  expect_warning({
    find_orthologs(genes = c(1232, 7315, 7316), from = c("4896", "7227"))
  })

  expect_error(
    find_orthologs(genes = c(1232), from = "1")
  )
})

test_that("Test Special cases", {
  # Test NA

  # Test No orthologs

  # Test ...
  res = find_orthologs(genes = c(1232, 7316))
  res %>%
    filter(symbol == "CCR3") %>%
    pull(to_symbol) %>%
    expect_equal("No Ortholog Data")

  res %>%
    filter(symbol == "UBC" & confidence == "high") %>%
    pull(to_symbol) %>%
    expect_equal("Ubi-p63E")

  expect_warning({
    find_orthologs(genes = c(1232, 7315, 7316), from = c("4896", "7227"))
  })

  expect_error(
    find_orthologs(genes = c(1232), from = "1")
  )
})



test_that("Test use in tibble", {
  df = tibble(
    genes = c("UBC", "APOE", "what?")
  ) %>%
    mutate(
      out = map(genes, find_orthologs)
    ) %>%
    unnest(out)
  expect_equal(ncol(df), 28L)

  # All case should return data.frame of the same shape
  df = tibble(
    genes = c(
      "APOE", # No ortho
      "UBC", # Yes ortho
      NA # NA
      )
  ) %>%
    mutate(out = map(genes, find_orthologs))
  expect_equal(colnames(df$out[[1]]), colnames(df$out[[2]]))
  expect_equal(colnames(df$out[[1]]), colnames(df$out[[3]]))

})
