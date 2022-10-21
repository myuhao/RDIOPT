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
