test_that("Convert a single gene", {
  expect_equal(
    convert_to_entrezID("MOBP"),
    c("MOBP" = "4336")
  )

})


test_that("Detect entrezID", {

  expect_equal(
    is_entrezID(c("1234", "3334")),
    TRUE
  )

  expect_equal(
    suppressWarnings(is_entrezID(c("1234", "APOE"))),
    FALSE
  )

  expect_warning(
    is_entrezID("MAPT", silent = FALSE)
  )


})

test_that("Convert a vector of genes", {
  expect_equal(
    convert_to_entrezID(c("MOBP", "APOE")),
    c("MOBP" = "4336", APOE = "348")
  )

  expect_equal(
    convert_to_entrezID(c("MOBP", "APOE", "not a gene")),
    c("MOBP" = "4336", APOE = "348", "not a gene" = NA_character_)
  )
})

test_that("Convert a vector of genes with special cases", {
  expect_error(
    convert_to_entrezID(c("MOBP", NA, "APOE", NA))
  )

  expect_equal(
    convert_to_entrezID(NA),
    c("NA"=NA_character_)
  )
})
