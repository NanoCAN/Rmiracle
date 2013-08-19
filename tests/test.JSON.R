context("JSON import")

test_that("RPPA slides can be loaded from MIRACLE", {
	require(Rmiracle)
	slide <- rppa.load(slideIndex="19", baseUrl="http://10.149.64.8:8080/MIRACLE/spotExport/")
	expect_true(nrow(slide) > 1)
	expect_equal(nrow(slide), 7392)
})