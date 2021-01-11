f1 <- factor(letters[1:3])
f2 <- ordered(letters[1:3], levels = letters[1:4])
f2 <- factor(letters[1:3])

test_that("simple test",
  expect_equal(f1, f2)
)

