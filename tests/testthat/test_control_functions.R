context("Test control functions")
test_that("controlling functions produces error/warning messages properly.", {
  data(cervical)

  expect_warning(discreteControl(method = NULL))

  # Give an error for unavailable classifiers.
  expect_error(classify(data = data.trainS4, method = "unkown",
                        control = trainControl(method = "repeatedcv", classProbs = TRUE)))

  # 'method' can not be NULL
  expect_error(classify(data = data.trainS4, method = NULL,
                        control = trainControl(method = "repeatedcv", classProbs = TRUE)))

  # Reference cannot be a "numeric".
  expect_error(classify(data = data.trainS4, method = "rpart", ref = 2,
                        control = trainControl(method = "repeatedcv", classProbs = TRUE)))

  # Reference cannot be a "logical".
  expect_error(classify(data = data.trainS4, method = "rpart", ref = FALSE,
                        control = trainControl(method = "repeatedcv", classProbs = TRUE)))

  # Class of "data" should be "DESeqDataSet
  expect_error(classify(data = data.train, method = "rpart",
                        control = trainControl(method = "repeatedcv", classProbs = TRUE)))

  # warning:
  # expect_warning(classify(data = data.trainS4, method = "rpart", normalize = "TMM", ref = "T",
  #                         control = trainControl(method = "repeatedcv", number = 2, repeats = 2, classProbs = TRUE)))
})



# expect({
#   sum(rnorm(10) >= 0) == 5
# }, "\n\nCondition is not satisfied.")
