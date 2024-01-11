#' @title BaseCrossVal
#' @description Base class for crossval: kfold.
#'
#' @name BaseCrossVal
#'
#' This an abstract base class providing a common interface for cross-validation.
#'
#' @export
#'
#' @usage
#' \dontrun{
#'   # This is an abstract class, do not instantiate directly.
#' }
#'

library(R6)
BaseCrossVal <- R6Class("BaseCrossVal",
                        # abstract = TRUE,
                        public = list(
                          initialize = function(model, X, Y, param_dict, folds = 10, bootnum = 100) {
                            self$model <- model
                            self$X <- X
                            self$Y <- Y
                            self$param_dict <- param_dict
                            self$param_list <- as.list(ParameterGrid(param_dict))
                            self$folds <- folds
                            self$bootnum <- bootnum
                            self$num_param <- length(param_dict)
                          },

                          calc_ypred = function() {
                            # To be implemented in the derived class
                          },

                          calc_stats = function() {
                            # To be implemented in the derived class
                          },

                          run = function() {
                            # To be implemented in the derived class
                          },

                          plot = function() {
                            # To be implemented in the derived class
                          }
                        )
)
