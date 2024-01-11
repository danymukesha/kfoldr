#' @title kfold
#' @description
#'   Exhaustive search over param_dict calculating binary metrics.
#'
#' @name kfold
#'
#' @param model An object storing bootlist attributes in .model (e.g., modelPLS.model.x_scores_).
#' @param X Predictor variables, where n_samples is the number of samples and n_features is the number of predictors.
#' @param Y Response variables, where n_samples is the number of samples.
#' @param param_dict List of attributes to calculate and return bootstrap confidence intervals.
#' @param folds A positive integer, the number of folds used in the computation. Default is 10.
#' @param bootnum A positive integer, the number of bootstrap samples used in the computation for the plot. Default is 100.
#'
#' @details
#'   The `kfold` class implements an exhaustive search over the specified parameter combinations,
#'   conducting k-fold cross-validation for each combination. It calculates binary metrics,
#'   generates interactive Bokeh plots, and allows assessment of model performance through bootstrap resampling.
#'
#' @section Methods:
#' \describe{
#'   \item{Run}{Runs all necessary methods prior to plot.}
#'   \item{Plot}{Creates a R2/Q2 plot.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Example usage
#'   library(kfoldr)
#'   model <- YourModelClass
#'   X <- YourX
#'   Y <- YourY
#'   param_dict <- YourParamDict
#'   folds <- 10
#'   bootnum <- 100
#'   kfold_instance <- kfold$new(model, X, Y, param_dict, folds, bootnum)
#'   kfold_instance$run()
#'   kfold_instance$plot(metric = "r2q2", grid_line = TRUE)
#' }


library(ggplot2)
library(data.table)
library(ggdist)

# Define the kfold class
kfold <- R6Class("kfold",
                 inherit = list(BaseCrossVal),
                 public = list(
                   initialize = function(model, X, Y, param_dict, folds=10, bootnum=100) {
                     super$initialize(model=model, X=X, Y=Y, param_dict=param_dict, folds=folds, bootnum=bootnum)
                     self$crossval_idx <- createFolds(self$Y, k=folds, list=TRUE, returnTrain = FALSE)
                   },

                   calc_ypred = function() {
                     self$ypred_full <- list()
                     self$ypred_cv <- list()
                     for (params in self$param_list) {
                       params_i <- params
                       model_i <- self$model(params_i)
                       # Full
                       model_i$train(self$X, self$Y)
                       ypred_full_i <- model_i$test(self$X)
                       self$ypred_full <- append(self$ypred_full, list(ypred_full_i))
                       # CV (for each fold)
                       ypred_cv_i <- self$calc_cv_ypred(model_i, self$X, self$Y)
                       self$ypred_cv <- append(self$ypred_cv, list(ypred_cv_i))
                     }
                   },

                   calc_stats = function() {
                     stats_list <- list()
                     for (i in 1:length(self$param_list)) {
                       stats_full_i <- binary_metrics(self$Y, self$ypred_full[[i]])
                       stats_cv_i <- binary_metrics(self$Y, self$ypred_cv[[i]])
                       stats_full_i <- setNames(stats_full_i, paste0(names(stats_full_i), "full"))
                       stats_cv_i <- setNames(stats_cv_i, paste0(names(stats_cv_i), "cv"))
                       stats_cv_i[["R²"]] <- stats_full_i[["R²full"]]
                       stats_cv_i[["Q²"]] <- stats_cv_i[["R²cv"]]
                       stats_combined <- c(stats_full_i, stats_cv_i)
                       stats_list <- append(stats_list, list(stats_combined))
                     }
                     self$table <- self$format_table(stats_list)
                     self$table <- self$table[order(names(self$table))]
                   },

                   run = function() {
                     self$calc_ypred()
                     self$calc_stats()
                     if (self$bootnum > 1) {
                       self$calc_ypred_boot()
                       self$calc_stats_boot()
                     }
                   },

                   calc_ypred_boot = function() {
                     self$ytrue_boot <- list()
                     self$ypred_full_boot <- list()
                     self$ypred_cv_boot <- list()
                     for (i in 1:self$bootnum) {
                       bootidx_i <- sample(1:length(self$Y), length(self$Y), replace=TRUE)
                       newX <- self$X[bootidx_i, ]
                       newY <- self$Y[bootidx_i]
                       ypred_full_nboot_i <- list()
                       ypred_cv_nboot_i <- list()
                       for (params in self$param_list) {
                         model_i <- self$model(params)
                         # Full
                         model_i$train(newX, newY)
                         ypred_full_i <- model_i$test(newX)
                         ypred_full_nboot_i <- append(ypred_full_nboot_i, list(ypred_full_i))
                         # cv
                         ypred_cv_i <- self$calc_cv_ypred(model_i, newX, newY)
                         ypred_cv_nboot_i <- append(ypred_cv_nboot_i, list(ypred_cv_i))
                       }
                       self$ytrue_boot <- append(self$ytrue_boot, list(newY))
                       self$ypred_full_boot <- append(self$ypred_full_boot, list(ypred_full_nboot_i))
                       self$ypred_cv_boot <- append(self$ypred_cv_boot, list(ypred_cv_nboot_i))
                     }
                   },

                   calc_stats_boot = function() {
                     self$full_boot_metrics <- list()
                     self$cv_boot_metrics <- list()
                     for (i in 1:length(self$param_list)) {
                       stats_full_i <- list()
                       stats_cv_i <- list()
                       for (j in 1:self$bootnum) {
                         stats_full <- binary_metrics(self$ytrue_boot[[j]], self$ypred_full_boot[[j]][[i]])
                         stats_full_i <- append(stats_full_i, list(stats_full))
                         stats_cv <- binary_metrics(self$ytrue_boot[[j]], self$ypred_cv_boot[[j]][[i]])
                         stats_cv_i <- append(stats_cv_i, list(stats_cv))
                       }
                       self$full_boot_metrics <- append(self$full_boot_metrics, list(stats_full_i))
                       self$cv_boot_metrics <- append(self$cv_boot_metrics, list(stats_cv_i))
                     }
                   },

                   calc_cv_ypred = function(model_i, X, Y) {
                     ypred_cv_i <- rep(NA, length(Y))
                     for (fold in 1:self$folds) {
                       train_idx <- self$crossval_idx[[fold]]
                       test_idx <- setdiff(1:length(Y), train_idx)
                       X_train <- X[train_idx, ]
                       Y_train <- Y[train_idx]
                       X_test <- X[test_idx, ]
                       model_i$train(X_train, Y_train)
                       ypred_cv_i[test_idx] <- model_i$test(X_test)
                     }
                     return(ypred_cv_i)
                   },

                   format_table = function(stats_list) {
                     table <- data.table(stats_list)
                     param_list_string <- sapply(self$param_list, function(x) toString(x))
                     setnames(table, param_list_string)
                     return(table)
                   },

                   plot = function(metric="r2q2", grid_line=TRUE) {
                     # Choose metric to plot
                     metric_title <- c("ACCURACY", "AUC", "F1-SCORE", "PRECISION", "R²", "SENSITIVITY", "SPECIFICITY")
                     metric_list <- c("acc", "auc", "f1score", "prec", "r2q2", "sens", "spec")
                     metric_idx <- which(metric_list == metric)

                     # Get full, cv, and diff
                     full <- self$table[, 2 * metric_idx + 2]
                     cv <- self$table[, 2 * metric_idx + 1]
                     diff <- abs(full - cv)
                     full_text <- names(self$table)[2 * metric_idx + 2]
                     cv_text <- names(self$table)[2 * metric_idx + 1]
                     diff_text <- paste0("DIFFERENCE (", full_text, " - ", cv_text, ")")

                     # Round full, cv, and diff for hovertool
                     full_hover <- round(full, 2)
                     cv_hover <- round(cv, 2)
                     diff_hover <- round(diff, 2)

                     # Get key, values (as string) from param_dict (key -> title, values -> x axis values)
                     key <- names(self$param_dict)
                     values <- unlist(self$param_dict)
                     values_string <- as.character(values)

                     # Create data frame
                     df <- data.frame(full=full, cv=cv, diff=diff, full_hover=full_hover, cv_hover=cv_hover, diff_hover=diff_hover, values_string=values_string)

                     # Create ggplot
                     p1 <- ggplot(df, aes(x=cv, y=diff)) +
                       geom_line(size=2, color="black", alpha=0.25) +
                       geom_point(aes(x=cv, y=diff, label=values_string), size=7, color="green", alpha=0.7) +
                       geom_text(aes(x=cv, y=diff, label=values_string), size=3, vjust=1.5, hjust=0) +
                       labs(title=paste0(diff_text, " vs ", cv_text), x=cv_text, y=diff_text) +
                       theme_minimal() +
                       theme(text = element_text(size=8))

                     p2_title <- paste0(full_text, " & ", cv_text, " vs no. of components")
                     p2 <- ggplot(df, aes(x=values_string)) +
                       geom_line(aes(y=full), size=2, color="red") +
                       geom_point(aes(y=full), size=8, color="red", fill="white") +
                       geom_line(aes(y=cv), size=2, color="blue") +
                       geom_point(aes(y=cv), size=8, color="blue", fill="white") +
                       labs(title=p2_title, x="components", y="Value") +
                       theme_minimal() +
                       theme(text = element_text(size=8)) +
                       theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))

                     # Combine plots
                     plot_grid(p1, p2, align="hv", rel_widths=c(1, 1))
                   }
                 )
)

# # Create an instance of kfold
# kfold_instance <- kfold$new(model=YourModelClass, X=YourX, Y=YourY, param_dict=YourParamDict, folds=10, bootnum=100)
#
# # Run the methods
# kfold_instance$run()
#
# # Plot the results
# kfold_instance$plot(metric="r2q2", grid_line=TRUE)
