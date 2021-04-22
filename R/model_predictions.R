#' Get predictions from a linear model or generalized linear model
#'
#' If no data are provided, modelPredictions returns a numeric vector predicted
#' values for the sample, functioning as a simple wrapper for fitted.values().
#' If a dataframe with new values for Xs are provided, modelPredictions adds
#' predicted values and SEs for these new data to the dataframe using predict()
#' from car package.
#'
#' If Data=NULL, returns a numeric vector of predicted values for sample.
#' If Data are provided, adds four new columns at the front of the dataframe
#' These variables are named Predicted (prediced value), CILo (lower bound of
#' - 1 SE from Predicted), CIHi (upper bound of + 1 SE), and SE (Standard error
#' of predicted value). NOTE: For GLM, +-1 SE are calculated on the link scale
#' and then converted to the response scale (which will be asymetric) if
#' type = response. If Label is not NULL, than Label is appended to end of these
#'  four variable names.
#'
#' @param model Name of the model that you'll use to generate predictions
#' @param data a dataframe containing cases for predictions. Must include all
#' regressors from model. Default is NULL with predictions returned for the current sample.
#' @param label A string label to append to variable names for predicted values,
#' CIs and SE. Default is NULL with no append
#' @param type 'response' or 'link'. Used only for glm objects. see predict()
#' @return Returns a dataframe of model predictions, standard errors, and
#' confidence intervals.
#' @export
model_predictions <- function (model, data = NULL, label = NULL, type = "response")
{
  if (is.null(data) & class(model)[1] == "lm") {
    return(fitted.values(model))
  }
  else {
    if (is.null(label)) {
      PredictName = "Predicted"
      CILoName = "CILo"
      CIHiName = "CIHi"
      SEName = "SE"
    }
    else {
      PredictName = paste0("Predicted", label)
      CILoName = paste0("CILo", label)
      CIHiName = paste0("CIHi", label)
      SEName = paste0("SE", label)
    }
    Predictions = matrix(data = NA, nrow = nrow(data), ncol = 4,
                         dimnames = list(1:nrow(data), c(PredictName, CILoName,
                                                         CIHiName, SEName)))
    if (class(model)[1] == "lm") {
      CILevel = 1 - 2 * pt(c(1), df = model$df.residual,
                           lower.tail = FALSE)
      Predictions[, 1:3] = predict(model, newdata = data,
                                   interval = "confidence", level = CILevel)
      Predictions[, 4] = Predictions[, 1] - Predictions[,
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if (class(model)[1] == "glm") {
      tmpPred = predict(model, newdata = data, type = "link",
                        se.fit = TRUE)
      upr <- tmpPred$fit + tmpPred$se.fit
      lwr <- tmpPred$fit - tmpPred$se.fit
      fit <- tmpPred$fit
      if (type == "response") {
        fit <- model$family$linkinv(fit)
        upr <- model$family$linkinv(upr)
        lwr <- model$family$linkinv(lwr)
      }
      Predictions[, 1] = fit
      Predictions[, 2] = lwr
      Predictions[, 3] = upr
      Predictions[, 4] = Predictions[, 1] - Predictions[,
                                                        2]
      Predictions = as.data.frame(Predictions)
    }
    if ((class(model)[1] == "lmerMod") || (class(model)[1] ==
                                           "glmerMod")) {
      Predictions[, c(1, 4)] = predictSE(model, data, se.fit = TRUE,
                                         type = type, level = 0, print.matrix = TRUE)
      Predictions[, 2] = Predictions[, 1] - Predictions[,
                                                        4]
      Predictions[, 3] = Predictions[, 1] + Predictions[,
                                                        4]
    }
    if (any(names(data) == PredictName) || any(names(data) ==
                                               CILoName) || any(names(data) == CIHiName) || any(names(data) ==
                                                                                                SEName)) {
      warning("Variable names (Predicted, CILo, CIHi, SE with label PostFix) used in data.  These variables removed before merging in predicted values")
      data[, c(PredictName, CILoName, CIHiName, SEName)] = list(NULL)
    }
    data = data.frame(Predictions, data)
    return(data)
  }
}
