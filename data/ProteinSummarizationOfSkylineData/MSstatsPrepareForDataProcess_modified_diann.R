library(survival)

#set working directory to current folder of R script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#' Check if a protein can be summarized with TMP
#' @param input data.table
#' @param remove50missing if TRUE, proteins with more than 50% missing values
#' in all runs will not be summarized
#' @return data.table
#' @keywords internal 
.isSummarizable = function(input, remove50missing) {
  n_obs_run = RUN = NULL
  
  if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0)) {
    msg = paste("Can't summarize for protein", unique(input$PROTEIN),
                "because all measurements are missing or censored.")
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg)
    return(NULL)
  }
  
  if (all(is.na(input$n_obs) | input$n_obs == 0)) {
    msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                "because all measurements are missing or censored.")
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg)
    return(NULL)
  } 
  
  if (all(input$n_obs == 1 | is.na(input$n_obs))) {
    msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                "because features have only one measurement across MS runs.")
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg)
    return(NULL)
  }
  
  if (all(is.na(input$newABUNDANCE) | input$newABUNDANCE == 0) | nrow(input) == 0) {
    msg = paste("After removing features which has only 1 measurement,",
                "Can't summarize for protein", unique(input$PROTEIN), 
                "because all measurements are missing or censored.")
    getOption("MSstatsMsg")("INFO", msg)
    getOption("MSstatsLog")("INFO", msg)
    return(NULL)
  }
  
  missing_runs = setdiff(unique(input$RUN), 
                         unique(input[n_obs_run == 0 | is.na(n_obs_run), RUN]))
  if (length(missing_runs) > 0 & length(intersect(missing_runs, as.character(unique(input$RUN))))) { 
    input = input[n_obs_run > 0 & !is.na(n_obs_run), ]
  }
  
  if (remove50missing) {
    if (all(input$prop_features <= 0.5 | is.na(input$prop_features))) {
      msg = paste("Can't summarize for protein", unique(input$PROTEIN), 
                  "because all runs have more than 50% missing values and",
                  "are removed with the option, remove50missing=TRUE.")
      getOption("MSstatsMsg")("INFO", msg)
      getOption("MSstatsLog")("INFO", msg)
      return(NULL)
    }
  }
  input
}

#' Fit Tukey median polish
#' @param input data.table with data for a single protein
#' @param is_labeled logical, if TRUE, data is coming from an SRM experiment
#' @inheritParams MSstatsSummarize
#' @return data.table
#' @keywords internal
.runTukey = function(input, is_labeled, censored_symbol, remove50missing) {
  Protein = RUN = newABUNDANCE = NULL
  
  if (nlevels(input$FEATURE) > 1) {
    tmp_result = .fitTukey(input)
  } else { 
    if (is_labeled) {
      tmp_result = .adjustLRuns(input, TRUE)
    } else {
      tmp_result = input[input$LABEL == "L", 
                         list(RUN, LogIntensities = newABUNDANCE)]
    }
  }
  tmp_result[, Protein := unique(input$PROTEIN)]
  tmp_result
}

#' Fit tukey median polish for a data matrix
#' @inheritParams .runTukey
#' @return data.table
#' @keywords internal
.fitTukey = function(input) {
  LABEL = RUN = newABUNDANCE = NULL
  
  features = as.character(unique(input$FEATURE))
  wide = data.table::dcast(LABEL + RUN ~ FEATURE, data = input,
                           value.var = "newABUNDANCE", keep = TRUE)
  tmp_fitted = MSstats:::median_polish_summary(as.matrix(wide[, features, with = FALSE]))
  wide[, newABUNDANCE := tmp_fitted]
  tmp_result = wide[, list(LABEL, RUN, newABUNDANCE)]
  
  if (data.table::uniqueN(input$LABEL) == 2) {
    tmp_result = .adjustLRuns(tmp_result)
  }
  tmp_result[, list(RUN, LogIntensities = newABUNDANCE)]
}



#' @importFrom data.table uniqueN
#' @importFrom survival survreg Surv
#' @keywords internal
.fitSurvival = function(input) {
  FEATURE = RUN = NULL
  
  missingness_filter = is.finite(input$newABUNDANCE)
  n_total = nrow(input[missingness_filter, ])
  n_features = data.table::uniqueN(input[missingness_filter, FEATURE])
  n_runs = data.table::uniqueN(input[missingness_filter, RUN])
  is_labeled = data.table::uniqueN(input$LABEL) > 1
  countdf = n_total  < n_features + n_runs - 1
  
  # TODO: set.seed here?
  set.seed(100)
  if (is_labeled) {
    if (length(unique(input$FEATURE)) == 1) {
      # with single feature, not converge, wrong intercept
      # need to check
      fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ RUN + ref,
                    data = input, dist = "gaussian")
    } else {
      if (countdf) {
        fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ RUN + ref,
                      data = input, dist = "gaussian")
      } else {
        fit = survreg(Surv(newABUNDANCE, cen, type='left') ~ FEATURE + RUN + ref,
                      data = input, dist = "gaussian")
      }
    }
  } else {
    if (n_features == 1L) {
      fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ RUN,
                    data = input, dist = "gaussian")    
    } else {
      if (countdf) {
        fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ RUN,
                      data = input, dist = "gaussian")
      } else {
        fit = survreg(Surv(newABUNDANCE, cen, type = "left") ~ FEATURE + RUN,
                      data = input, dist = "gaussian")
      }
    }  
  }
  fit
}


#' Get predicted values from a survival model
#' @param input data.table
#' @return numeric vector of predictions
#' @importFrom stats predict
#' @keywords internal
.addSurvivalPredictions = function(input) {
  LABEL = NULL
  
  survival_fit = .fitSurvival(input[LABEL == "L", ])
  predict(survival_fit, newdata = input)
}


.preProcessIntensities2 <-  function(input, log_base) {
  INTENSITY = ABUNDANCE = NULL
  
  if (any(!is.na(input$INTENSITY) & input$INTENSITY < 1, na.rm = TRUE)) {
    # n_smaller_than_1 = sum(!is.na(input$INTENSITY) & input$INTENSITY < 1, 
    #                        na.rm = TRUE)
    # input[, INTENSITY := ifelse(!is.na(INTENSITY) & INTENSITY < 1, 
    #                             1, INTENSITY)]
    # msg = paste("** There are", n_smaller_than_1, 
    #             "intensities which are zero or less than 1.",
    #             "These intensities are replaced with 1",
    #             collapse = " ")
    
    n_equal_zero = sum(!is.na(input$INTENSITY) & input$INTENSITY == 0, 
                           na.rm = TRUE)
    input[, INTENSITY := ifelse(!is.na(INTENSITY) & INTENSITY == 0, 
                                NA, INTENSITY)]
    msg = paste("** There are", n_equal_zero, 
                "intensities which are zero",
                collapse = " ")
    
    getOption("MSstatsLog")("INFO", msg)
    getOption("MSstatsMsg")("INFO", msg)
  } 
  input[, ABUNDANCE := log(INTENSITY, log_base)]
  getOption("MSstatsLog")("INFO",
                          paste("Logarithm transformation with base",
                                log_base,
                                "is done",
                                collapse = " "))
}


MSstatsPrepareForDataProcess2 <- function (input, log_base, fix_missing) {
  input = MSstats:::.checkDataValidity(input, fix_missing = fix_missing)
  input = MSstats:::.updateColumnsForProcessing(input)
  .preProcessIntensities2(input, log_base)
  input = MSstats:::.makeFactorColumns(input)
  input
}




#######################################
# raw = DDARawData
# 
# n <- 100
# idx <- sample(1:nrow(raw), n)
# 
# raw[idx, 10] <- runif(n, 0.000, 0.998)

# msstats
library(MSstats)


method = "TMP"
cens = "NA"
impute = TRUE

DIANN_results <- read.csv("Skyline_DIANN_GPF_msstats.csv", na.strings = "#N/A") 


raw <- SkylinetoMSstatsFormat(DIANN_results)

MSstatsConvert::MSstatsLogsSettings(FALSE)
input = MSstatsPrepareForDataProcess2(raw, 2, NULL)
input = MSstatsNormalize(input, FALSE)
input = MSstatsMergeFractions(input)
input = MSstatsHandleMissing(input, "TMP", FALSE, "NA", 0.999)
input = MSstatsSelectFeatures(input, "all")
processed = getProcessed(input)
input = MSstatsPrepareForSummarization(input, method, impute, cens, FALSE)
input_split = split(input, input$PROTEIN)
summarized = MSstatsSummarize(input_split, method, impute, cens, FALSE, TRUE)
output = MSstatsSummarizationOutput(input, summarized, processed,
                                    method, impute, cens)

saveRDS(output$ProteinLevelData, "DIANN_GPF_msstats_noNorm_noImpu.rda")
