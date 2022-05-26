#' The GetSurvivalProbe function
#' @description Get probes whose methylation state is predictive of patient survival
#' @param EpiMixResults List of objects returned from the EpiMix function
#' @param TCGA_CancerSite String indicating the TCGA cancer code (e.g. 'LUAD')
#' @param clinical.data (If the TCGA_CancerSite is specified, this parameter is optional) Dataframe with survival information. Must contain at least three columns: 'sample.id', 'days_to_death', 'days_to_last_follow_up'.
#' @param raw.pval.threshold numeric value indicting the raw p value threshold for selecting the survival predictive probes. Survival time is compared by log-rank test. Default: 0.05
#' @param p.adjust.method character string indicating the statistical method for adjusting multiple comparisons, can be either of 'holm', 'hochberg', 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', 'none'. Default: 'fdr'
#' @param adjusted.pval.threshold numeric value indicting the adjusted p value threshold for selecting the survival predictive probes. Default: 0.05
#' @param OutputRoot path to save the output. If not null, the return value will be saved as 'Survival)Probes.csv'.
#' @return a dataframe with probes whose methylation state is predictive of patient survival and the p value.
#' @export
#' @examples
#' \donttest{
#' library(survival)
#'
#' data('Sample_EpiMixResults_miRNA')
#'
#' survival.CpGs <- GetSurvivalProbe(EpiMixResults = Sample_EpiMixResults_miRNA,
#'                                      TCGA_CancerSite = 'LUAD')
#'
#' }
#'

GetSurvivalProbe <- function(EpiMixResults, TCGA_CancerSite = NULL, clinical.data = NULL,
    raw.pval.threshold = 0.05, p.adjust.method = "none", adjusted.pval.threshold = 0.05,
    OutputRoot = "") {

    if (!requireNamespace("survival")) {
        message("This function requires the 'survival' package.")
        return(invisible())
    }

    if (is.null(TCGA_CancerSite) & is.null(clinical.data)) {
        stop("Please provide the value for either TCGA_CancerSite or clinical.data")
    }

    FunctionalPairs <- EpiMixResults$FunctionalPairs

    # Group the genes associated with the same CpG together
    ProbeAnnotation <- FunctionalPairs %>%
        dplyr::group_by(.data$Probe) %>%
        dplyr::mutate(Genes = paste(.data$Gene, collapse = ";")) %>%
        dplyr::select(.data$Probe, .data$Genes)

    # Download clinical data
    if (!is.null(TCGA_CancerSite)) {
        cat(paste0("Downloading patient clinical data for ", TCGA_CancerSite, "\n"))
        clinical.data.directory <- get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite,
            saveDir = tempdir(), dataFileTag = "Merge_Clinical.Level_1", )
        clinical.data <- data.table::fread(paste0(clinical.data.directory, TCGA_CancerSite,
            ".merged_only_clinical_clin_format.txt"))
        clinical.data <- as.matrix(clinical.data)
        rownames(clinical.data) <- clinical.data[, 1]
        clinical.data <- clinical.data[, -1]
        clinical.data <- clinical.data[-1, ]
        survival.info <- data.frame(sample.id = toupper(clinical.data["patient.bcr_patient_barcode",
            ]), days_to_death = as.numeric(clinical.data["patient.days_to_death",
            ]), days_to_last_follow_up = as.numeric(clinical.data["patient.days_to_last_followup",
            ]))
    } else {
        survival.info <- clinical.data
    }

    # Construct a dataframe for survival information. Status: 1 = censored, 2 =
    # dead
    survival.info <- survival.info %>%
        dplyr::mutate(time = dplyr::if_else(is.na(.data$days_to_death), .data$days_to_last_follow_up,
            .data$days_to_death), status = dplyr::if_else(is.na(.data$days_to_death),
            1, 2))
    survival.info <- survival.info[!is.na(survival.info$time), ]
    survival.info <- survival.info[survival.info$time != 0, ]

    # Survival analysis
    meth.state <- getMethStates(EpiMixResults, EpiMixResults$MethylationDrivers)
    DMvalues <- EpiMixResults$MethylationStates
    DMvalues <- TCGA_GENERIC_CleanUpSampleNames(DMvalues, 12)

    cat("Finding survival-associated CpGs\n")
    target.probes <- rownames(DMvalues)
    iterations <- length(target.probes)
    Probe <- State <- pval <- hazard.ratio <- low.conf <- high.conf <- character(0)
    for (i in seq(1:iterations)) {
        target.probe <- state <- meth.values <- abnormal <- normal <- mixture.group <- target.survival <- sdf <- P.value <- sde <- HR <- low.cl <- high.cl <- NULL
        target.probe <- target.probes[i]
        state <- meth.state[target.probe]
        DM.value <- DMvalues[target.probe, ]
        if (state == "Hypo") {
            abnormal <- names(DM.value)[which(DM.value < 0)]
            normal <- names(DM.value)[which(DM.value == 0)]
        } else if (state == "Hyper") {
            abnormal <- names(DM.value)[which(DM.value > 0)]
            normal <- names(DM.value)[which(DM.value == 0)]
        } else {
            abnormal <- names(DM.value)[which(DM.value < 0)]
            normal <- names(DM.value)[which(DM.value > 0)]
        }

        if (length(abnormal) < 10 | length(normal) < 10) {
            (next)()
        }
        normal <- data.frame(sample.id = normal, State = 1)
        abnormal <- data.frame(sample.id = abnormal, State = 2)
        mixture.group <- rbind(normal, abnormal)
        target.survival <- merge(survival.info, mixture.group)
        sdf <- survival::survdiff(Surv(time, status) ~ State, data = target.survival)
        P.value <- 1 - stats::pchisq(sdf$chisq, length(sdf$n) - 1)
        sde <- survival::coxph(Surv(time, status) ~ State, data = target.survival)
        HR <- stats::coef(summary(sde))[, 2]
        low.cl <- exp(stats::confint(sde))[, 1]
        high.cl <- exp(stats::confint(sde))[, 2]
        Probe <- c(Probe, target.probe)
        State <- c(State, state)
        pval <- c(pval, P.value)
        hazard.ratio <- c(hazard.ratio, HR)
        low.conf <- c(low.conf, low.cl)
        high.conf <- c(high.conf, high.cl)
    }

    adjusted.pval <- p.adjust(pval, method = p.adjust.method)
    survival.results <- data.frame(Probe = Probe, State = State, HR = hazard.ratio,
        lower.Cl = low.conf, higher.Cl = high.conf, p.value = pval, adjusted.p.value = adjusted.pval)

    # Add annotation to CpGs
    survival.results <- merge(survival.results, ProbeAnnotation)
    survival.results <- survival.results %>%
        dplyr::select(.data$Probe, .data$Genes, .data$State, .data$HR, .data$lower.Cl,
            .data$higher.Cl, .data$p.value, .data$adjusted.p.value) %>%
        dplyr::distinct() %>%
        dplyr::filter(.data$p.value < raw.pval.threshold) %>%
        dplyr::filter(.data$adjusted.p.value < adjusted.pval.threshold) %>%
        dplyr::arrange(.data$adjusted.p.value, .data$Genes, .data$HR, decreasing = TRUE)
    rownames(survival.results) <- NULL
    cat("Found", nrow(survival.results), "survival predictive CpGs\n")

    if (OutputRoot != "" & length(survival.results) > 0) {
        cat(paste0("Saving the result to ", OutputRoot, "/", "Survival_CpGs.csv"))
        utils::write.csv(survival.results, paste0(OutputRoot, "/", "Survival_CpGs.csv"),
            row.names = FALSE)
    }
    return(survival.results)
}

#' EpiMix_PlotSurvival function
#' @description function to plot Kaplan-meier survival curves for patients with different methylation state of a specific probe.
#' @param EpiMixResults List of objects returned from the EpiMix function
#' @param plot.probe Character string with the name of the probe
#' @param TCGA_CancerSite TCGA cancer code (e.g. 'LUAD')
#' @param clinical.df (If the TCGA_CancerSite parameter has been specified, this parameter is optional) Dataframe with survival information. Must contain at least three columns: 'sample.id', 'days_to_death', 'days_to_last_follow_up'.
#' @param font.legend numeric value indicating the font size of the figure legend. Default: 16
#' @param font.x numeric value indicating the font size of the x axis label. Default: 16
#' @param font.y numeric value indicating the font size of the y axis label. Default: 16
#' @param font.tickslab numeric value indicating the font size of the axis tick label. Default: 14
#' @param legend numeric vector indicating the x,y coordinate for positioning the figure legend. c(0,0) indicates bottom left, while c(1,1) indicates top right. Default: c(0.8,0.9). If 'none', legend will be removed.
#' @param show.p.value logic indicating whether to show p value in the plot. P value was calculated by log-rank test.  Default: TRUE.
#'
#' @return Kaplan-meier survival curve showing the survival time for patients with different methylation states of the probe.
#' @export
#' @examples
#' \donttest{
#' library(survival)
#' library(survminer)
#'
#' data(Sample_EpiMixResults_miRNA)
#'
#' EpiMix_PlotSurvival(EpiMixResults = Sample_EpiMixResults_miRNA,
#'                     plot.probe = 'cg00909706',
#'                     TCGA_CancerSite = 'LUAD')
#' }
#'

EpiMix_PlotSurvival <- function(EpiMixResults, plot.probe, TCGA_CancerSite = NULL,
    clinical.df = NULL, font.legend = 16, font.x = 16, font.y = 16, font.tickslab = 14,
    legend = c(0.8, 0.9), show.p.value = TRUE) {

    if (!requireNamespace("survival")) {
        message("This function requires the 'survival' package.")
        return(invisible())
    }

    if (!requireNamespace("survminer")) {
        message("This function requires the 'survival' package.")
        return(invisible())
    }

    Classifications <- EpiMixResults$Classifications
    Classifications <- TCGA_GENERIC_CleanUpSampleNames(Classifications, 12)

    mixture.group <- Classifications[plot.probe, ]
    mixture.group <- data.frame(sample.id = names(mixture.group), State = mixture.group)
    rownames(mixture.group) <- NULL

    # Download clinical data
    if (!is.null(TCGA_CancerSite)) {
        file_path <- paste0("gdac_20160128/gdac.broadinstitute.org_", TCGA_CancerSite,
            ".Merge_Clinical.Level_1.2016012800.0.0/", TCGA_CancerSite, ".merged_only_clinical_clin_format.txt")
        if (!file.exists(file_path)) {
            cat(paste0("Downloading patient clinical data for ", TCGA_CancerSite))
            clinical.data.directory <- get_firehoseData(TCGA_acronym_uppercase = TCGA_CancerSite,
                dataFileTag = "Merge_Clinical.Level_1", saveDir = tempdir())
            clinical.data <- data.table::fread(paste0(clinical.data.directory, TCGA_CancerSite,
                ".merged_only_clinical_clin_format.txt"))
        } else {
            clinical.data <- data.table::fread(file_path)
        }
        clinical.data <- as.matrix(clinical.data)
        rownames(clinical.data) <- clinical.data[, 1]
        clinical.data <- clinical.data[, -1]
        clinical.data <- clinical.data[-1, ]
        survival.info <- data.frame(sample.id = toupper(clinical.data["patient.bcr_patient_barcode",
            ]), days_to_death = as.numeric(clinical.data["patient.days_to_death",
            ]), days_to_last_follow_up = as.numeric(clinical.data["patient.days_to_last_followup",
            ]))
    } else {
        survival.info <- clinical.df
    }

    # Construct a dataframe for survival information. Status: 1 = censored, 2 =
    # dead
    survival.info <- survival.info %>%
        mutate(time = dplyr::if_else(is.na(.data$days_to_death), .data$days_to_last_follow_up,
            .data$days_to_death), status = dplyr::if_else(is.na(.data$days_to_death),
            1, 2))
    survival.info <- survival.info[!is.na(survival.info$time), ]
    survival.info <- survival.info[survival.info$time != 0, ]
    target.survival <- merge(survival.info, mixture.group)
    target.survival <- target.survival[order(target.survival$State), ]
    survival <- survminer::ggsurvplot(survminer::surv_fit(survival::Surv(time, status) ~
        State, data = target.survival), pval = show.p.value, legend.title = "mixture component",
        legend.labs = unique(target.survival$State), xlab = "Days", ylab = "Overall survival probability",
        palette = RColorBrewer::brewer.pal(8, "Set1")[1:length(unique(target.survival$State))],
        font.legend = c(font.legend, "bold"), font.x = c(font.x, "bold"), font.y = c(font.y,
            "bold"), font.tickslab = c(font.tickslab), legend = legend)
    return(survival)
}
