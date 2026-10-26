###############################################################################
#  build_IST_dataset.R
#
#  PURPOSE
#    Rebuild our analysis dataset from the original, publicly released
#    International Stroke Trial (IST) data, one documented step at a time.
#    Anyone can run this script and get exactly the same files we used.
#
#  WHERE THE DATA COMES FROM
#    Study:     International Stroke Trial (IST), 1991 to 1996.
#               19,435 patients with acute stroke, randomised in a 2 x 3
#               factorial design: aspirin (yes/no) x subcutaneous heparin
#               (none / low dose 5,000 IU twice daily / medium dose 12,500 IU
#               twice daily).
#    Released:  Sandercock P, Niewada M, Czlonkowska A. International Stroke
#               Trial database (version 2). University of Edinburgh,
#               Department of Clinical Neurosciences, 2011.
#    DOI:       https://doi.org/10.7488/ds/104
#    Web page:  https://datashare.ed.ac.uk/handle/10283/124
#    Licence:   Open Data Commons Attribution License (ODC-By). Free to use
#               and share as long as the source is cited (see citations at
#               the bottom of this header).
#    Papers:    Data description: Sandercock et al. (2011) "The International
#               Stroke Trial database", Trials 12:101.
#               Main results: International Stroke Trial Collaborative Group
#               (1997), Lancet 349:1569-1581.
#
#  HOW WE DOWNLOADED IT
#    1. Go to https://datashare.ed.ac.uk/handle/10283/124
#    2. Click "Download all files" to get the zip file DS_10283_124.zip
#       (about 51 MB). It contains the data, the variable list (codebook),
#       the trial protocol and the data paper.
#    3. Put the zip in the same folder as this script.
#
#    Inside the zip we use ONE data file: IST_corrected.csv.
#    It is "version 2" of the database. Version 1 had errors in some date
#    variables introduced during anonymisation; the authors fixed them in
#    version 2 (see "IST version 2 additional commentary.pdf" in the zip).
#    The codebook we follow is IST_variables.csv / IST_variables.pdf.
#
#  WHAT THIS SCRIPT PRODUCES (saved in the same folder)
#    IST_data_matrix_full.csv   19,435 rows x 67 columns (all useful variables)
#    IST_data_matrix.csv        19,435 rows x 21 columns (core subset)
#    IST_build_log.txt          record of the input file, checks and outputs
#
#  CODING RULES USED THROUGHOUT
#    * One row = one patient.
#    * Yes/No questions:   Y -> 1,  N -> 0,  anything else -> missing (blank).
#      ("U" = unknown and "C" = can't assess are treated as missing.)
#    * Nothing is imputed or made up. Missing stays missing.
#
#  HOW TO CITE IN THE REPORT
#    Sandercock P, Niewada M, Czlonkowska A (2011). International Stroke Trial
#    database (version 2) [dataset]. University of Edinburgh.
#    https://doi.org/10.7488/ds/104
###############################################################################


## ---------------------------------------------------------------------------
## STEP 0. Set the folder that holds this script and the downloaded zip
## ---------------------------------------------------------------------------
## Change this path to your own folder if needed.
project_dir <- "~/Desktop/Columbia Organized/Coursework/Fall 2026/Advance DA"
if (dir.exists(project_dir)) setwd(project_dir)
cat("Working folder:", getwd(), "\n")


## ---------------------------------------------------------------------------
## STEP 1. Unzip the download and locate the data file
## ---------------------------------------------------------------------------
zip_file  <- "DS_10283_124.zip"        # the file downloaded from DataShare
data_file <- "IST_raw/IST_corrected.csv"

if (!file.exists(data_file)) {
  if (!file.exists(zip_file)) {
    message("Could not find ", zip_file, " in this folder. Please pick it.")
    zip_file <- file.choose()
  }
  unzip(zip_file, exdir = "IST_raw")    # extracts all files into IST_raw/
}
stopifnot(file.exists(data_file))

# Fingerprint of the input file. If someone else downloads the same file,
# they should get this same code, which proves it is the identical data.
input_md5 <- unname(tools::md5sum(data_file))
cat("Input file:", data_file, "\nMD5 fingerprint:", input_md5,
    "(expected c233a72eef0aa30b23569b2017ed17dd)\n")


## ---------------------------------------------------------------------------
## STEP 2. Read the raw data
## ---------------------------------------------------------------------------
## The file was saved with Windows-1252 text encoding (a few hospital
## comments contain accented letters). If you read it with a plain
## read.csv(), R stops silently after 146 rows. So we first convert the text
## to standard UTF-8, then read it.
raw_lines <- readLines(data_file, warn = FALSE, encoding = "bytes")
raw_lines <- iconv(raw_lines, from = "CP1252", to = "UTF-8", sub = "?")
raw <- read.csv(text = raw_lines, stringsAsFactors = FALSE,
                na.strings = c("", "NA"))

cat("Raw data:", nrow(raw), "patients x", ncol(raw), "variables\n")
stopifnot(nrow(raw) == 19435)   # the trial randomised 19,435 patients
stopifnot(ncol(raw) == 112)     # the released file has 112 variables


## ---------------------------------------------------------------------------
## STEP 3. Helper functions for recoding
## ---------------------------------------------------------------------------
# Y -> 1, N -> 0, everything else (U, C, blank) -> NA
yes_no <- function(x) ifelse(x == "Y", 1L, ifelse(x == "N", 0L, NA_integer_))

# Keep a value only if it is one of the allowed codes, otherwise NA
keep_codes <- function(x, allowed) ifelse(x %in% allowed, x, NA)


## ---------------------------------------------------------------------------
## STEP 4. Fix the two coding differences between the pilot and main trial
## ---------------------------------------------------------------------------
## (a) Heparin dose. Codebook for RXHEP: "M/L/N [M is coded as H=high in
##     pilot]". The 245 pilot patients coded "H" received the same 12,500 IU
##     dose as "M" in the main trial, so we recode H to M.
heparin_code <- ifelse(raw$RXHEP == "H", "M", raw$RXHEP)
cat("Heparin codes after recoding H to M:\n"); print(table(heparin_code))

## (b) Pilot phase. Codebook for RATRIAL: "Atrial fibrillation (Y/N); not
##     coded for pilot phase - 984 patients". We flag those 984 patients so
##     analysts can include or exclude them on purpose.
is_pilot <- as.integer(is.na(raw$RATRIAL))
stopifnot(sum(is_pilot) == 984)


## ---------------------------------------------------------------------------
## STEP 5. Build the analysis matrix
##         Each line: new_name = recode(RAW_NAME)   # codebook meaning
## ---------------------------------------------------------------------------
six_month_status <- keep_codes(raw$OCCODE, 1:4)   # OCCODE 0 and 9 = unknown

ist <- data.frame(

  ## --- Identifiers -------------------------------------------------------
  id       = seq_len(nrow(raw)),        # row number (no patient identifiers exist)
  hospital = raw$HOSPNUM,               # HOSPNUM: hospital number
  country  = raw$COUNTRY,               # COUNTRY: abbreviated country code
  pilot    = is_pilot,                  # 1 = pilot phase patient (see Step 4b)

  ## --- Randomised treatment (what the patient was ASSIGNED) ---------------
  aspirin      = yes_no(raw$RXASP),                       # RXASP: trial aspirin allocated (Y/N)
  heparin      = as.integer(heparin_code != "N"),         # RXHEP: any trial heparin (L or M)
  heparin_dose = match(heparin_code, c("N", "L", "M")) - 1L,  # 0 none, 1 low, 2 medium

  ## --- Patient at randomisation ------------------------------------------
  age           = raw$AGE,                                # AGE: age in years
  female        = as.integer(raw$SEX == "F"),             # SEX: M/F
  sbp           = raw$RSBP,                               # RSBP: systolic blood pressure (mmHg)
  af            = yes_no(raw$RATRIAL),                    # RATRIAL: atrial fibrillation
  consciousness = match(raw$RCONSC, c("F", "D", "U")) - 1L,  # RCONSC: 0 fully alert, 1 drowsy, 2 unconscious

  ## --- The stroke and its timing -----------------------------------------
  hours_to_treatment = raw$RDELAY,                        # RDELAY: hours from stroke to randomisation
  on_waking          = yes_no(raw$RSLEEP),                # RSLEEP: symptoms noted on waking
  stroke_type        = raw$STYPE,                         # STYPE: TACS/PACS/POCS/LACS/OTH (Oxford classification)
  ct_before          = yes_no(raw$RCT),                   # RCT: CT scan before randomisation
  infarct_on_ct      = yes_no(raw$RVISINF),               # RVISINF: infarct visible on CT
  def_face           = yes_no(raw$RDEF1),                 # RDEF1: face deficit
  def_arm            = yes_no(raw$RDEF2),                 # RDEF2: arm/hand deficit
  def_leg            = yes_no(raw$RDEF3),                 # RDEF3: leg/foot deficit
  def_dysphasia      = yes_no(raw$RDEF4),                 # RDEF4: dysphasia (speech)
  def_hemianopia     = yes_no(raw$RDEF5),                 # RDEF5: hemianopia (vision)
  def_visuospatial   = yes_no(raw$RDEF6),                 # RDEF6: visuospatial disorder
  def_brainstem      = yes_no(raw$RDEF7),                 # RDEF7: brainstem/cerebellar signs
  def_other          = yes_no(raw$RDEF8),                 # RDEF8: other deficit
  weekday            = raw$DAYLOCAL,                      # DAYLOCAL: estimated local day of week (1-7)
  hour_local         = keep_codes(raw$HOURLOCAL, 0:23),   # HOURLOCAL: local hour (99 = unknown)

  ## --- Medication before entering the trial -------------------------------
  heparin_prior_24h = yes_no(raw$RHEP24),                 # RHEP24: heparin within 24 h before
  aspirin_prior_3d  = yes_no(raw$RASP3),                  # RASP3: aspirin within 3 days before

  ## --- Final diagnosis of the stroke (decided after randomisation) --------
  dx_ischaemic     = yes_no(toupper(raw$DDIAGISC)),       # DDIAGISC: ischaemic stroke (3 entries were lowercase "u")
  dx_haemorrhagic  = yes_no(raw$DDIAGHA),                 # DDIAGHA: haemorrhagic stroke
  dx_indeterminate = yes_no(raw$DDIAGUN),                 # DDIAGUN: indeterminate stroke
  dx_not_stroke    = yes_no(raw$DNOSTRK),                 # DNOSTRK: not a stroke

  ## --- Treatment actually received (adherence and other drugs) -----------
  days_on_trial_drug    = raw$ONDRUG,                     # ONDRUG: days on trial treatment
  compliant_aspirin     = yes_no(raw$CMPLASP),            # CMPLASP: compliant for aspirin
  compliant_heparin     = yes_no(raw$CMPLHEP),            # CMPLHEP: compliant for heparin
  nontrial_antiplatelet = yes_no(raw$DAP),                # DAP: non-trial antiplatelet drug
  nontrial_sc_heparin   = yes_no(raw$DSCH),               # DSCH: non-trial subcutaneous heparin
  nontrial_iv_heparin   = yes_no(raw$DIVH),               # DIVH: non-trial intravenous heparin
  other_anticoagulant   = yes_no(raw$DOAC),               # DOAC: other anticoagulants
  thrombolysis          = yes_no(raw$DTHROMB),            # DTHROMB: thrombolysis
  steroids              = yes_no(raw$DSTER),              # DSTER: steroids
  aspirin_long_term_at_discharge = yes_no(raw$DASPLT),    # DASPLT: discharged on long-term aspirin

  ## --- Outcomes within 14 days (the trialists' own derived indicators) ---
  dead_14d                        = raw$ID14,             # ID14: died within 14 days
  recurrent_ischaemic_14d         = raw$ISC14,            # ISC14: recurrent ischaemic stroke
  brain_bleed_14d                 = raw$H14,              # H14: cerebral bleed / haemorrhagic stroke
  recurrent_unknown_type_14d      = raw$NK14,             # NK14: recurrent stroke of unknown type
  any_stroke_14d                  = raw$STRK14,           # STRK14: any recurrent stroke
  haemorrhagic_transformation_14d = raw$HTI14,            # HTI14: bleeding into the original infarct
  major_bleed_14d                 = raw$TRAN14,           # TRAN14: major non-cerebral bleed
  any_extracranial_bleed_14d      = raw$NCB14,            # NCB14: any non-cerebral bleed
  pe_14d                          = raw$PE14,             # PE14: pulmonary embolism
  dvt_14d                         = raw$DVT14,            # DVT14: deep vein thrombosis
  dead_or_stroke_14d = as.integer(raw$ID14 == 1 | raw$STRK14 == 1),  # composite we built
  discharged_alive   = yes_no(raw$DALIVE),                # DALIVE: discharged alive from hospital

  ## --- Outcomes at 6 months ------------------------------------------------
  outcome_6mo           = six_month_status,               # OCCODE: 1 dead, 2 dependent, 3 not recovered, 4 recovered
  dead_or_dependent_6mo = ifelse(is.na(six_month_status), NA,
                                 as.integer(six_month_status <= 2)),  # the trial's primary outcome
  dead_6mo      = yes_no(raw$FDEAD),                      # FDEAD: dead at 6-month follow-up
  recovered_6mo = yes_no(raw$FRECOVER),                   # FRECOVER: fully recovered (survivors only)
  dependent_6mo = yes_no(raw$FDENNIS),                    # FDENNIS: dependent (survivors only)
  home_6mo = ifelse(raw$FPLACE %in% c("A", "B"), 1L,      # FPLACE: A home, B relative's home
             ifelse(raw$FPLACE %in% c("C", "D", "E"), 0L, NA)),  # C/D/E care home or hospital
  cause_of_death = keep_codes(ifelse(!is.na(raw$FDEADC), raw$FDEADC, raw$DDEADC), 1:8),
                                                          # FDEADC/DDEADC: 1 initial stroke, 2 recurrent ischaemic,
                                                          # 3 recurrent haemorrhagic, 4 pneumonia, 5 coronary,
                                                          # 6 pulmonary embolism, 7 other vascular, 8 non-vascular

  ## --- Survival -------------------------------------------------------------
  days_to_death_or_censor = raw$TD,                       # TD: time of death or censoring (days)
  died                    = raw$DIED,                     # DIED: 1 died, 0 did not

  ## --- Trialists' predicted risks (baseline severity scores) ---------------
  pred_dead_or_dep_6mo = raw$EXPDD,                       # EXPDD: predicted prob. dead/dependent at 6 months
  pred_dead_6mo        = raw$EXPD6,                       # EXPD6: predicted prob. dead at 6 months
  pred_dead_14d        = raw$EXPD14                       # EXPD14: predicted prob. dead at 14 days
)

## Variables we deliberately did NOT carry over:
##   free-text comments (DMAJNCHX, DSIDEX, DNOSTRKX, DDEADX, FDEADX),
##   administrative dates (FU1_RECD, FU2_DONE, FU1_COMP, RDATE),
##   raw discharge-form fields duplicated by the derived 14-day indicators
##   (DRSISC, DRSH, DRSUNK, DPE and their dates),
##   DEAD1 to DEAD8 (combined into cause_of_death), and dose-detail fields
##   (DASP14, DLH14, DMH14, DHH14) summarised by the compliance variables.


## ---------------------------------------------------------------------------
## STEP 6. Check the result against the published trial
## ---------------------------------------------------------------------------
## If any check fails, the script stops, so a successful run means all passed.
checks <- c(
  "19,435 patients"                    = nrow(ist) == 19435,
  "9,720 aspirin / 9,715 no aspirin"   = all(table(ist$aspirin) == c(9715, 9720)),
  "Heparin 9,718 none / 4,861 low / 4,856 medium"
                                       = all(table(ist$heparin_dose) == c(9718, 4861, 4856)),
  "984 pilot patients"                 = sum(ist$pilot) == 984,
  "599 haemorrhagic strokes"           = sum(ist$dx_haemorrhagic, na.rm = TRUE) == 599,
  "420 not a stroke"                   = sum(ist$dx_not_stroke, na.rm = TRUE) == 420,
  "Dead or dependent 6 mo, aspirin 62.2%"
      = round(100 * mean(ist$dead_or_dependent_6mo[ist$aspirin == 1], na.rm = TRUE), 1) == 62.2,
  "Dead or dependent 6 mo, no aspirin 63.5%"
      = round(100 * mean(ist$dead_or_dependent_6mo[ist$aspirin == 0], na.rm = TRUE), 1) == 63.5
)
print(data.frame(check = names(checks), passed = checks), row.names = FALSE)
stopifnot(all(checks))
## The 6-month figures match the published result (Lancet 1997):
## 62.2% dead or dependent with aspirin vs 63.5% without.


## ---------------------------------------------------------------------------
## STEP 7. Save the datasets and a log
## ---------------------------------------------------------------------------
write.csv(ist, "IST_data_matrix_full.csv", row.names = FALSE, na = "")

core_columns <- c("id", "country", "pilot", "aspirin", "heparin", "heparin_dose",
                  "age", "female", "sbp", "af", "consciousness",
                  "hours_to_treatment", "stroke_type", "dx_ischaemic",
                  "dead_14d", "recurrent_ischaemic_14d", "brain_bleed_14d",
                  "major_bleed_14d", "dead_or_dependent_6mo", "dead_6mo",
                  "outcome_6mo")
write.csv(ist[, core_columns], "IST_data_matrix.csv", row.names = FALSE, na = "")

missing_counts <- colSums(is.na(ist))
log_lines <- c(
  paste("Built:", format(Sys.time(), "%Y-%m-%d %H:%M")),
  "Source: International Stroke Trial database v2, https://doi.org/10.7488/ds/104 (ODC-By)",
  paste("Input file:", data_file, " MD5:", input_md5),
  paste("Raw dimensions:", nrow(raw), "x", ncol(raw)),
  paste("Full matrix:", nrow(ist), "x", ncol(ist), "-> IST_data_matrix_full.csv"),
  paste("Core matrix:", nrow(ist), "x", length(core_columns), "-> IST_data_matrix.csv"),
  "All checks against the published trial passed:", paste(" ", names(checks)),
  "Missing values per column (blank cells):",
  paste(" ", names(missing_counts[missing_counts > 0]), missing_counts[missing_counts > 0])
)
writeLines(log_lines, "IST_build_log.txt")
cat("\nDone. Wrote IST_data_matrix_full.csv, IST_data_matrix.csv and IST_build_log.txt\n")
