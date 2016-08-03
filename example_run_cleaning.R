
source("~/Dropbox/GIT/SEMFR_CLEAN_THE_DATA/1_cleaning_modules.R")

# initial design of the matrix
year_matrix = NULL
generate_year_counts(dat, stage="initial", show=T)

# maternal and fetal ID problems
dat = fun_momID(dat)
dat = fun_kidID(dat)
dat = fun_momkidID(dat)  # should be 4 073 790  remaining 

# nonvalid pregnancies
dat = fun_multipregs(dat)
dat = fun_deadborn(dat)   # remaining 3 958 944
dat = fun_spont1990(dat) # also removes pregs with AR<1990
dat = fun_CSections(dat)

# maternal and fetal problems
dat = fun_matPrecond(dat) # remaining 3 652 319
dat = fun_ICDcodes(dat)  
dat = fun_mom_nordic(dat)

# phenotype problems
dat = fun_GAdating(dat)
dat = fun_GAmiss(dat)
dat = fun_MHmiss(dat)

fun_visualize_exclusions_by_year(year_matrix)
