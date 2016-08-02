
source("~/Dropbox/GIT/SEMFR_CLEAN_THE_DATA/1_cleaning_modules.R")

# initial design of the matrix
if(exists("year_matrix")==TRUE) rm(year_matrix)
generate_year_counts(dat,show=T)

# maternal and fetal ID problems
dat = fun_momID(dat); generate_year_counts(dat,show=F)
dat = fun_kidID(dat); generate_year_counts(dat,show=F)
dat = fun_momkidID(dat); generate_year_counts(dat,show=F)  # should be 4 073 790  remaining 

# nonvalid pregnancies
dat = fun_multipregs(dat); generate_year_counts(dat,show=F)
dat = fun_deadborn(dat); generate_year_counts(dat,show=F)   # remaining 3 958 944
dat = fun_spont1990(dat); generate_year_counts(dat,show=F) # also removes pregs with AR<1990
dat = fun_CSections(dat); generate_year_counts(dat,show=F)

# maternal and fetal problems
dat = fun_matPrecond(dat); generate_year_counts(dat,show=F) # remaining 3 652 319
dat = fun_ICDcodes(dat)  ; generate_year_counts(dat,show=F)
dat = fun_mom_nordic(dat); generate_year_counts(dat,show=F)

# phenotype problems
dat = fun_GAdating(dat) ; generate_year_counts(dat,show=F)
dat = fun_GAmiss(dat)   ; generate_year_counts(dat,show=F)
dat = fun_MHmiss(dat)   ; generate_year_counts(dat,show=F)

fun_visualize_exclusions_by_year(year_matrix)