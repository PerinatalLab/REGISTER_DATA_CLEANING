
# convert PARITY to PARITY_F
# per_child parity to per-pregnancy parity
# SEMFR cleaning scripts

# 2017 April 18, Jonas B.
# stored at ~/Dropbox/GIT/SEMFR_CLEAN_THE_DATA/2_renumber_parity_to_parityF.R

# load packages
library(lubridate)
library(dplyr)

recalculateParity = function(dat,thr_d_low,thr_d_upp) {
        
cat("\n\t RECONSTRUCT PREGNNACY-WISE PARITY VALUES (aka PARITET_F)
\t INPUTS:
\t dat: AFTER eliminating mothers without IDs and mom-kid duplications
\t thr_d_low: min number of days-between-bitrhs for two births to be considered of the same pregnancy (recommended = 31)
\t thr_d_upp: max number of days-between-births for mother's current & subsequent parities to be masked (recommended = 300)
\t OUTPUTS:
\t a vector with new parity (PARITET_F) values")

# for restoring order
dat$sq = seq(nrow(dat)) # to restore the original order in the end
        

###############
############### input check
###############

required_columns = c("AR","BFODDAT","BORDF2","BORDNRF2","PARITET","lpnr_mor","lpnr_BARN","sq")
if (any(!required_columns %in% colnames(dat))) { warning("necessary columns were not found!"); break }
if (!thr_d_low %in% 0:100)  warning("very strange threshold thr_d_low was chosen!")
if (thr_d_upp %in% 0:200)  warning("very strange threshold thr_d_upp was chosen!")
if (thr_d_upp %in% 500:100000)  warning("very inefficient threshold thr_d_upp was chosen!")
if (thr_d_low >= thr_d_upp)  warning("thr_d_low is higher than thr_d_upp!")

###############
############### meta prep
###############


# trim for easier handling and preview
dat = dat[,c("AR","BFODDAT","BORDF2","BORDNRF2","PARITET","lpnr_mor","lpnr_BARN","sq")]


###############
############### fix the birthdates
###############

# preview problems that will be solved:
#table(nchar(dat$BFODDAT))  # should be 8, otherwise - we'll fix here
#dat[which(nchar(dat$BFODDAT)==7),]
#dat[which(nchar(dat$BFODDAT)==6),][1:10,]
#dat[which(nchar(dat$BFODDAT)==4),][1:10,]

dat$year = as.integer(substr(dat$BFODDAT,1,4)) # to be compared to BFODDAT
dat$month = as.character(substr(dat$BFODDAT,5,6)) # must remain character string to preserve "0" in the beginning etc
dat$day = as.character(substr(dat$BFODDAT,7,8))   # --//--


# years
#sum(dat$AR != dat$year)
#table(dat$year,useNA = "a")
#all(unique(dat$year) %in% 1973:2013)  # all OK
if (any(!unique(dat$year) %in% 1973:2013)) warning("nonstandard years encountered - needs manual revision!")

# months
#table(dat$month,useNA = "a")
ok_months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
bad_rows = which(! dat$month %in% ok_months)
dat$month[bad_rows] = NA

# days
#table(dat$day,useNA = "a")
ok_days = as.character(c("01","02","03","04","05","06","07","08","09",10:31))
bad_rows = which(! dat$day %in% ok_days)
dat$day[bad_rows] = NA


# April has only 30 days!
dat$day[which((dat$month=="04")&(dat$day=="31"))] = "30" # none found
# June has only 30 days!
dat$day[which((dat$month=="06")&(dat$day=="31"))] = "30" # two instances
# September has only 30 days!
dat$day[which((dat$month=="09")&(dat$day=="31"))] = "30" # two instances
# November has only 30 days!
dat$day[which((dat$month=="11")&(dat$day=="31"))] = "30" # none found

#sum(is.na(dat$year))
#sum(is.na(dat$month))
#sum(is.na(dat$day))

#### impute missing months and days (years are never missing)

table(mnthNA = is.na(dat$month),IDna = is.na(dat$lpnr_BARN))
table(mnthNA = is.na(dat$month),dayNA = is.na(dat$day))


### impute birthdate in problematic pregnancies

###
### scenarion where birth MONTH is missing:
###

mids = unique(dat$lpnr_mor[which(is.na(dat$month))])
#rez = NULL
for (mid in mids) {
        
        rows_selected = which(dat$lpnr_mor==mid)
        tmp = dat[rows_selected,c("lpnr_BARN","lpnr_mor","PARITET","year","month","day")]
        
        if(nrow(tmp)>1) {  # standard scenario
                
        empty_rix = which( is.na(tmp$month))
        filld_rix = which(!is.na(tmp$month))
        if (any(is.na(tmp$day[filld_rix]))) warning(paste("in mid",mid,"- reference have no days"))
        yrs_prb = unique(tmp$year[empty_rix])
        yrs_okk = unique(tmp$year[filld_rix])
        if (length(yrs_prb)>1) warning(paste("in mid",mid,"- more than one problematic year found"))
        if (length(yrs_okk)==0) warning(paste("in mid",mid,"- no good (reference-)years were found"))
        
        # helpers
        ok_months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
        ok_days = as.character(c("01","02","03","04","05","06","07","08","09",10:31))
        
        # grid for problematic year (without exact birthdate)
        grd = expand.grid(yrs_prb,ok_months,ok_days) # grid
        grd$date = ymd(apply(grd,1,function(x) paste(x,collapse="")),quiet = T)
        grd = grd[which(!is.na(grd$date)),]
        
        # grid for good years (with exact birthdate)
        ref = tmp[filld_rix,c("year","month","day")]
        ref$date = ymd(apply(ref,1,function(x) paste(x,collapse="")),quiet = T)
        
        compare = expand.grid(prob=grd$date,ref=ref$date)
        compare$time_difs = as.numeric(compare$prob - compare$ref)
        solution = compare[which(abs(compare$time_difs)==min(abs(compare$time_difs)))[1],]
        
        filler = as.character(solution$prob)  #  use the closest possible date to actual birth with the date
        yr = substr(filler,1,4)
        mo = substr(filler,6,7)
        dy = substr(filler,9,10)
        dat$month[rows_selected[empty_rix]] = mo
        dat$day[rows_selected[empty_rix]] = dy
        
        #print(tmp)
        #cat("\n")
        #rez = rbind(rez,data.frame(mid,mindif=solution$time_difs,stringsAsFactors = F))
        rm(empty_rix,filld_rix,yrs_prb,yrs_okk,ok_days,ok_months,grd,ref,compare,solution,filler,yr,mo,dy) # mnt_prb
        } else {
                dat$month[rows_selected] = "06"
                dat$day[rows_selected] = "15"
        }
        rm(tmp,rows_selected)
}

#hist(rez$mindif,breaks=100,col="grey")
#sort(rez$mindif)

 
###
### scenarion where birth DAY is missing:
###

mids = unique(dat$lpnr_mor[which(is.na(dat$day))])
#rez = NULL
for (mid in mids) {
        rows_selected = which(dat$lpnr_mor==mid)
        tmp = dat[rows_selected,c("lpnr_BARN","lpnr_mor","PARITET","year","month","day")]
        
        if(nrow(tmp)>1) {  # standard scenario
                
        empty_rix = which( is.na(tmp$day))
        filld_rix = which(!is.na(tmp$day))
        
        if (any(is.na(tmp$day[filld_rix]))) warning(paste("in mid",mid,"- reference have no days"))
        if (any(is.na(tmp$month[filld_rix]))) warning(paste("in mid",mid,"- reference have no months"))
        yrs_prb = unique(tmp$year[empty_rix])
        mnt_prb = unique(tmp$month[empty_rix])
        yrs_okk = unique(tmp$year[filld_rix])
        
        if (length(yrs_prb)>1) warning(paste("in mid",mid,"- more than one problematic year found"))
        if (length(mnt_prb)>1) warning(paste("in mid",mid,"- more than one problematic month found"))
        if (length(yrs_okk)==0) warning(paste("in mid",mid,"- no good (reference-)years were found"))
        
        # helpers
        #ok_months = c("01","02","03","04","05","06","07","08","09","10","11","12") 
        ok_days = as.character(c("01","02","03","04","05","06","07","08","09",10:31))
        
        # grid for problematic year (without exact birthdate)
        grd = expand.grid(yrs_prb,mnt_prb,ok_days) # grid
        grd$date = ymd(apply(grd,1,function(x) paste(x,collapse="")),quiet = T)
        grd = grd[which(!is.na(grd$date)),]
        
        # grid for good years (with exact birthdate)
        ref = tmp[filld_rix,c("year","month","day")]
        ref$date = ymd(apply(ref,1,function(x) paste(x,collapse="")),quiet = T)
        
        compare = expand.grid(prob=grd$date,ref=ref$date)
        compare$time_difs = as.numeric(compare$prob - compare$ref)
        solution = compare[which(abs(compare$time_difs)==min(abs(compare$time_difs)))[1],]
        
        filler = as.character(solution$prob)  #  use the closest possible date to actual birth with the date
        dy = substr(filler,9,10)
        dat$day[rows_selected[empty_rix]] = dy
        
        #tmp = dat[rows_selected,c("lpnr_BARN","lpnr_mor","PARITET","year","month","day")]
        #print(tmp)
        #cat("\n")
        #rez = rbind(rez,data.frame(mid,mindif=solution$time_difs,stringsAsFactors = F))
        
        rm(empty_rix,filld_rix,yrs_prb,mnt_prb,yrs_okk,ok_days,grd,ref,compare,solution,filler,dy)
        } else {
                dat$day[rows_selected] = "15"
        }
        rm(rows_selected,tmp)
}

#hist(rez$mindif,breaks=100,col="grey")
#sort(rez$mindif)


# doublecheck
#sum(is.na(dat$year)) # should be 0
#sum(is.na(dat$month)) # should be 0
#sum(is.na(dat$day)) # if not 0, then execute the line below

# final touch
dat$day[which(is.na(dat$day))] = "15" # assign default. won't cause problems


# reformat as date
dat$bfoddat = paste(dat$year,dat$month,dat$day,sep="")
dat$bfoddat = ymd(dat$bfoddat,quiet = T)

if (any(is.na(dat$bfoddat))) warning("dat$bfoddat contains missing dates (row ~ 215)")

###############
############### split the data
###############

multi_moms = unique(dat$lpnr_mor[which(duplicated(dat$lpnr_mor))])
multimom_rix = which(dat$lpnr_mor %in% multi_moms)
singlmom_rix = which(!dat$lpnr_mor %in% multi_moms)

multi = dat[multimom_rix,] # mothers with at least two children
singl = dat[singlmom_rix,] # mothers with only one child in the mfr



###############
############### once-in-the-register mothers (simple case)
###############

####   "singl" mothers keep their original parity values
singl$PARITET_F = singl$PARITET  # F = forlossning (pregnancy-wise-parity)

# idea for the future:
# table(singl$AR,par234567.. = singl$PARITET>1) # potential reaosn of excluding PARITET>1 & AR> ~1990
# table(singl$AR,par345678.. = singl$PARITET>2) # potential reaosn of excluding PARITET>2 & AR> ~1990
# these values are less trustworthy, as there might have been twin pregs outside mfr follow-up 


###############
############### mothers with multiple children/pregnancies
###############

############ the problem here is that some birthdays of children are too close to each other
####   "multi" mothers first get their parity values re-assigned and then get pruned


#thr_d_low = 31 # min number of days-between-bitrhs for two births to be considered of the same pregnancy
#thr_d_upp = 300 # max number of days-between-births for mother's subsequent parities not to be included (-> unreliable)
multi = group_by(multi, lpnr_mor) %>% arrange(lpnr_mor,bfoddat) %>%
        mutate(bddif = bfoddat - lag(bfoddat)) %>% 
        mutate(dumtst = is.na(bddif)|bddif>thr_d_low) %>% 
        mutate(parity=cumsum(dumtst)) %>% ungroup()  # takes 2 minutes for "multi" moms
multi = as.data.frame(multi)
#head(multi)

# preview distribution of day-difference for between-births in same mom
#multi$bddif = as.numeric(multi$bddif)
#hist(multi$bddif,breaks=100,col="grey")
#abline(v=280,col="red")

# note that declared twin-pregs do not cover some inferred twin-pregs
#ix = which(multi$bddif<3)
#table(multi$BORDF2[ix])
#table(multi$BORDF2[ix-1])

# example of problematic parity assignment (not clear whether it was the same pregnancy or not)
#ix = which( (multi$bddif>thr_d_low)&(multi$bddif<thr_d_upp) )
#bad_mids = unique(multi$lpnr_mor[ix])
#bad_rix = which(multi$lpnr_mor %in% bad_mids)
#multi[bad_rix,][1:10,]

# mask some parities as unreliable
ix = which( (multi$bddif>thr_d_low)&(multi$bddif<thr_d_upp) )
ix = ix-1 # one parity above (currently it will be named as different parity value)
prob = multi[ix,] # problematic rows of the multi dataset
parities = 1:max(multi$parity) # all possible parities
mids_prts_bad = NULL
for (i in 1:nrow(prob)) {
        mid = prob$lpnr_mor[i]
        prt = prob$parity[i]
        #parities = 1:max(multi$parity[which(multi$lpnr_mor==mid)]) # all parities for this mom
        prts_sel = parities[which(parities >= prt)] # selected: all parities above (including) this parity
        mid_prt = paste(mid,prts_sel,sep="_")
        mids_prts_bad = c(mids_prts_bad,mid_prt)
        rm(mid,prt,prts_sel,mid_prt)
}

# assign mid-parity identification codes
multi$mids_prts = paste(multi$lpnr_mor,multi$parity,sep="_")

# mask unreliable parity assignments
bad_rix = which(multi$mids_prts %in% mids_prts_bad)
multi$parity[bad_rix] = NA

###
### correct inferred parities for previous pregnancies not followed-up in mfr
###

# estimate the minimum declared parity for each mother
df = group_by(multi,lpnr_mor) %>% summarise(minP=min(PARITET)) %>% ungroup()
df = as.data.frame(df)

# attach estimates to the original
multi = merge(multi,df,by="lpnr_mor",all.x=T)

# correct new parities based on the minimum declared parities
#head(multi)
multi$PARITET_F = multi$parity + multi$minP - 1  # important to substract "1"

# preview of some problematic (but resolved) cases
#table(multi$PARITET,multi$PARITET_F)
#mid = multi$lpnr_mor[which((multi$PARITET==1)&(multi$PARITET_F==5))]
#multi[which(multi$lpnr_mor %in% mid),]

# report on rows with ambiguous (and thus - masked) new parities
cat("\n\n\t RESULTS:")
cat("\n\t",sum(is.na(multi$PARITET_F)),"new parities will be assigned NA values (due to short interpregnancy interval)")
cat("\n\t but otherwise these rows will not be removed")

###############
############### combine output and finish
###############

tmp = rbind( multi[,c("sq","PARITET_F")], singl[,c("sq","PARITET_F")])
tmp = tmp[order(tmp$sq),]
#head(tmp); nrow(tmp)

if (nrow(tmp) != nrow(dat)) warning("input and output lengths do not match!")

return(tmp$PARITET_F)
} # end of function


# dplyr
#?lag
#?lead
# cumsum
# row_number











####### sandbox:

############### obsolete version for determining per-pregnancy parity
# library(digest)
# fun = function(matID,BD,thr_d) {
#         # data formating
#         tt = data.frame(mid=matID,bdt=BD,sq=seq(length(matID)),stringsAsFactors = F)
#         # 1) initial cycle
#         tt = tt %>% arrange(mid, bdt) %>%
#                 group_by(mid) %>% mutate(par=row_number(),dup=duplicated(par))
#         tt = group_by(tt, mid) %>%
#                 arrange(bdt) %>%
#                 mutate(dif = bdt - lag(bdt)) %>% ungroup()
#         tt = as.data.frame(arrange(tt,mid,bdt))
#         tt$dif = as.numeric(tt$dif) # since these are time differences
#         ix = which((tt$dup==FALSE)&(tt$dif<thr_d))
#         tt$par[ix] = tt$par[ix] - 1
#         # 2) looping cycles
#         hsh = digest(tt$par)  # summarize all parities
#         count = 0
#         repeat {
#                 count = count + 1 # counting and reporting cycles (should be 4)
#                 print(count)
#                 tt = tt[,-grep("^dup$",colnames(tt))]
#                 tt = tt %>% arrange(mid, bdt) %>%
#                         group_by(mid) %>% mutate(dup=duplicated(par))
#                 tt = as.data.frame(arrange(tt,mid,bdt))
#                 # update
#                 ix = which((tt$dup==FALSE)&(tt$dif<thr_d))
#                 tt$par[ix] = tt$par[ix] - 1
#                 hsh_new = digest(tt$par)
#                 if (hsh == hsh_new) { break } else {hsh = hsh_new}
#         }
#         tt = tt[,-grep("^dup$|^dif$",colnames(tt))]
#         # 3) parity simplificator
#         tt$cid = paste(tt$mid,tt$par,sep="_")
#         key = data.frame(cid = unique(tt$cid),stringsAsFactors = F)
#         lst = strsplit(key$cid,"_")
#         dd = do.call(rbind.data.frame, lst)
#         colnames(dd) = c("mid","par")
#         key = cbind(key,dd)        
#         key = key %>% arrange(mid, par) %>%
#                 group_by(mid) %>% mutate(par_new=row_number())
#         key = as.data.frame(key)
#         key = key[,c("cid","par_new")]
#         # merge
#         tt = merge(tt,key,by="cid",all=T)
#         tt = tt[order(tt$sq),]
#         return(tt$par_new)
# }
# execute:
# parities = fun(matID = multi$lpnr_mor,BD = multi$bfoddat,thr_d = 90)  # takes about 5 minutes
# toy dataset:
# tt = data.frame(mid=c(rep(1,10),rep(2,10)),
# bdt=c(1,10,10,10,20,30,30,40,50,50, 1,10,10,20,20,30,1,1,1,1),
# stringsAsFactors = F)

