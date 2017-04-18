

##############   CLEAN THE SEMFR DATA   ###################
# modules that do some part of the data cleaning procedure
# by Jonas B. 2016 July
# last update: 2016 Oct24

cat("NEW UPDATES!
    \t + more stringent parental origin check (through other pregs of the same mother)    
    \t + more universal parental origin script (replaces \"nordic\" and \"very_nordic\") 
    \t + previous CS 
    \t + current CS 
    \t + ICD code choice ooption
    \t + PROM
    \t + matHghQC flexible thresholds, value-transfer from other pregs, discordance filters
    \t + maternal Hgh-Wgh bivar filter with BMI-acording-to-age")

cat("\n UNFIXED ERRORS:
    \t fun_visualize_exclusions_by_year(year_matrix) duoda klaidingus skaichius - RESOLVED!")

cat("\n TO BE DONE:
    \t - previous CS: check MDIAG for icd10 O757 (Vaginal birth after previous cesarean)
    \t - exclude iatrogenics based on 658D and O755 codes
    \t - exclude Caesareans based on P034 (fetus and newborn affected by caesarean birth)
    \t - option to choose which maternal conditions to be used in filtering
    \t - twin cleaning based on the twin register")

cat("\n CLEANING FUNCTIONS:
    \t fun_momID - mothers must have an ID
    \t fun_kidID - children must have an ID (will exclude deadborns)
    \t fun_momkidID - mom-child pair must have a unique ID combination
    \t fun_multipregs - remove twins,tripplets are excluded
    \t fun_deadborn - remove stillbirths
    \t fun_matPrecond - remove mothers with diabetes, kidney disease, SLE etc etc 
    \t fun_spont1990 - use only spontaneous deliveries
    \t fun_previousCS - remove previous Caesarean section
    \t fun_currentCS - remove current Caesarean section (if elective)
    \t fun_noPROM - remove deliveries with premature rupture of membranes (ICD indicators)
    \t fun_ICDcodes - remove pregnancies with pregnancy-related medical problems
    \t fun_GAmiss - exclude rows with mising GestAge data
    \t fun_GAdating - use only specified method of GA dating
    \t fun_mHghQC - clean maternal height, exclude rows with missing values
    \t fun_mWghQC - recover maternal weight (does not exclude rows)
    \t fun_mBmiQC - bivariate plot of maternal height and weight, excludes rows
    \t fun_origin - parent(s) must be from selected countries")


library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(hexbin)
library(tidyr)


# this function documents how the exclusions affected each year of the register
generate_year_counts = function(dat, stage, show) {
        year_matrix_current = group_by(dat, AR) %>% summarize(rows=n())
        if(show) print(year_matrix_current)
        
        year_matrix_current$stage = stage
        year_matrix <<- bind_rows(year_matrix, year_matrix_current)
}

# duplicated, missing or nonsensical maternal IDs
fun_momID = function(dat) {
        n_row_before  = nrow(dat)
        cat("MOTHER IDs:
             - before the cleaning (first character of personal number):")
        print(table(substr(dat$lpnr_mor,1,1),useNA="a"))
        dat = dat[which(dat$lpnr_mor>0),] # needed for twin-removal (remaining 4 073 802)
        n_row_after = nrow(dat)
        cat("
            - after the cleaning (first character of personal number):")
        print(table(substr(dat$lpnr_mor,1,1),useNA="a"))
        cat("
            - before: ",n_row_before,", after: ",n_row_after,", removed: ",n_row_before-n_row_after,sep="")
        
        generate_year_counts(dat, "momID", F)
        dat
}

fun_kidID = function(dat) {
        cat("CHILD IDs: first character of personal number:")
        print(table(substr(dat$lpnr_BARN,1,1),useNA="a"))
        cat("nothing will be removed, as a child with NA personal number is a deadborn!")
        dat
}

fun_momkidID = function(dat) {
        # extract subsets of the dataset to be pruned and merged later
        tmp0 = dat[which( (!is.na(dat$lpnr_BARN)) & (!is.na(dat$lpnr_mor))),] # both mom and child have PersNumb
        tmp1 = dat[which( (is.na(dat$lpnr_BARN)) | (is.na(dat$lpnr_mor))),] # at lest one does not have PersNumb
        
        #  mother and child pair was entered twice:
        morbar_unqcode = paste(tmp0$lpnr_BARN,tmp0$lpnr_mor,sep="_") # unique code for mother-child pair
        dpl_codes = unique(morbar_unqcode[which(duplicated(morbar_unqcode))]) # codes that are duplicated
        dpl_rows = which(morbar_unqcode %in% dpl_codes) # rows that contain duplicated pairs
        unq_rows = which(! morbar_unqcode %in% dpl_codes) # rows that contain duplicated pairs
        
        cat("MOM-CHILD:
- observed number of mother-fetus pairs (rows):", nrow(dat),"
- number of pairs (rows) where child has no ID (deadborn):", sum(is.na(dat$lpnr_BARN)),"
- number of pairs (rows) where mother has no ID:", sum(is.na(dat$lpnr_mor)),"
- number of mother-fetus pairs, where both have non-missing personal numbers:", nrow(tmp0),"
- number of unique mom-child pairs that are entered more than once:", length(dpl_codes))
        
        # extract subsets of the dataset to be pruned and merged later
        tmp0_dpl = tmp0[dpl_rows,] # duplicated mother-child pairs
        tmp0_unq = tmp0[unq_rows,] # unique mother-child pairs
        
        # which columns should be used when determining which individuals/pairs should be kept
        cols_to_use = which(colnames(tmp0_dpl) %in% c("AR","MLANGD","GRDBS","MFODLAND","PARITET","MALDER","lpnr_BARN","lpnr_mor"))
        
        cat("
    running pruning scipt which leaves only one row per mom-child pair: the row with most info")
        
        recovered = NULL
        for (b_id in unique(tmp0_dpl$lpnr_BARN)) {
                sub = tmp0_dpl[which(tmp0_dpl$lpnr_BARN==b_id),] 
                sub = sub[which((!is.na(sub$MLANGD))&(!is.na(sub$GRDBS))),] # strictly must be present
                ix = apply(sub[,cols_to_use],1,function(x) sum((is.na(x))|(x=="")))
                ix = which(ix == min(ix))
                if (length(ix)>1) ix = sample(ix,1)
                recovered = rbind(recovered,sub[ix,])
        }
        tmp = rbind(tmp0_unq,recovered,tmp1)
        tmp = tmp[sample(nrow(tmp),replace=F),]
        cat("
    number of remaining rows = ",nrow(tmp)) # remaining 4 073 790
        
        generate_year_counts(tmp, "momkidID", F)
        tmp
}

fun_multipregs = function(dat) {
        # TWINS/multiplets should be removed first, only then remove covar-missing individuals
        # I will (ambiguously) use the same-YEAR-of delivery as an indicator for MULTIPREG pregnancy
        
        cat("MULTIPREGS:
    twins,triplets etc will be determined using two methods:
    1) based on motherID-PregYear duplicates
    2) based on BORDF2 and BORDNRF2 variables        

values for variable Year (first character in value fied):
")
        print(table(substr(dat$AR,1,1),useNA="a"))
        
        
        momYear_code = paste(dat$AR,dat$lpnr_mor,sep="_") # same order as in "dat" !
        multipreg_indic = sort(unique(momYear_code[duplicated(momYear_code)]))
        dat$multiplet = 0
        dat$multiplet[which(momYear_code %in% multipreg_indic)]=1
        rm(momYear_code, multipreg_indic)
        
        cat("
concordance between inferred and declared multiplet status:
")
        print(table(declared=dat$BORDF2,inferred=dat$multiplet,useNA="a"))
        
        #table(dat$BORDF2,dat$BORDNRF2,useNA="a")
        #table(dat$BORDNRF2,useNA="a")
        #table(dat$BORDF2,useNA="a")
        
        good_rix = which((dat$multiplet==0)&(dat$BORDF2==1)&(is.na(dat$BORDNRF2))) # non-twins (singletons)
        bad_rix = which((dat$multiplet==1)|(dat$BORDF2==2)|(!is.na(dat$BORDNRF2))) # warning with NAs! ***
        length(bad_rix) # report
        length(good_rix) # report
        mean(dat$GRDBS[ bad_rix],na.rm=T)  # report
        mean(dat$GRDBS[good_rix],na.rm=T)  # report
        
        cat("
number of twin/triplet+ mom-child pairs (rows) removed: ",length(bad_rix),"
number of singleton children (rows) remaining: ",length(good_rix)) # should be 3 971 491 remaining
        
        # remove multipregs
        dat = dat[good_rix,]
        rm(good_rix,bad_rix)
        
        generate_year_counts(dat, "multipregs", F)
        dat
}

fun_deadborn = function(dat) {
        cat("DEADBORN:
    remove stillbirths based on variable DODFOD and childID missingness.

concordance of childID missingness with DODFOD variable values:
")
        print(table(IDisNA=is.na(dat$lpnr_BARN),DODFOD=dat$DODFOD,useNA = "a"))
        cat("
    (DODFOD: 1 = dies before delivery, 2 = dies during delivery)")
        
        #table(dat$DODFOD) # 1= dies before delivery, 2= during delivery
        
        good_rix = which((is.na(dat$DODFOD))&(!is.na(dat$lpnr_BARN)))
        dat = dat[good_rix,]  
        
        cat("
number of rows remaining = ", nrow(dat)) 
        
        generate_year_counts(dat, "deadborn", F)
        dat 
}

fun_matPrecond = function(dat) {
                cat("MATERNAL PRECONDITIONS
    exclude pregnancies where mother has an ongoing medical condition,
    which is found in one of the MFR's special columns
    ")
                
                preconditions = c("NJURSJUK","DIABETES","EPILEPSI","ULCOLIT","SLE","HYPERTON") # left unused: "URINVINF", "ASTMA"
                
                # estimate prevalence of each condition, as well as effect size on GestAge
                tbl = NULL
                for (colname in preconditions) {
                        precond = as.numeric(!is.na(dat[,colname]))
                        n_ill = sum(precond==1)
                        coefs = coef(summary(lm(dat$GRDBS ~ precond)))
                        df = data.frame(condition=colname,N=n_ill,beta=coefs[2,1],pval=coefs[2,4])
                        tbl = rbind(tbl,df)
                        rm(df,precond,n_ill,coefs)
                }

#  silenced section below is now transfered into a separate module
#                # add  PREVIOUS C-SECTION
#                table(dat$TSECTIO,useNA="a")
#                precond = rep(0,nrow(dat))
#                precond[which(dat$TSECTIO==1)] = 1  # note, some previous CS might still remain (missing data)
#                coefs = coef(summary(lm(dat$GRDBS ~ precond))); coefs
#                df = data.frame(condition="TSECTIO",N=sum(precond==1),beta=coefs[2,1],pval=coefs[2,4])
#                tbl = rbind(tbl,df); rm(df)
                
                cat("
    pregnancies with the following maternal medical conditions will be excluded from the MFR:
    ")
                print(tbl)
                cat("(beta and pval are estimated as an effect on gestational age)
    ")
                
                
                
                # remove all pregnancies with preconditions from the data
                sub = dat[,preconditions]
                
                #  silenced line below is now transfered into a separate module
                # sub$prevCS = NA; sub$prevCS[which(precond==1)] = 1
                
                bad_rows = apply(sub,1,function(x) sum(!is.na(x))); table(bad_rows); sum(bad_rows>0)
                
                cat("
    number of pregnancies with a specific number of maternal medical conditions:
    ")
                print(table(bad_rows))
                
                dat = dat[which(bad_rows==0),]; dim(dat) 
                
                cat("
number of rows selected for exclusion:", sum(bad_rows>0),"
number of rows remaining:",nrow(dat))
                
                generate_year_counts(dat, "matPrecond", F)
                dat
}

fun_spont1990 = function(dat) {
        
        cat("SPONTANEOUS:
            based on variables FLSPONT and FLINDUKT
            (this script will eliminate all pregnancies before year 1990!)")
        
        #df = group_by(dat,AR) %>% summarise(n=n(),s = sum(!is.na(FLINDUKT))) %>% ungroup()
        #t(df)
        #df = group_by(dat,AR) %>% summarise(n=n(),s = sum(!is.na(FLSPONT))) %>% ungroup()
        #t(df)
        cat("
            variable concordance for the years 1973-1989 :
            ")
        sub = dat[which(dat$AR<1990),]
        print(table(spont = sub$FLSPONT,induct = sub$FLINDUKT,useNA="a"))
        
        cat("
            variable concordance for the years 1990-2012 :
            ")
        sub = dat[which(dat$AR>=1990),]
        print(table(spont = sub$FLSPONT,induct = sub$FLINDUKT,useNA="a"))
        
        good_rix = which((dat$FLSPONT==1)&(is.na(dat$FLINDUKT)))
        
        
        cat("after cleaning:
            - pregnancy-year distribution:
            ")
        print(table(dat[good_rix,"AR"]))
        
        cat(" - number of rows removed:",nrow(dat)-length(good_rix),"
            - number of rows remaining:",length(good_rix))
        
        cat("
            NOTE that separate files could be created:
            - when AR>1989 and FLINDUKT is missing and FLSPONT is missing;
            - when AR<1990 and FLINDUKT is missing and FLSPONT is missing;
            ")
        
        dat = dat[good_rix,]
        generate_year_counts(dat, "spont1990", F)
        dat
}

fun_previousCS = function(dat) {
        
        # TO BE DONE:  ICD10 "O757" Vaginal birth after previous cesarean - check MDIAG
        
        cat("CAESAREAN SECTION: \n \t exclusion for previous CS \n\n")
        cat("\t 1) pregnancies with a previous CS (TSECTIO) summary: \n")        
        print(table(dat$TSECTIO,useNA="a"))
        cat("\n pregnancies with a previous CS (TSECTIO) split by year: \n")        
        print(table(dat$TSECTIO,dat$AR,useNA="a"))
        
        cat("\n exclusions will be done if TSECTIO=1 (yes) \n")        
        
        bad_rows = which(dat$TSECTIO==1)
        
        cat("\n in total ",length(bad_rows)," rows will be removed \n")        
        dat = dat[-bad_rows,]
        cat("\n in total ",nrow(dat)," are remaining \n")        
        
        generate_year_counts(dat, "previousCSection", F)
        dat

}

fun_currentCS = function(dat) {
        cat("CAESAREAN SECTION: \n \t exclusion for current CS \n\n")
        cat("\t 1) pregnancies with a current CS (SECFORE) summary: \n")        
        print(table(dat$SECFORE,useNA="a"))
        cat("\n pregnancies with a current CS (SECFORE) split by year: \n")        
        print(table(dat$SECFORE,dat$AR,useNA="a"))
        
        cat("\n \t 2) pregnancies with a current CS (ELEKAKUT) summary: \n")        
        print(table(dat$ELEKAKUT,useNA="a"))
        cat("\t (1 = elective, 2 = acute) \n")        
        
        cat("\n pregnancies with a current CS (ELEKAKUT) split by year: \n")        
        print(table(dat$ELEKAKUT,dat$AR,useNA="a"))
        cat("\n note that before 1999 there is no data! \n")        
        
        cat("\n exclusions will be done if SECFORE=1 (yes) or ELEKAKUT=1 (elective) \n")        
        print(table(SECFORE=dat$SECFORE,ELEKAKUT=dat$ELEKAKUT,useNA="a"))
        
        bad_rows = which((dat$ELEKAKUT==1)|(dat$SECFORE==1))
        
        cat("\n in total ",length(bad_rows)," rows will be removed \n")        
        dat = dat[-bad_rows,]
        cat("\n in total ",nrow(dat)," are remaining \n")        
        
        generate_year_counts(dat, "currentCSection", F)
        dat
}



fun_noPROM = function(dat) {
cat("PROM deliveries: \n \t exclude deliveries that start with premature rupture of membranes \n\n")

cat( "Exclusion codes based on ICD-08 system (1969-1986):
\t ^76910  - ruptura praematura membranae fetus")
cat( "Exclusion codes based on ICD-09 system (1987-1997):
\t ^658B  - för tidig hinnbristning (för tidig fostervattenavgång)
\t ^658C  - fördröjd förlossning efter spontan eller ospecificerad hinnbristning
\t ^658D -  fördröjd förlossning efter artificiell hinnsprängning ")
cat( "Exclusion codes based on ICD-10 system (1997-xxxx):
\t ^O42  -  PROM 
\t ^O755 - Fördröjd förlossning efter artificiell hinnsprängning
\t ^O756 - Fördröjd förlossning efter spontan eller icke specificerad hinnbristning (3cat: A,B,X) ")

icd08  = "^76910" #  ruptura praematura membranae fetus
icd09b = "^658B" #  för tidig hinnbristning (för tidig fostervattenavgång)
icd09c = "^658C" #  fördröjd förlossning efter spontan eller ospecificerad hinnbristning
icd09d = "^658D" #  fördröjd förlossning efter artificiell hinnsprängning ")
icd10a = "^O42"  # PROM
icd10b = "^O755" # Fördröjd förlossning efter artificiell hinnsprängning
icd10c = "^O756" # Fördröjd förlossning efter spontan eller icke specificerad hinnbristning (3cat: A,B,X)

icd_codes = c(icd08,icd09b,icd09c,icd09d,icd10a,icd10b,icd10c)    
        

tbl = rixs = cnt = NULL # table of effect sizes, exclusion rows, counts per year
        for (i in 1:length(icd_codes)) {
                icd_code = icd_codes[i]
                print(icd_code)
                
                # which rows (pregnancies) have a problem
                rix =  NULL; for(j in grep("^MDIAG|^BDIAG",colnames(dat))) rix=c(rix,grep(icd_code,dat[,j]))
                
                # year-split of incidence (cumulative)
                phe = factor(rep(0,nrow(dat)),levels = c(0,1))
                
                if(length(rix)>0) {
                phe[unique(rix)] = 1
                coefs = coef(summary(lm(dat$GRDBS ~ phe)))
                df = data.frame(condition=icd_code,N=length(rix),beta=round(coefs[2,1],1),pval=coefs[2,4])
                } else {
                df = data.frame(condition=icd_code,N=length(rix),beta=NA,pval=NA)
                }
                
                cnt = rbind(cnt, table(dat$AR,phe)[,2]) # count                
                tbl = rbind(tbl,df)
                rixs = c(rixs,rix)
                rm(df,phe,rix)
        }        
        rixs = sort(unique(rixs))
        
        cnt = as.data.frame(t(cnt))
        colnames(cnt) = icd_codes
        print(cnt)
        print(tbl)

cat("in total ",length(rixs), " rows will be removed")
        dat = dat[-rixs,]
cat(" there are ",nrow(dat), " rows remaining") 
        
        generate_year_counts(dat, "noPROM", F)
        dat
}

# other name - pregCompl
fun_ICDcodes = function(dat,icd_use) {
        # second argument is a vector of logical TRUE/FALSE values for ICD codes. example:
        # icd_use = c(T,T,T,T,T,T,T, T,T,T)
        if(length(icd_use)!=10) {
                warning("the length of icd_use is not equal to 10")
                break
        }
        
cat("\t NOTE:
\t - icd08  codes were used in 1969-1986 
\t - icd09  codes were used in 1987-1996
\t - icd10 codes were used in 1997+ 
\t - there is a significant nonsynonymous overlap between ICD-8 and ICD-9 \n")
        
cat("\n\t IMPORTANT:
\t - only icd10 & icd09 will be used -> all pregs before 1987 will be EXCLUDED!
\t - icd09 and icd10 codes do not match perfectly
\t - icd08 codes are not implemented here \n") 


# dangerous! temporary! is it worth?!
dat = dat[which(dat$AR>=1987),]

cat("\n\t YEAR-REMOVAL STAGE: 
    \t - in total",length(bad_rows),"rows were removed 
    \t - and",nrow(dat),"rows are left remaining \n") 

icd_codes = c("^O40",      # icd10: Abnormally large amount of amniotic fluid
              "^O410",     # icd10: Oligohydramnion
              "^O41[2-9]", # icd10: Other amnion and membrane problems
              "^O43",      # icd10: Placentalt transfusionssyndrom, Missbildning, Patologiskt fastsittande, accreta/increta/percreta
              "^O44",      # icd10: Placenta praevia, types
              "^O45",      # icd10: Premature separation of the placenta (Placenta abruptio)
              "^O46",      # icd10: Bleeding before delivery
              "^641",      # icd09: Placenta previa, Abruptio placentae, Hemorrhage in pregnancy
              "^657",      # icd09: Polyhydramnios
              "^658")      # icd09: Oligohydramnios, PROM, other amniotic cavity and membranes problems

icd_short = c("icd10_MuchAmnio",
              "icd10_OligAmnio",
              "icd10_AmnioMemb",
              "icd10_PlacGrowth",
              "icd10_PlacPrevia",
              "icd10_PlacAbrupt",
              "icd10_Bleeding",
              "icd09_PreviaAbruptio",
              "icd09_MuchAmnio",
              "icd09_OligAmnioPROM")

icd_long = c("Abnormally large amount of amniotic fluid",
             "Oligohydramnion",
             "Other amnion/membrane problems",
             "Abnormality in the placenta types",
             "Placenta praevia 2types",
             "Premature separation of the placenta (abruptio)",
             "Bleeding before delivery",
             "PlacentaPrevia,AbruptioPlacentae,Hemorrhage in pregnancy",
             "Polyhydramnios",
             "Oligohydramnios,PROM,oth.amnion&membr.problms")
        
# other ideas:
# "^656[0,5,7]|^7622"
# ^641  Blödning i sen graviditet, för tidig avlossning av moderkakan och T föreligg ande moderkaka
# Haemorrhagia in graviditate posteriore, abruptio placentae et placenta praevia 
# ^642  Hypertoni som komplikation till graviditet, förlossning och barnsängs (-9 d)
        
#generate a report of exclusions to be used:        
rep = data.frame(ICD_codes = icd_codes,use_filter=icd_use,short=icd_short,description=icd_long)

cat("\n These ICD10 codes will be used as exclusion criteria: \n")
print(rep)

# leave only those which were marked to be used
rep = rep[which(rep$use_filter),]

tbl = rixs = cnt = NULL # table of effect sizes, exclusion rows, counts per year
for (i in 1:nrow(rep)) {
        icd_code = rep$ICD_codes[i]
        #print(as.character(icd_code))
        
        # which rows (pregnancies) have a problem
        rix = cds = NULL  # rix = row indexes; cds = actual codes
        for(j in grep("^MDIAG|^BDIAG",colnames(dat))) {
                rix_tmp = grep(icd_code,dat[,j])
                rix=c(rix,rix_tmp)
                cds = c(cds,unique(dat[rix_tmp,j]))
                rm(rix_tmp)
        }
        # report the discovered matching patterns
        print(sort(unique(cds)))
        
        # year-split of incidence (cumulative)
        phe = rep(0,nrow(dat)); phe[unique(rix)]=1
        cnt = rbind(cnt, table(dat$AR,phe)[,2]) # count
        
        coefs = coef(summary(lm(dat$GRDBS ~ phe)))
        df = data.frame(code=rep$ICD_codes[i],condition=rep$short[i],N=length(rix),
                        beta=round(coefs[2,1],1),pval=coefs[2,4],stringsAsFactors = F)
        tbl = rbind(tbl,df)
        rixs = unique(c(rixs,rix))
        rm(df,coefs,phe,rix)
}

rixs = unique(rixs)

cnt = as.data.frame(t(cnt))
colnames(cnt) = rep$short
    
cat("ICD codes found for each year: \n")
print(cnt)

cat("ICD codes and their effect sizes on GestAge: \n")
        print(tbl)
        
cat("\n\t ICD-based exclusion stage: 
    \t - in total",length(rixs), "rows will be removed \n")
        dat = dat[-rixs,]
cat("\t - and",nrow(dat),"rows are left remaining \n") 

        generate_year_counts(dat, "ICDcodes", F)
        return(dat)
}


fun_GAmiss = function(dat) {
        cat("MISSING GESTATIONAL AGE:
            ")
        good_rows = which(!is.na(dat$GRDBS))
        cat("in total ",nrow(dat)-length(good_rows),"rows will be removed")
        dat = dat[good_rows,]
        cat("in total ",nrow(dat),"left remaining")
        
        generate_year_counts(dat, "GAmiss", F)
        dat
}

fun_GAdating = function(dat,ok_codes) {
        #ok_codes = c(1,2,5,6,7,8,10)
        
        cat("RELIABLE GESTATIONAL AGE DATING METHOD \n \n")
        cat("\t these codes will be used as reliable method indicators: \n")
        print(paste("GRMETOD = ",paste(ok_codes,collapse=", "),sep=""))
        
        bad_rows = which(! dat$GRMETOD %in% ok_codes)
        cat("in total ",length(bad_rows),"rows will be removed")
        dat = dat[-bad_rows,]
        cat("in total ",nrow(dat),"left remaining")
        
        generate_year_counts(dat, "GAdating", F)
        dat
}


fun_mHghQC = function(dat,low,upp,discordantClean,transferHeight) {
        # USAGE:
        # dat -                 data frame with columns [AR,MLANGD,lpnr_mor] in any order
        # low -                 lower threshold for excluding height values (e.g., 135 cm)
        # upp -                 upper threshold for excluding height values (e.g., 205 cm)
        # discordantClean -     logical (TRUE/FALSE), whether mothers with multiple pregnancies 
        #                       and very discordant height values between them should be swiped out
        # transferHeight -      logical (TRUE/FALSE), whether missing height values should be 
        #                       transfered from other pregnancies of the same mother (+/- 1 year)
        
        cat("REMOVE MISSING OR UNRELIABLE MATERNAL HEIGHT, TRANSFER FROM OTHER PREGS: \n")
        cat("\t - note: values with missing height will be removed \n")
        cat("\t - note: MatHgh only available for years 1982-2012, so years 1973-1981 will be lost! \n")
        cat("\t - note: should be run only AFTER removing rows with missing mother ID \n")
        
        ########
        ######## 1) DATA PREPARATION AND REFORMATING
        ########
        
        # save the original column order
        orig_clnms = colnames(dat)
        
        # save the original row-wise order
        dat$original_order = seq(nrow(dat)) # to be restored after this block is done
        
        # set missing/emtpty height values to NA
        dat$MLANGD[which(dat$MLANGD=="")]=NA # maternal height must be present for adjustments
        
        # make all values numeric
        dat$MLANGD = as.numeric(dat$MLANGD) # ensure that values are numerical
        
        # for final reporting
        n_missing_vals_start = sum(is.na(dat$MLANGD)) # for comparison at the end
        
        ########
        ######## 2) SIMPLE CLEANING OF EXTREME VALUES
        ########
        
        cat(" initial range of maternal height: ",paste(range(dat$MLANGD,na.rm=T),collapse="-"),"cm \n")
        
        sd_low = round(mean(dat$MLANGD,na.rm=T)-sd(dat$MLANGD,na.rm=T)*2,0)
        sd_upp = round(mean(dat$MLANGD,na.rm=T)+sd(dat$MLANGD,na.rm=T)*2,0)
        cat(" +/- 2SD range of maternal height: ",paste(sd_low,sd_upp,sep="-"),"cm \n")
        
        sd_low = round(mean(dat$MLANGD,na.rm=T)-sd(dat$MLANGD,na.rm=T)*3,0)
        sd_upp = round(mean(dat$MLANGD,na.rm=T)+sd(dat$MLANGD,na.rm=T)*3,0)
        cat(" +/- 3SD range of maternal height: ",paste(sd_low,sd_upp,sep="-"),"cm \n")
        
        sd_low = round(mean(dat$MLANGD,na.rm=T)-sd(dat$MLANGD,na.rm=T)*4,0)
        sd_upp = round(mean(dat$MLANGD,na.rm=T)+sd(dat$MLANGD,na.rm=T)*4,0)
        cat(" +/- 4SD range of maternal height: ",paste(sd_low,sd_upp,sep="-"),"cm \n")
        
        cat("\t selected range for inclusion:",paste(low,upp,sep="-"),"cm \n")
        
        # report
        cat("\t\t",sum(dat$MLANGD<low,na.rm=T),"[too low] height values set to missing \n")
        cat("\t\t",sum(dat$MLANGD>upp,na.rm=T),"[too high] height values set to missing \n")
        
        dat$MLANGD[which(dat$MLANGD<low)]=NA # suspicious-height threshold determined by Julius (also, no sib-similarity below this thr)
        dat$MLANGD[which(dat$MLANGD>upp)]=NA # almost incredible values
        
        ########
        ######## 3) DELETE VERY DICORDANT HEIGHTS
        ########
        
        if (discordantClean==TRUE) {
                
                # when maternal height values from different pregnancies are too discordant
                cat("\t selected option:  discordantClean=TRUE \n")
                
                # assign spec identifiers
                dat$ids_hghs = paste(dat$lpnr_mor,dat$MLANGD,sep="_")
                
                repeat {
                        sub = dat[which(!is.na(dat$MLANGD)),c("lpnr_mor","MLANGD")]
                        
                        # prepare the working dataframe
                        df = as.data.frame(group_by(sub,lpnr_mor) %>% 
                                                   summarise(n=n(),mnH=mean(MLANGD),
                                                             miH=min(MLANGD),mxH=max(MLANGD)) %>% ungroup())
                        # range of height values for each mother
                        df$dif = df$mxH - df$miH
                        
                        # first tier - when there are only two height values
                        bad_moms_1_temp = df[which( (df$n == 2) & (df$dif >= 10) ),]  # PARAMETER !
                        bad_moms_1_tempA = bad_moms_1_temp[,c("lpnr_mor","miH")]
                        bad_moms_1_tempB = bad_moms_1_temp[,c("lpnr_mor","mxH")]
                        colnames(bad_moms_1_tempA)=c("lpnr_mor","MLANGD")
                        colnames(bad_moms_1_tempB)=c("lpnr_mor","MLANGD")
                        bad_moms_1 = rbind(bad_moms_1_tempA,bad_moms_1_tempB)
                        
                        # cleanup
                        rm(bad_moms_1_temp,bad_moms_1_tempA,bad_moms_1_tempB)
                        
                        # second tier - when there are three or more values              
                        cand_mom_ids = unique(df$lpnr_mor[which((df$n>2)&(df$dif>10))])   # PARAMETER !
                        sub = sub[which(sub$lpnr_mor %in% cand_mom_ids),] # reduce to improve speed
                        
                        # find which height value to remove for each mom
                        cum = list() # cummulator
                        num = 0 # counter
                        for (id in cand_mom_ids) {
                                num = num + 1
                                vls = sub[which(sub$lpnr_mor == id),"MLANGD"]
                                vls = sort(vls[which(!is.na(vls))]) # very important to sort !
                                # estimate left and right two-value distances
                                dif_low = vls[2] - vls[1]
                                dif_upp = vls[length(vls)] - vls[length(vls)-1]
                                min_val = vls[1]
                                max_val = vls[length(vls)]
                                nm_vals = length(vls)
                                cum[[num]] = data.frame(id,nm_vals,min_val,max_val,dif_low,dif_upp,stringsAsFactors = F)
                                rm(vls,dif_low,dif_upp,min_val,max_val,nm_vals)
                        }
                        
                        cum = do.call(rbind.data.frame, cum) # combine list to dataframe
                        
                        # cleanup
                        rm(df,id,cand_mom_ids,num)
                        
                        if (nrow(cum)>0) {
                                bad_moms_2 = cum[which(cum$dif_low > cum$dif_upp),c("id","min_val")]
                                bad_moms_3 = cum[which(cum$dif_low < cum$dif_upp),c("id","max_val")]
                                colnames(bad_moms_2)=c("lpnr_mor","MLANGD")
                                colnames(bad_moms_3)=c("lpnr_mor","MLANGD")
                        } else {
                                bad_moms_2 = bad_moms_3 = NULL
                        }
                        
                        cand_bad_ids = unique(cum$id[which(cum$dif_low == cum$dif_upp)])
                        bad_moms_4 = sub[which(sub$lpnr_mor %in% cand_bad_ids),]
                        colnames(bad_moms_4)=c("lpnr_mor","MLANGD")
                        
                        # combine
                        bad_moms_all = rbind(bad_moms_1,bad_moms_2,bad_moms_3,bad_moms_4)
                        bad_ids_hghs = unique(paste(bad_moms_all$lpnr_mor,bad_moms_all$MLANGD,sep="_"))
                        
                        
                        bad_hgh_rix = which(dat$ids_hghs %in% bad_ids_hghs)
                        
                        
                        if (length(bad_hgh_rix)==0) break
                        
                        # set to missing
                        dat$MLANGD[bad_hgh_rix] = NA  # set to missing (as unreliable)
                        
                        # report
                        cat("\t\t",length(bad_hgh_rix),"height values set to missing \n")
                        # cleanup
                        rm(bad_hgh_rix,bad_ids_hghs,cand_bad_ids,sub,cum)
                        rm(bad_moms_1,bad_moms_2,bad_moms_3,bad_moms_4,bad_moms_all)
                        
                }
        }     
        
        ########
        ######## 3) TRANSFER VALUES FROM OTHER PREGNENCIES
        ########
        
        # when mother has many pregnancies, at some pregnancies she might have missing height data.
        # in such instances, height values could be transfered from her other pregnancies.
        
        # transferHeight function
        if(transferHeight==TRUE) {
                cat("\t selected option:  transferHeight=TRUE \n")
                # note: height transfer will not be done for birth years <1982 even though it is possible
                
                # unique identifier of mother and her pregnancy year
                dat$midAR = paste(dat$lpnr_mor,dat$AR,sep="_")
                
                # which years have absolutely no height data
                tbl = as.matrix(table(is.na(dat$MLANGD),dat$AR))
                noHgh_years = as.numeric(names(which(tbl[1,]==0)))
                
                # two types of rows of interest (source of data, and the target)
                filed_height_rix = which( (!is.na(dat$MLANGD))&(!dat$AR %in% noHgh_years)) # with height in good years
                empty_height_rix = which( (is.na(dat$MLANGD))&(!dat$AR %in% noHgh_years)&
                                                  (dat$lpnr_mor %in% dat$lpnr_mor[filed_height_rix])) # no height in good years with potential to restore
                filed_height_rix = which( (!is.na(dat$MLANGD))&(!dat$AR %in% noHgh_years)&
                                                  (dat$lpnr_mor %in% dat$lpnr_mor[empty_height_rix])) # refine: only useful
                
                
                # rows
                wrk_rix = c(filed_height_rix, empty_height_rix) # working rows used in manipulations
                oth_rix = seq(nrow(dat))[which(! seq(nrow(dat)) %in% wrk_rix)] # other rows
                # datasets
                wrk_dat = dat[wrk_rix,] # data that WILL be used in manipulations
                oth_dat = dat[oth_rix,] # data that will NOT be used in manipulations
                rm(wrk_rix,oth_rix)
                
                # extract two subsets (target and source) to facilitate the running time
                noH = dat[empty_height_rix,c("AR","lpnr_mor")]  # no height
                isH = dat[filed_height_rix,c("AR","lpnr_mor","MLANGD")] # with height
                #sum(is.na(isH$MLANGD)) # should be 0
                
                # cleanup
                rm(filed_height_rix,empty_height_rix,noHgh_years,tbl)
                
                fun_transf = function(wrk_dat,noH,isH,ofs,ok_dif) {
                        # create a list with recovered maternal heights
                        yrs = sort(unique(noH$AR),decreasing = T) # years of interest
                        #ofs = ofs # mother with plus/minus this number of years is assumed to have the same height)
                        #ok_dif = ok_dif # allowed difference between two heights of the same mother
                        lst = list() # accumulator of imputed (transfered) height values
                        cnt = 0 # counter
                        for (yr in yrs) {
                                cnt = cnt + 1 # reset the counter
                                noH_momIDs = unique(noH$lpnr_mor[which(noH$AR==yr)]) # these moms have no height values
                                # subset a dataframe that contains the same mothers of similar age, but WITH height data:
                                tmp = isH[which((isH$AR %in% (yr-ofs):(yr+ofs))&(isH$lpnr_mor %in% noH_momIDs)),]
                                df = group_by(tmp,lpnr_mor) %>% summarise(n=n(),mnH=mean(MLANGD),miH=min(MLANGD),mxH=max(MLANGD)) %>% ungroup()
                                #table(df$n) # how many height values how many mothers have
                                df$mnH = round(df$mnH,0) # to avoid non-integer values
                                df$dif = df$mxH - df$miH # difference in max and min reported heights
                                #table(df$dif) # the number of mothers with specific range of height values
                                #sum(is.na(df$mnH)) # should be 0
                                df$mnH[which(df$dif>ok_dif)] = NA # unreliable values
                                if(nrow(df)>0) {
                                        df = data.frame(midAR=paste(df$lpnr_mor,yr,sep="_"),mnH=df$mnH,stringsAsFactors = F)
                                        lst[[cnt]] = df
                                }
                                rm(noH_momIDs,tmp,df)
                        }
                        rm(cnt,ofs,yrs,yr)
                        newH = do.call(rbind.data.frame, lst) # combine list to dataframe
                        tmp = merge(wrk_dat,newH,by="midAR",all.x=T) 
                        # rows where original missing height should be replaced with height from other pregnancies
                        repl_rix = which((is.na(tmp$MLANGD))&(!is.na(tmp$mnH))) #  replacement rows
                        tmp[repl_rix,"MLANGD"] = tmp[repl_rix,"mnH"] # replace missing values
                        tmp = tmp[order(tmp$original_order),] # restore order
                        tmp = tmp[,-grep("mnH",colnames(tmp))]
                        cat("\t\t",length(repl_rix),"missing height values recovered \n")
                        return(tmp)
                        rm(tmp,newH,repl_rix)
                } # end of inner function
                
                # transfer values using different offsets. couple of times!
                wrk_dat = fun_transf(wrk_dat,noH,isH,ofs=5,ok_dif=2)  # settings!
                wrk_dat = fun_transf(wrk_dat,noH,isH,ofs=0,ok_dif=2)  # settings!
                
                # possible parameter settings (for offest in birth year of the child):
                # offset = 0,  recovered height values (rows) = 1281
                # offset = 1,  recovered height values (rows) = 32628
                # offset = 2,  recovered height values (rows) = 124431
                # offset = 3,  recovered height values (rows) = 187935
                # offset = 4,  recovered height values (rows) = 225131
                # offset = 5,  recovered height values (rows) = 247558
                # offset =20,  recovered height values (rows) = 292315
                # but avoid using very large offest, not to transfer height of very young females to very tall ones
                #if ofs=20 and then ofs=0 ->  52 more values are recovered
                
                # reassemble the original dataframe with transfered height values
                dat = rbind(wrk_dat,oth_dat)
                
                # cleanup
                rm(wrk_dat,oth_dat,noH,isH)
                
        } # end of if
        
        n_missing_vals_end = sum(is.na(dat$MLANGD)) # for comparison 
        n_difference = n_missing_vals_start - n_missing_vals_end
        gained_removed = ifelse(n_difference>0,"gained","removed")
        cat(" - in total,",n_difference,"height values will be",gained_removed," \n")
        dat = dat[which(!is.na(dat$MLANGD)),]
        cat(" - in total,",nrow(dat),"rows are left \n")
        
        generate_year_counts(dat, "mHghQC", F)
        
        # restore original order
        dat = dat[order(dat$original_order),]
        dat = dat[,orig_clnms]
        rm(orig_clnms)
        return(dat)
}








fun_mWghQC = function(dat,duringPregTransf) {
        # USAGE:
        # dat -  data frame with columns [AR,MVIKT,MVIKTFV,lpnr_mor] in any order
        # duringPregTransf - should MVIKTFV values be used to guess missing MVIKT
        
        cat("TRANSFER MISSING MATERNAL WEIGHT FROM OTHER PREGS (and cols): \n")
        cat("\t - note: values with missing weight will NOT be removed \n")
        cat("\t - note: MatWgh not available for years 1990-1991 \n")
        cat("\t - note: should be run only AFTER removing rows with missing mother ID \n")
        
        ########
        ######## 1) META
        ########
        
        # save the original column order
        orig_clnms = colnames(dat)
        
        # save the original row-wise order
        dat$original_order = seq(nrow(dat)) # to be restored after this block is done
        
        # for final reporting
        n_missing_vals_start = sum(is.na(dat$MVIKT)) # for comparison at the end
        
        ########
        ######## 2) TRANSFER VALUES FROM OTHER PREGNENCIES
        ########
        
        # when mother has many pregnancies, at some pregnancies she might have missing weight data.
        # in such instances, weight values could be transfered from her other pregnancies.
        
        # unique identifier of mother and her pregnancy year
        dat$midAR = paste(dat$lpnr_mor,dat$AR,sep="_")
        
        # two types of rows of interest (source of data, and the target)
        filed_weight_rix = which(!is.na(dat$MVIKT)) # with height in good years
        empty_weight_rix = which( (is.na(dat$MVIKT))&(dat$lpnr_mor %in% dat$lpnr_mor[filed_weight_rix])) # no weight but with potential to restore
        filed_weight_rix = which( (!is.na(dat$MVIKT))&(dat$lpnr_mor %in% dat$lpnr_mor[empty_weight_rix])) # refine: only useful
        
        # rows
        wrk_rix = c(filed_weight_rix, empty_weight_rix) # working rows used in manipulations
        oth_rix = seq(nrow(dat))[which(! seq(nrow(dat)) %in% wrk_rix)] # other rows
        # datasets
        wrk_dat = dat[wrk_rix,] # data that WILL be used in manipulations
        oth_dat = dat[oth_rix,] # data that will NOT be used in manipulations
        rm(wrk_rix,oth_rix)
        
        # extract two subsets (target and source) to facilitate the running time
        noW = dat[empty_weight_rix,c("AR","lpnr_mor")]  # no weight
        isW = dat[filed_weight_rix,c("AR","lpnr_mor","MVIKT")] # with weight
        #sum(is.na(isW$MVIKT)) # should be 0
        
        # cleanup
        rm(filed_weight_rix,empty_weight_rix)
        
        fun_transf = function(wrk_dat,noW,isW,ofs,ok_dif) {
                # create a list with recovered maternal weights
                yrs = sort(unique(noW$AR),decreasing = T) # years of interest
                #ofs = ofs # mother with plus/minus this number of years is assumed to have the same weight)
                #ok_dif = ok_dif # allowed difference between two weights of the same mother
                lst = list() # accumulator of imputed (transfered) height values
                cnt = 0 # counter
                for (yr in yrs) {
                        cnt = cnt + 1 # reset the counter
                        noW_momIDs = unique(noW$lpnr_mor[which(noW$AR==yr)]) # these moms have no weight values
                        # subset a dataframe that contains the same mothers of similar age, but WITH weight data:
                        tmp = isW[which((isW$AR %in% (yr-ofs):(yr+ofs))&(isW$lpnr_mor %in% noW_momIDs)),]
                        df = group_by(tmp,lpnr_mor) %>% summarise(n=n(),mnW=mean(MVIKT),miW=min(MVIKT),mxW=max(MVIKT)) %>% ungroup()
                        #table(df$n) # how many weight values how many mothers have
                        df$mnW = round(df$mnW,0) # to avoid non-integer values
                        df$dif = df$mxW - df$miW # difference in max and min reported weights
                        #table(df$dif) # the number of mothers with specific range of weight values
                        #sum(is.na(df$mnW)) # should be 0
                        df$mnW[which(df$dif>ok_dif)] = NA # unreliable values
                        if(nrow(df)>0) {
                                df = data.frame(midAR=paste(df$lpnr_mor,yr,sep="_"),mnW=df$mnW,stringsAsFactors = F)
                                lst[[cnt]] = df
                        }
                        rm(noW_momIDs,tmp,df)
                }
                rm(cnt,ofs,yrs,yr)
                newW = do.call(rbind.data.frame, lst) # combine list to dataframe
                tmp = merge(wrk_dat,newW,by="midAR",all.x=T) 
                # rows where original missing weight should be replaced with weight from other pregnancies
                repl_rix = which((is.na(tmp$MVIKT))&(!is.na(tmp$mnW))) #  replacement rows
                tmp[repl_rix,"MVIKT"] = tmp[repl_rix,"mnW"] # replace missing values
                tmp = tmp[order(tmp$original_order),] # restore order
                tmp = tmp[,-grep("mnW",colnames(tmp))]
                cat("\t\t",length(repl_rix),"missing weight values recovered \n")
                return(tmp)
                rm(tmp,newW,repl_rix)
        } # end of inner function
        
        # transfer values using different offsets. couple of times!
        wrk_dat = fun_transf(wrk_dat,noW,isW,ofs=5,ok_dif=10)  # settings!
        wrk_dat = fun_transf(wrk_dat,noW,isW,ofs=0,ok_dif=5)  # settings!
        
        # reassemble the original dataframe with transfered height values
        dat = rbind(wrk_dat,oth_dat)
        
        # cleanup
        rm(wrk_dat,oth_dat,noW,isW)
        
        n_missing_vals_end = sum(is.na(dat$MVIKT)) # for comparison 
        n_difference = n_missing_vals_start - n_missing_vals_end
        cat(" - in total,",n_difference,"weight values will be recovered \n")
        
        
        if (duringPregTransf==TRUE) {
                
                cat("\t selected option:  duringPregTransf=TRUE \n")
                #note: 1982-1989 MVIKTF values vere truncated to be <100 kg
                
                fun_hexbplot = function(X,Y,LOG,XLAB,YLAB) {
                        h=hexbin(Y~X); x=h@xcm; y=h@ycm; s=h@count
                        plot(x,y,pch=1,cex=log(s,LOG)+0.2,xlab=XLAB,ylab=YLAB,
                             xlim=range(X,na.rm=T),ylim=range(Y,na.rm=T))
                }
                fun_hexbplot(dat$MVIKT,dat$MVIKTFV,200,"matWgh_beforePreg (MVIKT)",
                             "matWgh_duringPreg (MVIKTFV)")
                abline(-15,1); abline(50,1)
                
                # temporary modify the dataset 
                dat$par3fctr = dat$PARITET
                dat$par3fctr[which(dat$par3fctr>3)] = 3
                dat$par3fctr = as.factor(dat$par3fctr)
                dat$twinsing = dat$BORDF2
                dat$twinsing = as.factor(dat$twinsing)
                
                # model on the dataset without obvious outliers
                full_ix = which( (dat$MVIKTFV<(dat$MVIKT+50))&(dat$MVIKTFV>(dat$MVIKT-15)))
                sub = dat[full_ix,]
                mod = lm(MVIKT ~ poly(MVIKTFV,2) + poly(AR,2) + poly(MALDER,2) + par3fctr + twinsing, data = sub)
        
                # split original data into "full" and "nonfull"        
        rix = which((!is.na(dat$PARITET))&(!is.na(dat$AR))&(!is.na(dat$MVIKTFV))&(!is.na(dat$MALDER))&(!is.na(dat$BORDF2)))
        part1 = dat[rix,]
        part2 = dat[-rix,]
        
        # predict MVIKT values where possible, return to original data structure
        ys = as.numeric(predict(mod,newdata = part1))
        part1$pred_MVIKT = ys
        part2$pred_MVIKT = NA
        dat = rbind(part1,part2)
                
        fun_hexbplot(dat$MVIKT,dat$pred_MVIKT,200,"matWgh_beforePreg_real","matWgh_beforePreg_predicted")
        abline(0,1)
        

        n_vals_gained = sum((is.na(dat$MVIKT))&(!is.na(dat$pred_MVIKT)))
        vals_gained_rix = which((is.na(dat$MVIKT))&(!is.na(dat$pred_MVIKT)))
        cat("\t\t",n_vals_gained,"maternal weight values were recovered using MVIKTFV \n")
        
        hist(dat$pred_MVIKT[vals_gained_rix],breaks=100,col="grey",
             xlab="new (recovered) MVIKT values",main="")
        
        dat$MVIKT[vals_gained_rix] = as.integer(round(dat$pred_MVIKT[vals_gained_rix],0))
        
        } # end of if
        
        
        #generate_year_counts(dat, "mWghQC", F)  # this was only a helper cleaning, thus no need for reporting
        
        # restore original order
        dat = dat[order(dat$original_order),]
        dat = dat[,orig_clnms]
        rm(orig_clnms)
        return(dat)
}

head(dat)


fun_mBmiQC = function(dat) {
        cat("MATERNAL HEIGHT-WEIGHT BIVARIATE ANOMALIES (BMI): \n")
        cat("\t note: ",round(mean(is.na(dat$MVIKT))*100),"% of WEIGHT values are missing \n",sep="")
        
        ########
        ######## 1) META
        ########
        
        # save the original column order
        orig_clnms = colnames(dat)
        
        # save the original row-wise order
        dat$original_order = seq(nrow(dat)) # to be restored after this block is done
        
        # a function which plots bivariate plot without the time consumtion
        fun_hexbplot = function(X,Y,LOG,XLAB,YLAB) {
                h=hexbin(Y~X); x=h@xcm; y=h@ycm; s=h@count
                plot(x,y,pch=1,cex=log(s,LOG)+0.2,xlab=XLAB,ylab=YLAB,
                     xlim=range(X,na.rm=T),ylim=range(Y,na.rm=T))
        }
        
        ########
        ######## 2) ESTIMATE BMI AND ITS THRESHOLDS
        ########
        
        # BMI
        dat$bmi = dat$MVIKT / (dat$MLANGD / 100)^2
        
        # TRESHOLDS BY AGE
        ages = 14:48  # ages with reasonable numbers of individuals
        prbs = c(0.01,0.05,seq(0.1,0.9,0.1),0.95,0.99) # percentiles
        qnts = matrix(NA,nr=length(ages),nc=length(prbs))
        for (i in 1:length(ages)) {
                age = ages[i]
                bmis = dat$bmi[which(dat$MALDER==age)]
                qnts[i,] = as.numeric(quantile(bmis,probs = prbs,na.rm = T))
        }
        
        # VISUALIZE
        plot(NA,xlim=range(ages),ylim=c(min(qnts)-2,max(qnts)+2),
             xlab="maternal age",ylab="maternal BMI",main="determine thresholds for BMI")
        for (i in 1:length(prbs)) points(qnts[,i]~ages,type="l",col="grey")
        text(x = rep(min(ages),length(prbs)),y = qnts[1,],labels = prbs,pos = 4,cex = 0.5)
        text(x = rep(max(ages),length(prbs)),y = qnts[nrow(qnts),],labels = prbs,pos = 2,cex = 0.5)
        
        # smooth extreme quantiles
        mod_upp = loess(qnts[,length(prbs)]~ages,span = 0.5) # 99 percentile
        mod_low = loess(qnts[,1]~ages,span = 0.5) # 1 percentile
        df = data.frame(ages,stringsAsFactors = F)
        df$bmi_upp = predict(mod_upp,df) 
        df$bmi_low = predict(mod_low,df) 
        # draw smoothed percentiles 1 and 99
        points(df$bmi_upp ~ df$ages,type="l",col="red")
        points(df$bmi_low ~ df$ages,type="l",col="red")
        
        # add offset
        df$bmi_upp = df$bmi_upp + 2
        df$bmi_low = df$bmi_low - 2
        # draw smoothed percentiles 1 and 99  with offset
        points(df$bmi_upp ~ df$ages,type="l",col="blue")
        points(df$bmi_low ~ df$ages,type="l",col="blue")
        
        # add thresholds to the data-frame
        dat = merge(dat,df,by.x="MALDER",by.y="ages",all.x=T)
        
        # note, that very extreme maternal age will get no threshold, thus will not be removed
        
        
        ########
        ######## 3) CLEANING
        ########
        
        fun_hexbplot(dat$MLANGD,dat$MVIKT,200,"Maternal height","Maternal weight")
        
        # hard-set bmi threshold (ONLY FOR COMPARISON)
        high_bmi_ix = which(dat$bmi>40)
        low_bmi_ix = which(dat$bmi<15)
        points(MVIKT~MLANGD,data=dat[high_bmi_ix,],col="orange",pch=19,cex=1)
        points(MVIKT~MLANGD,data=dat[low_bmi_ix,],col="cornflowerblue",pch=19,cex=1)
        rm(high_bmi_ix,low_bmi_ix)
        
        # mixup mistakes (height = weight)
        bad0 = which(dat$MLANGD == dat$MVIKT)
        points(MVIKT~MLANGD,data=dat[bad0,],col="green",pch=19,cex=1.2)
        
        # hardset outlier thresholds
        
        abline(-475,4,col="blue")
        bad1=which(  dat$MVIKT > (dat$MLANGD*4 - 475))
        points(MVIKT~MLANGD,data=dat[bad1,],col="red",pch=19,cex=0.7)

        abline(1100,-5,col="blue")
        bad2=which(  dat$MVIKT > (dat$MLANGD*(-5) + 1100))
        points(MVIKT~MLANGD,data=dat[bad2,],col="red",pch=19,cex=0.7)
        
        abline(-730,4,col="blue")
        bad3=which(  dat$MVIKT < (dat$MLANGD*4 - 730))
        points(MVIKT~MLANGD,data=dat[bad3,],col="red",pch=19,cex=0.7)
        
        ## age-dependent BMI thresholds and outliers
        bad_upp = which(dat$bmi>dat$bmi_upp)
        bad_low = which(dat$bmi<dat$bmi_low)
        points(y=dat$MVIKT[bad_upp],x=dat$MLANGD[bad_upp],pch=19,cex=0.3)
        points(y=dat$MVIKT[bad_low],x=dat$MLANGD[bad_low],pch=19,cex=0.3)
        
        # acummulate all bad rows
        bad_row_indx = unique(c(bad0,bad1,bad2,bad3,bad_upp,bad_low))
        
        # remove bad rows
        cat("\t",length(bad_row_indx),"rows will be removed \n")
        if(length(bad_row_indx)>0) dat=dat[-bad_row_indx,]
        cat("\n rows remaining: ", nrow(dat))        
        
        # report
        generate_year_counts(dat, "mBmiQC", F)
        
        # restore original order
        dat = dat[order(dat$original_order),]
        dat = dat[,orig_clnms]
        rm(orig_clnms)
        return(dat)
}

# other ways to identify overweight / low weight
#F50 R630 - anorexy
#E66 - obesity
#sub = dat[grep("^E66",dat$MDIAG2),c("MLANGD","MVIKT")]
#points(x=sub$MLANGD,y=sub$MVIKT,pch=19,col="blue")


fun_origin = function(dat,orig,countryBlocks,parentVariabl,strict) {
        # orig = original uncleaned dataset
        # strict = LOGICAL. if TRUE - check other pregs of the same mother/father to be consistent
        # example values
        # countryBlocks="nordic"
        # parentVariabl=c("MFODLAND","MNAT","FNAT")
        
        # possible country blocks and their countries
        sweden = "SVERIGE"
        nordic = c("SVERIGE","FINLAND","NORGE","DANMARK")
        europe = c(nordic,"POLEN","TYSKLAND","FRANKRIKE","ESTLAND","ISLAND","UKRAINA","SPANIEN","GREKLAND")
        
        # which one was selected?
        countries = get(countryBlocks)
        
        # report
        cat("PARENTAL ORIGIN: \n")
        cat("\t selected countries:",countries,"\n")
        cat("\t selected variables:",parentVariabl,"\n")
        cat("\t (note, that strictly ALL selected variables must match selected countries)\n")
        cat("\t (if no value exists for at least one variable - pregnancy gets excluded!)\n")
        
        if (length(parentVariabl)>1) {
                good_rix = which(apply(dat[,parentVariabl],1,function(x) all(x %in% countries)&all(!is.na(x))))
        }
        
        if (length(parentVariabl)==1) {
                x = dat[,parentVariabl]
                good_rix = which( (x %in% countries) & (!is.na(x)) )
                rm(x)
        }
        
        cat("\t in total",length(good_rix),"rows will remain \n")
        cat("\t in total",nrow(dat)-length(good_rix),"rows will be removed \n")
        
        dat = dat[good_rix,]
        
        if (length(parentVariabl)>1) {
                cat("\n \t report on concordance: \n")
                cmbn = combn(length(parentVariabl),2)
                for (j in 1:ncol(cmbn)) {
                        var1 = parentVariabl[cmbn[1,j]]
                        var2 = parentVariabl[cmbn[2,j]]
                        cat("\n \t rows=",var1,", columns=",var2," \n")
                        print(table(dat[,var1],dat[,var2]))
                        rm(var1,var2)
                }
                rm(cmbn)
        }
        
        ### extra precaution: what if mother in her other pregnancy reported diffrent nationality?
        if (strict) {
                #colnms = parentVariabl[which(parentVariabl %in% c("MFODLAND","MNAT"))]
                tmp = orig[which(orig$lpnr_mor %in% dat$lpnr_mor ),c("lpnr_mor","MFODLAND","MNAT")]
                ix1 = which((!tmp$MFODLAND %in% countries)&(tmp$MFODLAND!=""))
                ix2 = which((!tmp$MNAT %in% countries)&(tmp$MNAT!=""))
                if ( all(c("MFODLAND","MNAT") %in% parentVariabl) ) {
                        bad_ids = unique(c(tmp$lpnr_mor[ix1],tmp$lpnr_mor[ix2])) # refers to orginal object
                        bad_rows = which(dat$lpnr_mor %in% bad_ids) # refers to the object being cleaned
                        cat("\n ALSO, due to the \"strict\" filter flag,",length(bad_rows),"pregs will be excluded \n")
                        if (length(bad_rows)>0) dat = dat[-bad_rows,]         
                        rm(bad_ids,bad_rows)
                }
                if ( (sum(parentVariabl=="MFODLAND")==1) & (sum(parentVariabl=="MNAT")==0) )  {
                        bad_ids = unique(tmp$lpnr_mor[ix1]) # refers to orginal object
                        bad_rows = which(dat$lpnr_mor %in% bad_ids) # refers to the object being cleaned
                        cat("\n ALSO, due to the \"strict\" filter flag,",length(bad_rows),"pregs will be excluded \n")
                        if (length(bad_rows)>0) dat = dat[-bad_rows,]
                        rm(bad_ids,bad_rows)
                }
                if ( (sum(parentVariabl=="MNAT")==1) & (sum(parentVariabl=="MFODLAND")==0) )  {
                        bad_ids = unique(tmp$lpnr_mor[ix2]) # refers to orginal object
                        bad_rows = which(dat$lpnr_mor %in% bad_ids) # refers to the object being cleaned
                        cat("\n ALSO, due to the \"strict\" filter flag,",length(bad_rows),"pregs will be excluded \n")
                        if (length(bad_rows)>0) dat = dat[-bad_rows,]
                        rm(bad_ids,bad_rows)
                }
                rm(tmp,ix1,ix2)
        }
        rm(countries,good_rix)
        generate_year_counts(dat, "parent_origin", F)
        return(dat)
} 

fun_visualize_exclusions_by_year = function(year_matrix) {
        year_changes = arrange(year_matrix, AR) %>%
                mutate(rowsb = lag(rows), change = rows/rowsb) %>%
                filter(stage!="initial", change>0)
        
        breaks = c(0,0.1,0.3,0.5,0.7,0.8,  0.9,0.91,0.92,0.93,0.94 ,0.95,0.96,0.97,0.98,0.99,1)
        n_cols = length(breaks)-1
        n_blus = floor(n_cols/3)
        n_grns = n_blus
        n_reds = n_cols - n_blus - n_grns
        
        blus = brewer.pal(n_blus+2,"Blues")[2:(n_blus+1)]
        grns = rev(brewer.pal(n_grns+2, "Greens")[2:(n_grns+1)])
        reds = rev(brewer.pal(n_reds+2, "Reds")[2:(n_reds+1)])
        colrs = c(reds,grns,blus)
        
        process_labels = function(l){
                ll = regmatches(l, gregexpr("[0-9.]+", l))
                ll = sapply(ll, function(x) paste(as.numeric(x)*100, collapse="-"))
                paste(ll, "%")
        }
        
        print(ggplot(year_changes, aes(x=AR, y = factor(stage, level=unique(stage)))) +
                geom_point(size=5, aes(col = cut(change, breaks=breaks))) +
                scale_color_manual(values = colrs, guide = guide_legend(title="% of samples \n passing this step"),
                                   labels = process_labels) +
                ylab("step") + xlab("year") + theme_bw())
        
        year_matrix = mutate(year_matrix, AR = as.factor(AR), stage=as.factor(stage)) %>%
                complete(AR, stage, fill = list(rows=0))
        
        year_changes = group_by(year_matrix, AR) %>%
                summarize(final = min(rows), initial = max(rows)) %>%
                mutate(loss = paste(round((1-final/initial)*100,1), "%"))
        print.data.frame(year_changes)
}





# to consider:
# instead of removing preconditions/ICDcodes - add specific number of gest days

# to be done: 
#table(dat$ROK0)
#table(dat$ROK1)
#table(dat$ROK2)


# also add ICD-code based diagnosis... (needs a separate file with ICD codes)
#cond = read.table("~/Dropbox/PHD COURSES/2016/MEDICAL STATISTICS 2/Assignment/Project/pregnancy_complications_temporay_selection.txt",sep="\t",dec=",",h=T,stringsAsFactors = F)
#cond = cond[which(cond$decision != ""),]
# below - not exactly correct, time point of diagnossi not considered, icd file composed hastily
#icd_dat = read.table("~/Biostuff/SEMFR_DATA/diagICD_firstDate_aliIDs_verenaIDs_160512julius.csv",h=T,stringsAsFactors = F,sep=",")
#icd_dat$diag_3char = substr(icd_dat$diag,1,3)
#mom_ids = unique(icd_dat$id_verena[which(icd_dat$diag_3char %in% cond$dgn)]) # mothers to exclude
#sum(dat$lpnr_mor %in% mom_ids)  # marked fro exclusion (n=65 410)
#dat = dat[which(!dat$lpnr_mor %in% mom_ids),]
#dim(dat) # remaining 3 586 909
#rm(tbl,sub,icd_dat,bad_rows,colname,mom_ids,precond,preconditions,cond,coefs)

#.... congenital malformations, pre-eclampsia, 
#.....  not done yet




#######   CLEANING STEP 6:  UNRELIABLE GEST AGE
#df = group_by(dat,AR) %>% summarise(n=n(), mnga = mean(GRDBS,na.rm=T)) %>% ungroup()
#plot(df$mnga ~ df$AR,ylim = c(min(df$mnga,na.rm=T),max(df$mnga,na.rm=T)+1),
#     xlab="Birth year",ylab="Child's Mean Gestational Age at Birth (days)")
#mod = loess(df$mnga ~ df$AR)
#points(predict(mod,data.frame(AR=df$AR)) ~ df$AR,type="l")
#ys = rep(max(df$mnga,na.rm=T)+1,length(df$AR))
#text(df$AR,ys,round(log(df$n,10),0),cex=0.6)
#rm(mod,ys,df)

#######   CLEANING STEP 8:  MATERNAL AGE
#table(dat$MALDER)
#sum(!dat$MALDER %in% 14:45)
#sum(dat$MALDER %in% 14:45)
#dat = dat[which(dat$MALDER %in% 14:45),]
#dim(dat)
#rm(bad_rows,bad_rows_2,bad_rows_1,mod,ys,df)

#############################################
###########################   PREVIEW DATA
#############################################
#df = group_by(dat,AR) %>% summarise(n=n(), mnga = mean(GRDBS,na.rm=T)) %>% ungroup()
#plot(df$mnga ~ df$AR,ylim = c(min(df$mnga,na.rm=T),max(df$mnga,na.rm=T)+1),
#     xlab="Birth year",ylab="Child's Mean Gestational Age at Birth (days)")
#mod = loess(df$mnga ~ df$AR)
#points(predict(mod,data.frame(AR=df$AR)) ~ df$AR,type="l")
#ys = rep(max(df$mnga,na.rm=T)+1,length(df$AR))
#text(df$AR,ys,round(log(df$n,10),0),cex=0.6)
#rm(mod,ys,df)

##tmp = dat # preserve original
#dat = dat[which(dat$MFODLAND=="SVERIGE"),]
#dat = tmp; rm(tmp) # restore original

#### Maternal Height for each BirthYear
#df = group_by(dat,AR) %>% summarise(n=n(), mnhgh = mean(MLANGD)) %>% ungroup()
#plot(df$mnhgh ~ df$AR,ylim = c(min(df$mnhgh,na.rm=T),max(df$mnhgh,na.rm=T)+1),
#     xlab="Pregancy year",ylab="Mean of Maternal Height (cm)")
#mod = loess(df$mnhgh ~ df$AR)
#points(predict(mod,data.frame(AR=df$AR)) ~ df$AR,type="l")
#ys = rep(max(df$mnhgh,na.rm=T)+1,length(df$AR))
#text(df$AR,ys,round(log(df$n,10),0),cex=0.6)
#rm(mod,ys,df)

#tmp = dat # preserve original
#dat = dat[which(dat$PARITET==1),]
#dat = tmp; rm(tmp)  # restore original 

#### Maternal AGE for each BirthYear
#df = group_by(dat,AR) %>% summarise(n=n(), mnage = mean(MALDER)) %>% ungroup()
#plot(df$mnage ~ df$AR,ylim = c(min(df$mnage,na.rm=T),max(df$mnage,na.rm=T)+1),
#     xlab="Pregnancy year",ylab="Mean of Maternal Age (years)")
#mod = loess(df$mnage ~ df$AR,span = 0.3)
#points(predict(mod,data.frame(AR=df$AR)) ~ df$AR,type="l")
#ys = rep(max(df$mnage,na.rm=T)+1,length(df$AR))
#text(df$AR,ys,round(log(df$n,10),0),cex=0.6)
#rm(mod,ys,df)

#table(dat$KON)
#table(dat$PARITET)
#table(dat$MALDER)


#1973 - 1986]  ICD8
#1987 - 1996  ICD9
#1997 -  ICD10

