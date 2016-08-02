

##############   CLEAN THE SEMFR DATA   ###################
# modules that do some part of the data cleaning procedure
# by Jonas B. 2016 July


library(dplyr)
library(RColorBrewer)

# this function documents how the exclusions affected each year of the register
generate_year_counts = function(dat,show) {
        if(exists("year_matrix")) {
                tbl = table(dat$AR)
                tmp = data.frame(year=names(tbl),rows=as.numeric(tbl),stringsAsFactors = F)
                tbl = table(dat$AR)
                year_matrix = merge(year_matrix,tmp,by="year",all=T)
                year_matrix = year_matrix[order(year_matrix$year),]
                colnames(year_matrix)[-1] = paste("step",seq(length(colnames(year_matrix))-1),sep="")
        } else {
                tbl = table(dat$AR)
                year_matrix = data.frame(year=names(tbl),rows=as.numeric(tbl),stringsAsFactors = F)       
        }
        year_matrix <<- year_matrix # global assignment
        if(show==TRUE) return(year_matrix)
}
# run this function once BEFORE before cleaning the data, and then after each cleaning step
# if(exists("year_matrix")==TRUE) rm(year_matrix)
# generate_year_counts(dat,show=T)


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

# add  PREVIOUS C-SECTION
table(dat$TSECTIO,useNA="a")
precond = rep(0,nrow(dat))
precond[which(dat$TSECTIO==1)] = 1  # note, some previous CS might still remain (missing data)
coefs = coef(summary(lm(dat$GRDBS ~ precond))); coefs
df = data.frame(condition="TSECTIO",N=sum(precond==1),beta=coefs[2,1],pval=coefs[2,4])
tbl = rbind(tbl,df); rm(df)

cat("
    pregnancies with the following maternal medical conditions will be excluded from the MFR:
    ")
print(tbl)
cat("(beta and pval are estimated as an effect on gestational age)
    ")



# remove all pregnancies with preconditions from the data
sub = dat[,preconditions]; sub$prevCS = NA; sub$prevCS[which(precond==1)] = 1
bad_rows = apply(sub,1,function(x) sum(!is.na(x))); table(bad_rows); sum(bad_rows>0)

cat("
    number of pregnancies with a specific number of maternal medical conditions:
    ")
print(table(bad_rows))

dat = dat[which(bad_rows==0),]; dim(dat) 

cat("
number of rows selected for exclusion:", sum(bad_rows>0),"
number of rows remaining:",nrow(dat))

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
        
        dat[good_rix,]
}

fun_CSections = function(dat) {
        cat("CAESAREAN SECTION: \n \t additional exclusion for PREVIOUS CS and CURRENT CS \n\n")
        cat("\t 1) pregnancies with a previous CS (SECFORE) summary: \n")        
        print(table(dat$SECFORE,useNA="a"))
        cat("\n pregnancies with a previous CS (SECFORE) split by year: \n")        
        print(table(dat$SECFORE,dat$AR,useNA="a"))
        
        cat("\n \t 2) pregnancies with a current CS (ELEKAKUT) summary: \n")        
        print(table(dat$ELEKAKUT,useNA="a"))
        cat("\t (1 = elective, 2 = acute) \n")        
        
        cat("\n pregnancies with a current CS (ELEKAKUT) split by year: \n")        
        print(table(dat$ELEKAKUT,dat$AR,useNA="a"))
        cat("\n note that before 1999 there is no data! \n")        
        
        cat("\n exclusions will be done if SECFORE=1 or ELEKAKUT=1 \n")        
        
        bad_rows = which((dat$ELEKAKUT==1)|(dat$SECFORE==1))
        
        cat("\n in total ",length(bad_rows)," rows will be removed \n")        
        cat("\n in total ",nrow(dat)," are remaining \n")        
        dat = dat[-bad_rows,]
        dat
        
}

fun_ICDcodes = function(dat) {
        ### based on ICD codes reported in MFR
        
        icd_codes = c("^O40","^O41","^O42","^O43","^O44","^O45","^O46","^762","^641","^642")
        cat("ICD CODES:
            the following ICD codes will be searched in MDIAG and BDIAG columns:
            ",icd_codes,"
            ^O40 - Onormalt stor mängd fostervatten (-1d)
            ^O41	Oligohydramnion and Co  (+4 d)
            ^O42	För tidig hinnbristning, types... (-35 d)
            ^O43	Placentalt transfusionssyndrom, Missbildning, Patologiskt fastsittande, accreta/increta/percreta (-31 d)
            ^O44	Placenta praevia, types... (-23 d)
            ^O45	För tidig avlossning av placenta, types... (-29 d)  (abruptio)
            ^O46	Blödning före förlossningen, types... (-15 d)
            ^641  Blödning i sen graviditet, för tidig avlossning av moderkakan och T föreligg ande moderkaka
            Haemorrhagia in graviditate posteriore, abruptio placentae et placenta praevia 
            ^642  Hypertoni som komplikation till graviditet, förlossning och barnsängs (-9 d)
            ")
        
        tbl = rixs = NULL
        for (icd_code in icd_codes) {
                print(icd_code)
                rix = NULL
                for(j in grep("^MDIAG|^BDIAG",colnames(dat))) {
                        rix = c(rix,grep(icd_code,dat[,j]))
                }
                rix = unique(rix)
                phe = rep(0,nrow(dat)); phe[rix]=1
                coefs = coef(summary(lm(dat$GRDBS ~ phe)))
                df = data.frame(condition=icd_code,N=length(rix),beta=round(coefs[2,1],1),pval=coefs[2,4])
                tbl = rbind(tbl,df)
                rixs = c(rixs,rix)
                rm(df,coefs,phe,rix)
        }
        # preview results
        rixs = unique(rixs)
        
        cat("numbers of rows that are selected for exclusion due to ICD codes: 
            ")
        print(tbl)
        
        cat(" in total ",length(rixs), " will be removed") # selected for exclusion : 67 272
        dat = dat[-rixs,]
        cat(" there are ",nrow(dat), " rows remaining") 
        dat
}

fun_GAmiss = function(dat) {
        cat("MISSING GESTATIONAL AGE:
            ")
        good_rows = which(!is.na(dat$GRDBS))
        cat("in total ",nrow(dat)-length(good_rows),"rows will be removed")
        dat = dat[good_rows,]
        cat("in total ",nrow(dat),"left remaining")
        dat
}

fun_GAdating = function(dat) {
        cat("UNRELIABLE GESTATIONAL AGE DATING METHOD:
            ")
        bad_rows = which(dat$GRMETOD %in% c(0,3,4,9,11,12,13))
        cat("in total ",length(bad_rows),"rows will be removed")
        dat = dat[-bad_rows,]
        cat("in total ",nrow(dat),"left remaining")
        dat
}

fun_MHmiss = function(dat) {
        cat("MISSING (or unreasonable) MATERNAL HEIGHT:
            ")
        cat(" - range of present maternal height: ",paste(range(dat$MLANGD,na.rm=T),collapse="-"),"cm
            ")
        
        # maternal height (set to missing before deleting later)
        # PROBLEM: MatHgh is only present for years 1982-2012, so 1973-1981 are lost!
        # I will preserve the 1973-1981 data for analyses without adjustments
        dat$MLANGD[which(dat$MLANGD=="")]=NA # maternal height must be present for adjustments
        dat$MLANGD = as.numeric(dat$MLANGD)
        dat$MLANGD[which(dat$MLANGD<140)]=NA # suspicious-height threshold determined by Julius (also, no sib-similarity below this thr)
        dat$MLANGD[which(dat$MLANGD>210)]=NA # almost incredible values
        #ix=sample(nrow(dat),1e3,replace=F); plot(dat$MLANGD[ix]~dat$AR[ix]); rm(ix)
        #table(dat$AR[which(!is.na(dat$MLANGD))])
        
        cat(" - in total ",sum(is.na(dat$MLANGD)),"rows will be removed
            ")
        dat = dat[which(!is.na(dat$MLANGD)),]
        cat(" - in total ",nrow(dat),"left remaining
            ")
        cat(" - range of remaining maternal height: ",paste(range(dat$MLANGD,na.rm=T),collapse="-"),"cm
            ")
        
        dat
}

fun_very_nordic = function(dat) {
        
        # countries
        #tbl = table(dat$MFODLAND)
        #df = data.frame(country=names(tbl),cnt=as.numeric(tbl),stringsAsFactors = F)
        #df = df[order(df$cnt,decreasing = T),]
        #df[1:25,]
        #df[26:50,]
        
        nordic = c("SVERIGE","FINLAND","NORGE","DANMARK")
        europe = c(nordic,"POLEN","TYSKLAND","FRANKRIKE","ESTLAND","ISLAND","UKRAINA","SPANIEN","GREKLAND")
        
        cat("NORDIC PARENTS: \n \t (all indicators point that parents are nordic) \n")
        cat("\t nordic countries:",nordic,"\n")
        
        good_rix = which((dat$MFODLAND %in% nordic)&(dat$MNAT %in% nordic)&(dat$FNAT %in% nordic))
        #sub = dat[which((dat$MFODLAND %in% nordic)&(dat$MNAT %in% nordic)),]
        #sub = dat[which(dat$MFODLAND %in% nordic),]
        
        cat("\t in total",length(good_rix),"rows will remain \n")
        cat("\t in total",nrow(dat)-length(good_rix),"rows will be removed \n")
        cat("\t report on concordance: \n")
        print(table(momBRTH=dat$MFODLAND[good_rix],momNAT=dat$MNAT[good_rix]))
        print(table(momNAT=dat$MNAT[good_rix],dadNAT = dat$FNAT[good_rix]))
        
        dat[good_rix,]
} # very stringent

fun_mom_nordic = function(dat) {
        nordic = c("SVERIGE","FINLAND","NORGE","DANMARK")
        #europe = c(nordic,"POLEN","TYSKLAND","FRANKRIKE","ESTLAND","ISLAND","UKRAINA","SPANIEN","GREKLAND")
        
        cat("NORDIC PARENTS: \n \t (MFODLAND indicates nordic origin) \n")
        cat("\t nordic countries:",nordic,"\n")
        
        good_rix = which(dat$MFODLAND %in% nordic)
        #sub = dat[which((dat$MFODLAND %in% nordic)&(dat$MNAT %in% nordic)),]
        #sub = dat[which(dat$MFODLAND %in% nordic),]
        
        cat("\t in total",length(good_rix),"rows will remain \n")
        cat("\t in total",nrow(dat)-length(good_rix),"rows will be removed \n")
        cat("\t report on concordance: \n")
        
        dat[good_rix,]
} # not-so stringent

fun_visualize_exclusions_by_year = function(year_matrix) {
        years = as.numeric(year_matrix$year)
        ref_vals = year_matrix$step1
        m = as.matrix(year_matrix[,-(1:2)]) # all steps without initial numbers (2nd col)
        
        # define a function used to estimate binomial p-value
        p_binom = function(n_this_after,n_this_before,n_all_after,n_all_before) {
                binom.test(n_this_after,n_this_before,p=n_all_after/n_all_before)$p.value
        }
        plot_dat = NULL
        for (j in 1:ncol(m)) {
                m[which(is.na(m[,j])),j] = 0 # TEMPORARY ****
                
                if (j==1) { 
                        # estimate a p-value for each year
                        ps = NULL
                        for(i in 1:nrow(m)) {
                                p = p_binom(m[i,j],ref_vals[i],sum(m[,j]),sum(ref_vals))
                                ps = c(ps, p); rm(p)
                        }         
                        # estimate effect size ("fraction lost") for each year
                        tmp = data.frame(years, step=j, fr = m[,j]/ref_vals,p_binom = ps)
                        plot_dat = rbind(plot_dat,tmp); rm(tmp,ps)
                } else {
                        # estimate a p-value for each year
                        ps = rep(NA,nrow(m))
                        gix = which((m[,j-1]!=0)&(m[,j]!=0)) # good-rows' indexes
                        for(i in gix) { # for each year
                                p = p_binom(m[i,j],m[i,j-1],sum(m[gix,j]),sum(m[gix,j-1]))
                                if (is.logical(p)) p = 1
                                ps[i] = p; rm(p)
                        }
                        # estimate effect size ("fraction lost") for each year
                        tmp = data.frame(years, step=j, fr = m[,j] / m[,j-1],p_binom = ps)
                        plot_dat = rbind(plot_dat,tmp); rm(tmp,ps)
                }
        }
        
        library(RColorBrewer)
        breaks = c(0,0.1,0.3,0.5,0.7,0.8,  0.9,0.91,0.92,0.93,0.94 ,0.95,0.96,0.97,0.98,0.99,1)
        n_cols = length(breaks)-1
        n_blus = floor(n_cols/3)
        n_grns = n_blus
        n_reds = n_cols - n_blus - n_grns
        
        blus = brewer.pal(n_blus+2,"Blues")[2:(n_blus+1)]
        grns = rev(brewer.pal(n_grns+2, "Greens")[2:(n_grns+1)])
        reds = rev(brewer.pal(n_reds+2, "Reds")[2:(n_reds+1)])
        #blus = brewer.pal(n_blus,"Blues")
        #grns = rev(brewer.pal(n_grns, "Greens")
        #reds = rev(brewer.pal(n_reds, "Reds")
        
        colrs = c(reds,grns,blus)
        plot_dat$colrs = as.character(cut(plot_dat$fr,breaks = breaks,labels = colrs)) #brewer.pal(length(breaks)-1, "Spectral")
        plot_dat$colrs[which(is.na(plot_dat$colrs))] = NA
        plot(plot_dat$years,plot_dat$step,pch=19,cex=1.5,xlab="year",ylab="stage",xaxt="n",yaxt="n",
             ylim=c(0,ncol(m)+2),xlim = c(min(years),max(years)+15),col=plot_dat$colrs)
        ix = plot_dat$p_binom < 0.01/(40*10)
        points(plot_dat$years[ix],plot_dat$step[ix],pch=19,cex=0.3,col="white")
        axis(1,at = seq(1975,2010,5),labels = seq(1975,2010,5))
        axis(2,at = seq(ncol(m)),labels = seq(ncol(m)))
        
        legend(x = max(years)+6, y = ncol(m)+1,lwd=6,xjust = 0,title = "sample loss",
               legend = paste("<",round((1-breaks[-length(breaks)])*100,0),"%",sep=""),
               col = colrs,border = F,box.col = 0,y.intersp = 1.2,cex = 0.8)
        abline(v=1973:2012,lty=2,lwd=0.4,col="grey")
        abline(h=seq(ncol(m)),lty=2,lwd=0.4,col="grey")
        # add final summum of changes
        fin_left = m[,ncol(m)]/ref_vals
        fin_left_col = as.character(cut(fin_left,breaks = breaks,labels = colrs))
        fin_left_col[which(is.na(fin_left_col))] = NA
        points(years,rep(ncol(m)+1,length(years)),pch=19,cex=2.5,col=fin_left_col)
        text(rep(2013,ncol(m)),seq(ncol(m)),pos = 4,cex=0.5,
             labels = apply(m,2,function(x) sum(x,na.rm=T)))
        
        data.frame(years=years,loss = round((1-fin_left)*100,1),stringsAsFactors = F)
        
        
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

