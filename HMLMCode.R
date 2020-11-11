


##################################################################
###  Fits the general heterogeneous logit model described in
# Tutz, G. (2020) Uncertain Choices: the Heterogeneous Multinomial Logit Model,
#                  Sociological Methodology
##################################################################

# The function GHMNL fits a  more general model with two sorts of explanatory variables

##########################        Explanatory variables: distinguish between  

##########################        global variables (not category-specific) like gender,age
##########################        category-specific variables like price in transportation mode
##########################        if all variables are global specify any variable in nameshet, but use Indcats=0
##########################        in Tutz, G. (2020) Uncertain Choices: the Heterogeneous Multinomial Logit Model,
#                                    Sociological Methodology, all variables are glogal

###Basic arguments of the function:

#GHMNL <- function(dat, namesglob,Indglob,namescats, Indcats,nameshet,nameresp, k,penalty){
  
  ### dat:         data in long format (k rows for each observation ) 
  ### namesglob:   names of variables with global effects (subject-specific)
  ### Indglob >0:  global are included, otherwise ignored
  ### namescats:   names of category-specific variables (global parameters) 
  #                if no category-specific variables are in the data set, 
  #                choose any variable and set Indcats =0
  ### Indcats >0:  category-specific variables are included, otherwise ignored
  ### nameshet:    names of variables in heterogeneity term
  #                if nameshet <- c()  no heterogeneity term in the model
  ### nameresp:    name  of response variable
  ### k:           number of response category
  ### penalty:     number >= 0, if penalty=0 pure maximum likelihood
  #                             if penalty>0 ridge penalized maximum likelihood 
  
### Note: data have to be in long format (each observation provides k rows,  the response is indicated by a 1 in 
#                                        the corresponding row, 0 in all other rows)
#

# install relevant packages 
install.packages("dplyr")
install.packages("EffectStars2")
install.packages("matlib")
#install.packages("pracma")


# load relevant packages 
library("dplyr")
library("EffectStars2")
library("matlib")

##### Code for fitting

source("LikelihoodfunctionsOct2020.R")
#############


########################################
###Party choice  data (from package EffectStars2)
########################################

data(election)
summary(election)
head(election)
######


#############################
###### call by names
##############################
select <- dplyr::select #avoid clash with mass

###########################
#transform data into long format and select variables
k <- 5
datx <- cbind( election$Gender,election$West,election$Age)
datx[,1]<- datx[,1]-1  # gender, 1: male 0: female
datx[,2]<- datx[,2]-1  #west, 1:wes1 0:0therwise    

datxblow <- blowuphet(datx,k)
summary(datx)

datrespmult <- as.numeric(election$Partychoice)
blresp <- blowresp(as.matrix(datrespmult),k)

dat <- cbind(blresp,datxblow)
colnames(dat ) <- c("Partychoice", "Gender", "West","Age")
dat<- as.data.frame(dat)


##### dat contains data in long format
head(dat)
summary(dat)
##########################################



#######data fit: specification of variables

nameresp  <- c("Partychoice")           # response     
namesglob <- c("Gender", "West","Age")  # explanatory variables
namescats <- c("West")                  # heterogeneity variable
nameshet  <- c()                          # no heterogeneity

Indglob <- 1  # include predictors
Indcats <- 0  # no category-specific variables included

######################################
###model without heterogeneity

pen <- 0.0# no penalty
Tr <- GHMNL(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr


#######################################################
nameshet  <- c("Age")  #heterogeneity in age

Tr <- GHMNL(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr

##################################
nameshet  <- c("Gender", "West","Age") #heterogeneity in gender, west, age
 
Tr <- GHMNL(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr

######################
nameshet  <- c("West","Age") #heterogeneity in gender, west

Tr <- GHMNL(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr
###




 
################################################
#####contraceptive data
################################################

go <- read.table(cmc.data",sep=",")

summary(go)
colnames(go) <- c("wage", "wedu", "hedu", "child", "wrel","wwork","hocc","livind","medexpo", "contra")

n <- dim(go)[1]
for(i in 1:n){if(go$wedu[i] >1) {go$wedubin[i] <- 1}}
for(i in 1:n){if(go$wedu[i] < 2) {go$wedubin[i] <- 0}}

summary(go)

#1. Wife's age                     (numerical)
#   2. Wife's education               (categorical)      1=low, 2, 3, 4=high
#3. Husband's education            (categorical)      1=low, 2, 3, 4=high
#   4. Number of children ever born   (numerical)
#   5. Wife's religion                (binary)           0=Non-Islam, 1=Islam
#6. Wife's now working?            (binary)           0=Yes, 1=No
#   7. Husband's occupation           (categorical)      1, 2, 3, 4
#8. Standard-of-living index       (categorical)      1=low, 2, 3, 4=high
#9. Media exposure                 (binary)           0=Good, 1=Not good
#10. Contraceptive method used     (class attribute)  1=No-use 2=Long-term 3=Short-term

go$wedu <- as.factor(go$wedu) 
go$wedubin <- as.factor(go$wedubin) 
summary(go)

#########################################################
#transform data into long format and select variables

datx <- cbind( go$wage,go$child,go$wrel,go$wedu, go$hedu )
datxblow <- blowuphet(datx,k)

p <- 5
#pbas1 <- p+1
k <- 3

datrespmult <- as.numeric(go$contra)
blresp <- blowresp(as.matrix(datrespmult),k)

dat <- cbind(blresp,datxblow)
colnames(dat ) <- c("contra","agew", "child", "relw", "eduw", "eduh" )
dat<- as.data.frame(dat)

head(dat)
summary(dat)
###########################################

#####################################################
####### data fit: specification of variables
nameresp  <- c("contra")          
namesglob <- c("agew", "child", "relw", "eduw", "eduh")
namescats <- c("agew")
nameshet  <- c()    #no heterogeneity

Indglob <- 1  # include predictors
Indcats <- 0  # no category-specific variables included


######################################
###model without heterogeneity

pen <- 0.0# no penalty
Tr <- GHMNL(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr

######################################
###model with  heterogeneity
nameshet  <- c("eduw")  #heterogeneity in eduw

Tr <- fitHMLMCatspecN3(dat, namesglob,Indglob,namescats,Indcats,nameshet,nameresp,k,pen)
Tr




