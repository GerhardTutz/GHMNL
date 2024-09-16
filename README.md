 GHMNL

 
Fits the general heterogeneous logit model described in
Tutz, G. (2020) Uncertain Choices: the Heterogeneous Multinomial Logit Model,
                Sociological Methodology
The R function fits the General Heterogeneous Multinomial Logit Model (GHMNL), which can handle two types of covariates: global (chooser-specific) and category-specific (choice-specific) variables.

Related Publications:
Tutz, G. (2020). "Uncertain Choices: The Heterogeneous Multinomial Logit Model." Sociological Methodology 51(1): 86-111. doi:10.1177/0081175020979689
Mauerer, I. & Tutz, G. (2023). "Heterogeneity in General Multinomial Choice Models." Statistical Methods & Applications 32: 129–148. doi:10.1007/s10260-022-00642-5


 The function GHMNL fits a  more general model with two sorts of explanatory variables

        Explanatory variables: distinguish between  

        global variables (not category-specific) like gender,age
        category-specific variables like price in transportation mode
        if all variables are global specify any variable in nameshet, but use Indcats=0
        in Tutz, G. (2020) Uncertain Choices: the Heterogeneous Multinomial Logit Model,
                           Sociological Methodology, all variables are global
    
Basic arguments of the function:

GHMNL <- function(dat, namesglob,Indglob,namescats, Indcats,nameshet,nameresp, k,penalty)
  
    dat:         data in long format (k rows for each observation ) 
   
    namesglob:   names of variables with global effects (subject-specific)
   
    Indglob >0:  global are included, otherwise ignored
   
    namescats:   names of category-specific variables (global parameters) 
                 if no category-specific variables are in the data set, 
                 choose any variable and set Indcats =0
   
    Indcats >0:  category-specific variables are included, otherwise ignored
   
    nameshet:    names of variables in heterogeneity term
                  if nameshet <- c()  no heterogeneity term in the model
   
    nameresp:    name  of response variable
   
    k:           number of response category
   
    penalty:     number >= 0, if penalty=0 pure maximum likelihood
                             if penalty>0 ridge penalized maximum likelihood 
  
Note: data have to be in long format (each observation provides k rows,  the response is indicated by a 1 in 
                                        the corresponding row, 0 in all other rows)
