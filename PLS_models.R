#This file contains the code for both outcomes AMU and respiratory, changing the model from plsRlm to plsR and 
#the outcome from Respiratory to AMU is needed to run selection for AMU; now it is set for counts of respiratory treatments 

#Power analysis is at the bottom oft
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#________________Clear_code_for_dependent_PLS____________________
rm(list=ls()) 
library(plsRglm)
library(dplyr)
library(tibble)
library(Metrics)
library(rlist)
library(pls)

#loading data+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("Calves_Clean.Rdata")
load("FDF_sup.Rdata")

#loading functions+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
############
`%notin%` <- Negate(`%in%`)

round.off <- function (x, digits=0) 
{
  posneg = sign(x)
  z = trunc(abs(x) * 10 ^ (digits + 1)) / 10
  z = floor(z * posneg + 0.5) / 10 ^ digits
  return(z)
}

VIP_PLS <- function(object) {
  
  #SS1 <- c(object$Yloadings[1,])^2 * colSums(object$scores^2)
  #SS2 <- c(object$Yloadings[2,])^2 * colSums(object$scores^2)
  #SS<-SS1+SS2
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}



VIP_PLSRGLM <- function(object) {
  
  SS <- c(object$CoeffC)^2 * colSums(object$tt^2)
  Wnorm2 <- colSums(object$wwnorm^2)
  SSW <- sweep(object$wwnorm^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


## VIPjh returns the VIP of variable j with h components
VIPjh_PLSRGLM <- function(object, j, h) {
  
  
  b <- c(object$CoeffC)[1:h]
  T <- object$tt[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$wwnorm[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}


######
#End of functions loading
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#check for missings
collect_min<-c()
for(i in 1:length(FDF)){
  collect_min<-append(collect_min,sum(is.na(FDF[,i])))
}
FDF<-FDF[,-which((collect_min/nrow(FDF))>0)]
#ceiling(which(is.na(FDF))/nrow(FDF))
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

# Function to identify variables with same value >70%

collect<-c()      
for (i in 1:length(FDF)){
  if(any(is.na(FDF[,i]))==F){
    a<-FDF %>% 
      dplyr::group_by(FDF[,i]) %>%
      dplyr::summarise(count = n (),fr=n()*100/nrow(FDF))
    if (any(a[,3]>70)){
      collect<-append(collect,i)
    }
    
  }else {
    a<-table(FDF[-which(is.na(FDF[,i])),i])/sum(table(FDF[-which(is.na(FDF[,i])),i]))
    
    
    if (any(a>0.70)){
      collect<-append(collect,i)
    }
  }
}


# Final orientation
FDF<-FDF[,-c(collect,which(colnames(FDF) %in% c("H3_Begin_date_388_",
                                                "H3_End_date_389_",
                                                "B3ii_herd_size_124_",
                                                "B7_Places_Num_Stk_134_",
                                                "B7_Places_Num_Bbx_132_",
                                                "brok_Quant_C11iii",
                                                "mais_Quant_C11iii",
                                                "straw_Quant_C11iii",
                                                "muesli_Quant_C11iii",
                                                "C4i_Water_access_Bbx_174_",
                                                "C4i_Milk_equip_height_cm_Bbx_161_",
                                                "C4i_Water_equip_height_cm_ZkDr_194_",
                                                "C4i_Milk_equip_height_cm_ZkDr_171_",
                                                "H8i_Sick_calves_House",
                                                "C4i_Milk_equip_height_cm_Stk_166_",
                                                "C15_Vitamins_other_ZkDr_286_",
                                                "AIAO_farm_B3i_BIN",
                                                "AIAO_at_least_stal_B3i_BIN",
                                                "G2Aii_Health_check_time_Bbx_381_",
                                                "F5_Vent_syst_mech_reg_inlaat_BIN",
                                                "AIAO_stal_B3i_BIN",
                                                "Solid_feed_access_Bbx_C11_BIN",
                                                "G3_How_where_any_treatm_record_384_",
                                                "Stragglers_manag_H7_BIN",
                                                "I3ii_Bbx_clean_score_484_",
                                                "Stragglers_manag_H7_NUM",
                                                "Stragglers_manag_H7_BIN",
                                                "J2i_Motivation_566_",
                                                "J2ii_Busin_facility_improv_567_",
                                                "BB_with_other_ages_F5"   ,                  
                                                "BB_n_str_BIN_F5"  ,
                                                #"I5i_Fumigation_493_",
                                                "H8ii_Nonantibiotic_Naam_NSAID",
                                                "H8ii_Nonantibiotic_Naam_Vit",
                                                "I1i_How_med_given_Feed_Mixer",
                                                "A1_Province_4_Cattle_dens")))]


FDF<-FDF[c(1:(which(colnames(FDF)=="A2iii_Financ_Involv_7_")-1),which(colnames(FDF)=="A3ii_Tot_num_kalfplaces_8_"),which(colnames(FDF)=="A2iii_Financ_Involv_7_"),(which(colnames(FDF)=="A2iii_Financ_Involv_7_")+2):length(FDF))]
glimpse(FDF)


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Start of analysis
################
CV_choices<-rep(c(1,2,3,4,6,9),1)#sample(c(1,2,3,4,6,9),10,replace=T)


CONF_list_FIIN<-list()
FINAL_list_mod_ALL<-list()
for (FIN_RUN in CV_choices){
  olp<-ifelse(FIN_RUN==1,1,150)
  
  # First for loop to identify which model has lowest avg rmse by backselecting on vip 20 times 
  final_models_list<-list()
  for (RUNS in 1:olp){
    
    #_Cross_Validation
    ids<-sample(1:nrow(FDF),nrow(FDF),replace=F)
    test_ids<-split(ids, ceiling(seq_along(ids)/FIN_RUN))
    class(test_ids)
    
    test_id<-list()
    train_id<-list()
    for(i in 1:length(test_ids)){
      test_id[[i]]<-as.numeric(unlist(test_ids[i]))
      train_id[[i]]<-sort(ids)[-as.numeric(unlist(test_ids[i]))]
    }
    
    
    
    AVG_RMSE<-c()
    remove<-c()
    RMSE_InitDF<-NULL
    predictors_names_init<-colnames(FDF)[c(which(colnames(FDF)=="A2iii_Financ_Involv_7_"):length(FDF))]
    CONF_list<-list()
    
    for(G in 1:(length(predictors_names_init)-4)){
      
      
      VIPS_folds<-NULL
      RMSE_Init<-NULL
      for(i in 1:length(test_id)){
        
        
        
        #Y<-as.matrix(cbind(FDF$DDDA[train_id[[i]]],FDF$H5i_Respiratory_Tr[train_id[[i]]]))
        
        Y<-as.matrix(cbind(FDF$H5i_Respiratory_Tr[train_id[[i]] ])) #$DDDA
        
        if(is.null(remove)){
          predictors_names_init<-colnames(FDF)[c(which(colnames(FDF)=="A2iii_Financ_Involv_7_"):length(FDF))]
          
          #X<-as.matrix(FDF[train_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names)])
          
          Xtrain<-as.matrix(FDF[train_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names_init)])
          Xtest<-as.matrix(FDF[test_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names_init)])
          
          PLS_model_1<-plsRglm(Y,Xtrain,modele="pls-glm-poisson",scaleX = T,scaleY = T,verbose=F)
          fitted_Init<-predict.plsRglmmodel(PLS_model_1,Xtest,type="response")
          
          #PLS_model_1<-pls::plsr(I(Y) ~ I(Xtrain), center=T,scale=T,ncomp=2)
          #fitted_Init <- data.frame(predict(PLS_model_1,ncomp=2,
           #                      newdata = Xtest,
            #                     type="response"))
          #PLS_model_1 <- plsR(Y,Xtrain,scaleX = T,scaleY = T)
          #fitted_Init <- predict(PLS_model_1, newdata = Xtest,type="response")
          
        }else{ 
          predictors_names<- predictors_names_init[-which(predictors_names_init %in% (remove))]
          #X<-as.matrix(FDF[train_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names)])
          
          Xtrain<-as.matrix(FDF[train_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names)])
          Xtest<-as.matrix(FDF[test_id[[i]],c("A3ii_Tot_num_kalfplaces_8_",predictors_names)])
          
          PLS_model_1<-plsRglm(Y,Xtrain,modele="pls-glm-poisson",scaleX = T,scaleY = T,verbose=F)
          fitted_Init<-predict.plsRglmmodel(PLS_model_1,Xtest,type="response")
          
          #PLS_model_1<-pls::plsr(I(Y) ~ I(Xtrain), center=T,scale=T,ncomp=2)
          
          #fitted_Init <- data.frame(predict(PLS_model_1,ncomp=2,
            #                     newdata = Xtest,
             #                    type="response"))
          #PLS_model_1 <- plsR(Y,Xtrain,scaleX = T,scaleY = T)
          #fitted_Init <-  predict(PLS_model_1, newdata = Xtest,type="response")
          
          
        }
        
        
        
        VIPS_folds<-rbind(VIPS_folds,VIP_PLSRGLM(PLS_model_1)[2,])
        
        
        
        #RMSE_Init<-append(RMSE_Init,rmse(fitted_Init[,1],FDF$DDDA[test_id[[i]]])+
         #             rmse(fitted_Init[,2],FDF$H5i_Respiratory_Tr[test_id[[i]]]) )
        RMSE_Init<-append(RMSE_Init,rmse(fitted_Init,FDF$H5i_Respiratory_Tr[test_id[[i]]]) ) #$DDDA be careful here how the object fitted_Init is coming out of the different models
        
      }
      RMSE_InitDF<-rbind(RMSE_InitDF,RMSE_Init)
      AVG_RMSE<-append(AVG_RMSE,mean(RMSE_Init))
      
      # #///////////////////////////////////////////////////////
      # #Chcek of coef for conf effects
      # 
      # if(is.null(remove)){
      # PLS_model_CONF<-plsRglm(FDF$DDDA,
      #                         FDF[,c("A3ii_Tot_num_kalfplaces_8_",
      #                                predictors_names_init)],modele="pls-glm-gaussian",scaleX = T,scaleY = T)
      # }else{
      #   PLS_model_CONF<-plsRglm(FDF$DDDA,
      #                           FDF[,c("A3ii_Tot_num_kalfplaces_8_",
      #                                  predictors_names_init[-which(predictors_names_init %in% (remove))])],modele="pls-glm-gaussian",scaleX = T,scaleY = T)
      # 
      # }
      # CONF_list[[length(CONF_list)+1]]<-PLS_model_CONF$Coeffs
      # #///////////////////////////////////////////////////////
      
      to_remove<-  names(which.min(colMeans(VIPS_folds)[-1]))
      
      #remove<-append(remove,substr(to_remove, 10, nchar(to_remove)) )

      remove<-append(remove,to_remove )
      
      print(c(paste("/FIN_RUN for ALL_BEST for dif CV choices->",FIN_RUN),
              paste("//// RUNS for selecting BW models->",RUNS),
              paste("//// G for variable number removed", length(FDF),"->",G),
              paste("//// i for CV folds within CV choices->",i)))
      print(remove)
      
      
      
    }
    
    print(plot(AVG_RMSE))
    
    
    final_models_list[[length(final_models_list)+1]]<-c("A3ii_Tot_num_kalfplaces_8_",predictors_names_init[-which(predictors_names_init %in%remove[c(1:(which.min(AVG_RMSE)-1))])])
    
    
    # #///////////////////////////////////////////////////////
    # CONF_list<-data.frame(row.names(list.rbind(CONF_list)),list.rbind(CONF_list))
    # row.names(CONF_list)<-NULL
    # colnames(CONF_list)<-c("VAR","COEF")
    # CONF_list_FIIN[[length(CONF_list_FIIN)+1]]<-CONF_list
    # #///////////////////////////////////////////////////////
  }
  
  FINAL_list_mod_ALL<-append(FINAL_list_mod_ALL,final_models_list)
  
}


#From all selected models checking which is best for unseen data

unlisted_ALL<-sort(table(unlist(FINAL_list_mod_ALL)))

PREDS_list<-list()
for (i in unique(unlisted_ALL)){
  PREDS_list[[length(PREDS_list)+1]]<-names(which(sort(unlisted_ALL)>=i))
}

PREDS_list<-c(FINAL_list_mod_ALL,PREDS_list)

###___The check
top_model<-c()

for (o in c(1,rep(c(2,3,4,6,9),20))){
  
  #cross_validation
  ids<-sample(1:nrow(FDF),nrow(FDF),replace=F)
  test_ids<-split(ids, ceiling(seq_along(ids)/o))
  
  test_id<-list()
  train_id<-list()
  for(w in 1:length(test_ids)){
    test_id[[w]]<-as.numeric(unlist(test_ids[w]))
    train_id[[w]]<-sort(ids)[-as.numeric(unlist(test_ids[w]))]
  }
  
  
  
  rmseRES<-list()
  
  for (g in 1:length(PREDS_list)){
    rmse_C<-c()
    for (i in 1:length(test_id)){
      
      #Y<-as.matrix(cbind(FDF$DDDA[train_id[[i]]],FDF$H5i_Respiratory_Tr[train_id[[i]]]))
      
      Y<-as.matrix(cbind(FDF$H5i_Respiratory_Tr[train_id[[i]] ])) #$DDDA
      
      #PLSRGLM_model_Final<-plsRglm(FDF$DDDA[train_id[[i]]],FDF[train_id[[i]],PREDS_list[[g]]],modele="pls-glm-gaussian",scaleX = T,scaleY = T,verbose=F)
      #rmse_C<-append(rmse_C, rmse(FDF$DDDA[test_id[[i]] ],predict.plsRglmmodel(PLSRGLM_model_Final,FDF[test_id[[i]],PREDS_list[[g]] ],type="response")))
      
      
      Xtrain<-as.matrix(FDF[train_id[[i]],PREDS_list[[g]] ])
      Xtest<-as.matrix(FDF[test_id[[i]],PREDS_list[[g]] ])
      
      PLS_model_FINAL<-plsRglm(Y,Xtrain,modele="pls-glm-poisson",scaleX = T,scaleY = T,verbose=F)
      fitted_Init<-predict.plsRglmmodel(PLS_model_FINAL,Xtest,type="response")
      
      #PLS_model_FINAL<-pls::plsr(I(Y) ~ I(Xtrain), center=T,scale=T,ncomp=2)
      
      #fitted_Init <- data.frame(predict(PLS_model_FINAL,ncomp=2,
       #                                 newdata = Xtest,
        #                                type="response"))
      #PLS_model_FINAL <- plsR(Y,Xtrain,scaleX = T,scaleY = T)
      #fitted_Init <- predict(PLS_model_FINAL, newdata = Xtest,type="response")

      #rmse_C<-append(rmse_C,rmse(fitted_Init[,1],FDF$DDDA[test_id[[i]]])+
       #                   rmse(fitted_Init[,2],FDF$H5i_Respiratory_Tr[test_id[[i]]]) )
      
      
      rmse_C<-append(rmse_C,rmse(fitted_Init,FDF$H5i_Respiratory_Tr[test_id[[i]]]) ) #$DDDA
      
      
      
      
      
    }
    
    print(c(o,g,i))
    
    rmseRES[[length(rmseRES)+1]]<-mean(rmse_C)
    
  }
  
  top_model<-append(top_model,which.min(unlist(rmseRES)))
  
}

PREDS_list<-c(PREDS_list,names(which.max(table(top_model))),table(top_model))








#################################################################################################
#Power analysis using pwr.f2test()
#################################################################################################
#setting up bootstrap to identify significance of all the variables involved 
n_rep<-1
store_df<-NULL
while(n_rep<2001){
  
  ndim<-21
  # with this functino we generate multivariate random variables usign the covariance structure of the observed data
  #of the selected AMU model
  DF_boot<- mvrnorm(n=sum(185 , 64),mu=colMeans(cbind(Outcomes,Predictors)),Sigma=cov(cbind(Outcomes,Predictors))) )

  #This the bootstrap
for ( N_req in c(c(20,30),seq(40,270,10),276) ) {

  DF_boot_n<-DF_boot[sample(c(1:nrow(DF_boot)),N_req,replace=T),]
 
  
  
  
Outcomes_boot<-as.matrix(cbind(DF_boot_n[,1]))
Predictors_boot<-as.matrix(cbind(DF_boot_n[,-1]))

#CV 10 fold to reduce overfitting and compare best case with 36
#cross_validation
ids<-sample(1:nrow(DF_boot_n),nrow(DF_boot_n),replace=F)
test_ids<-split(ids, ceiling(seq_along(ids)/ifelse(N_req<=21,3,
                                                   round(length(ids)/10,0))))

test_id<-list()
train_id<-list()
for(w in 1:length(test_ids)){
  test_id[[w]]<-as.numeric(unlist(test_ids[w]))
  train_id[[w]]<-sort(ids)[-as.numeric(unlist(test_ids[w]))]
}




stor_rmse<-c()
stor_RSR<-c()
stor_coeflm1<-c()
stor_coeflm2<-c()
stor_PVlm1<-c()
stor_PVlm2<-c()
stor_efstdS1<-c()
stor_efstdS2<-c()
stor_power_test<-c()
stor_power_train<-c()
stor_R2_test<-c()
stor_R2_train<-c()
stor_f2_train<-c()
stor_f2_test<-c()
stor_RS1_pvalue<-c()
for(ui in 1:length(test_id)){
  if(length(unlist(test_id[ui]))>=3){
modpls <-NULL

    modpls <- plsR(Outcomes_boot[unlist(train_id[ui]) ],verbose=F,
               Predictors_boot[unlist(train_id[ui]),],scaleX = T, scaleY = T,nt=2)

# Sample real and predicted values
lm_model<-lm(Outcomes_boot[unlist(train_id[ui]) ] ~ modpls$tt[,1] + modpls$tt[,2])


train_LMcoef_1<-coefficients(lm_model)[2]
train_LMcoef_2<-coefficients(lm_model)[3]
train_LMpval_1<-coef(summary(lm_model))[2,4]
train_LMpval_2<-coef(summary(lm_model))[3,4]

effect_sizeD1<-train_LMcoef_1/sd(Outcomes_boot[unlist(train_id[ui]) ])
effect_sizeD2<-train_LMcoef_2/sd(Outcomes_boot[unlist(train_id[ui]) ])


#test
real_values <- DF_boot_n[unlist(test_id[ui]),1]

#train
real_values_train <- DF_boot_n[unlist(train_id[ui]),1]

if(length(unlist(test_id[ui]))==1){
  
  #test for R2 out of sample
  newwd<-data.frame(DF_boot_n[unlist(test_id[ui]),-1])
  colnames(newwd)[1]<-""
  predicted_values <- as.vector(predict.plsRmodel(modpls,newdata= t(newwd) ))
  
  
  #train fit for in sample R2
  neewdtrain<-data.frame(DF_boot_n[unlist(train_id[ui]),-1])
  colnames(neewdtrain)[1]<-""
  predicted_values_train <- as.vector(predict.plsRmodel(modpls,newdata= t(neewdtrain) ))
  
}else{
  #TEST
predicted_values <- as.vector(predict.plsRmodel(modpls,newdata=DF_boot_n[unlist(test_id[ui]),-1]))


#train
predicted_values_train <- as.vector(predict.plsRmodel(modpls,newdata=DF_boot_n[unlist(train_id[ui]),-1]))

}


Rsquared_test<-1 - (mean((real_values-predicted_values)^2) / mean((real_values - mean(real_values_train))^2))


Rsquared_train<-1 - (sum((real_values_train-predicted_values_train)^2) / sum((real_values_train - mean(real_values_train))^2))

power_train<-pwr.f2.test(u = 22, 
                         v = round.off(length(unlist(train_id[ui])) -  as.numeric(modpls$ic.dof[,1][3])),
                         f2 = (Rsquared_train/(1 - Rsquared_train)), #0.27
                         sig.level = 0.05)$power
if(Rsquared_test==-Inf | Rsquared_test<0 ){
  power_test<-0
  
 
}else{
tryCatch({

 # modpls$ic.dof[,1][3]
  
power_test<-pwr.f2.test(u = 22, 
                       v = round.off(length(unlist(train_id[ui])) -  as.numeric(modpls$ic.dof[,1][3])),
                       f2 =(Rsquared_test/(1 - Rsquared_test)), #0.27
                       sig.level = 0.05)$power




},error=function(e){power_test<-NA;break;print("OVERCOMING ERROR!!!!!!!!!")})
}

S1_pvalue<-ifelse(Rsquared_test>=R_S1,1,0)

R2_test<-Rsquared_test
R2_train<-Rsquared_train

f2_test<-Rsquared_test/(1 - Rsquared_test)
f2_train<-Rsquared_train/(1 - Rsquared_train)


stor_rmse<-append(stor_rmse,rmse(real_values,predicted_values))
stor_RSR<-append(stor_RSR,rmse(real_values,predicted_values)/sd(real_values))
stor_coeflm1<-append(stor_coeflm1,train_LMcoef_1)
stor_coeflm2<-append(stor_coeflm2,train_LMcoef_2)
stor_PVlm1<-append(stor_PVlm1,train_LMpval_1)
stor_PVlm2<-append(stor_PVlm2,train_LMpval_2)
stor_efstdS1<-append(stor_efstdS1,effect_sizeD1)
stor_efstdS2<-append(stor_efstdS2,effect_sizeD2)
stor_power_test<-append(stor_power_test,power_test)
stor_power_train<-append(stor_power_train,power_train)
stor_R2_test<-append(stor_R2_test,R2_test)
stor_R2_train<-append(stor_R2_train,R2_train)
stor_f2_test<-append(stor_f2_test,f2_test)
stor_f2_train<-append(stor_f2_train,f2_train)

stor_RS1_pvalue<-append(stor_RS1_pvalue,S1_pvalue)
  }
}


if(is.null(modpls)){
  znfg<-cbind(n_rep,N_req,rep(NA,12))
}else{
znfg<-cbind(n_rep,N_req,
        RMSE=mean(stor_rmse,na.rm=T),
     RSR=mean(stor_RSR,na.rm=T),
     lmC_1= mean(stor_coeflm1),
       lmC2=mean(stor_coeflm2),
pval1=mean(stor_PVlm1,na.rm=T),
pval2=mean(stor_PVlm2,na.rm=T),
stdC1=mean(stor_efstdS1,na.rm=T),
stdC2=mean(stor_efstdS2,na.rm=T),
PowerTest=mean(stor_power_test,na.rm=T),
PowerTrain=mean(stor_power_train,na.rm=T),
R2_test=mean(stor_R2_test,na.rm=T),
R2_train=mean(stor_R2_train,na.rm=T),
f2_test=mean(stor_f2_test,na.rm=T),
f2_train=mean(stor_f2_train,na.rm=T),
RS1_pvalue=mean(stor_RS1_pvalue,na.rm=T))

           #AIC= modpls$AIC[3],Rsq=modpls$R2[2])

}

store_df<-rbind(store_df,znfg)



}
  print(n_rep)
n_rep<-n_rep+1

}
