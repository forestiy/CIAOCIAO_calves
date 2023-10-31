#This file contians the code for both outcomes AMU and respiratory, changing the modle from plsRlm to plsR and 
the outcome from Respiratory to AMU is needed to run selection for AMU; now it is set for counts of respiratory treatments 
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

