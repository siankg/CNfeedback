##### Setup #####

rm(list = ls())

library(ncdf4)
library(terra)
library(latex2exp)
library(rnaturalearth)

dir_processed="" #location of processed files
figures_directory="" #location to save figures

setwd(dir_processed)

##### Load forcing files #####

CO2<-read.csv("Supplementary_Table_UoM_GHGConcentrations-1-1-0_annualmeans_v23March2017.csv")

Ndep<-rast("Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012.nc")*365.25*24*60*60 #kg N m-2 s-1 -> kg N m-2 yr-1
Ndep_total<-ncvar_get(nc_open("Ndep_input4MIPs_surfaceFluxes_C4MIP_MPI-B-1pctNdep-1-0_gn_185001-200012_total.nc"),"drynhx")*365.25*24*60*60/10^9 #kg N s-1 -> Tg N yr-1
latitudes_Ndep<-ncvar_get(nc_open("Ndep_zonsum_timsum.nc"),"lat")
Ndep_lat<-ncvar_get(nc_open("Ndep_zonsum_timsum.nc"),"drynhx") #kg N s-1
Ndep_lat_binned <- aggregate(Ndep_lat,by=list(cut(latitudes_Ndep,seq(-90,90,5))),mean,na.rm=TRUE) #kg N s-1

##### Load files #####

nc_files_netAtmosLandCO2Flux<-list.files(pattern=glob2rx("netAtmosLandCO2Flux_*total.nc"))
for(i in 1:length(nc_files_netAtmosLandCO2Flux)){
  nc_temp<-nc_open(nc_files_netAtmosLandCO2Flux[i])
  df<-ncvar_get(nc_temp,"netAtmosLandCO2Flux")*365.25*24*60*60/10^12
  name<-substr(nc_files_netAtmosLandCO2Flux[i],1,nchar(nc_files_netAtmosLandCO2Flux[i])-9)
  assign(name,df)
}

nc_files_fgco2<-list.files(pattern=glob2rx("fgco2_*total.nc"))
for(i in 1:length(nc_files_fgco2)){
  nc_temp<-nc_open(nc_files_fgco2[i])
  df<-ncvar_get(nc_temp,"fgco2")*365.25*24*60*60/10^12
  name<-substr(nc_files_fgco2[i],1,nchar(nc_files_fgco2[i])-9)
  assign(name,df)
}

nc_files_tas<-list.files(pattern=glob2rx("tas_*mean.nc"))
for(i in 1:length(nc_files_tas)){
  nc_temp<-nc_open(nc_files_tas[i])
  df<-ncvar_get(nc_temp,"tas")-273.15
  name<-substr(nc_files_tas[i],1,nchar(nc_files_tas[i])-8)
  assign(name,df)
}

nc_files_cVeg<-list.files(pattern=glob2rx("cVeg_*total.nc"))
for(i in 1:length(nc_files_cVeg)){
  nc_temp<-nc_open(nc_files_cVeg[i])
  df<-ncvar_get(nc_temp,"cVeg")/10^12
  name<-substr(nc_files_cVeg[i],1,nchar(nc_files_cVeg[i])-9)
  assign(name,df)
}

nc_files_nVeg<-list.files(pattern=glob2rx("nVeg_*total.nc"))
for(i in 1:length(nc_files_nVeg)){
  nc_temp<-nc_open(nc_files_nVeg[i])
  df<-ncvar_get(nc_temp,"nVeg")/10^12
  name<-substr(nc_files_nVeg[i],1,nchar(nc_files_nVeg[i])-9)
  assign(name,df)
}

nc_files_fN2O<-list.files(pattern=glob2rx("fN2O_*total.nc"))
for(i in 1:length(nc_files_fN2O)){
  nc_temp<-nc_open(nc_files_fN2O[i])
  df<-ncvar_get(nc_temp,"fN2O")*365.25*24*60*60/10^9
  name<-substr(nc_files_fN2O[i],1,nchar(nc_files_fN2O[i])-9)
  assign(name,df)
}
`fN2O_UKESM1-0-LL_1pctCO2`<-rep(NA,140)
`fN2O_UKESM1-0-LL_1pctCO2-bgc`<-rep(NA,140)
`fN2O_UKESM1-0-LL_1pctCO2-rad`<-rep(NA,140)
`fN2O_UKESM1-0-LL_1pctCO2Ndep`<-rep(NA,140)
`fN2O_UKESM1-0-LL_1pctCO2Ndep-bgc`<-rep(NA,140)

nc_files_fgn2o<-list.files(pattern=glob2rx("fgn2o_*total.nc"))
for(i in 1:length(nc_files_fgn2o)){
  nc_temp<-nc_open(nc_files_fgn2o[i])
  df<-ncvar_get(nc_temp,"fgn2o")*365.25*24*60*60/10^9
  name<-substr(nc_files_fgn2o[i],1,nchar(nc_files_fgn2o[i])-9)
  assign(name,df)
}
`fgn2o_MIROC-ES2L_1pctCO2-bgc`<-rep(NA,140)
`fgn2o_MIROC-ES2L_1pctCO2-rad`<-rep(NA,140)
`fgn2o_MIROC-ES2L_1pctCO2Ndep-bgc`<-rep(NA,140)
`fgn2o_UKESM1-0-LL_1pctCO2`<-rep(NA,140)
`fgn2o_UKESM1-0-LL_1pctCO2-bgc`<-rep(NA,140)
`fgn2o_UKESM1-0-LL_1pctCO2-rad`<-rep(NA,140)
`fgn2o_UKESM1-0-LL_1pctCO2Ndep`<-rep(NA,140)
`fgn2o_UKESM1-0-LL_1pctCO2Ndep-bgc`<-rep(NA,140)
`fgn2o_MPI-ESM1-2-LR_1pctCO2`<-rep(NA,140)
`fgn2o_MPI-ESM1-2-LR_1pctCO2-bgc`<-rep(NA,140)
`fgn2o_MPI-ESM1-2-LR_1pctCO2-rad`<-rep(NA,140)
`fgn2o_MPI-ESM1-2-LR_1pctCO2Ndep`<-rep(NA,140)
`fgn2o_MPI-ESM1-2-LR_1pctCO2Ndep-bgc`<-rep(NA,140)

nc_files_fBNF<-list.files(pattern=glob2rx("fBNF_*total.nc"))
for(i in 1:length(nc_files_fBNF)){
  nc_temp<-nc_open(nc_files_fBNF[i])
  df<-ncvar_get(nc_temp,"fBNF")*365.25*24*60*60/10^9
  name<-substr(nc_files_fBNF[i],1,nchar(nc_files_fBNF[i])-9)
  assign(name,df)
}

nc_files_fNnetmin<-list.files(pattern=glob2rx("fNnetmin_*total.nc"))
for(i in 1:length(nc_files_fNnetmin)){
  nc_temp<-nc_open(nc_files_fNnetmin[i])
  df<-ncvar_get(nc_temp,"fNnetmin")*365.25*24*60*60/10^9
  name<-substr(nc_files_fNnetmin[i],1,nchar(nc_files_fNnetmin[i])-9)
  assign(name,df)
}

##### Average files for time series #####

# Calculate C:N
for(model in c("MIROC-ES2L","UKESM1-0-LL","MPI-ESM1-2-LR")){
  for(experiment in c("1pctCO2","1pctCO2-bgc","1pctCO2-rad","1pctCO2Ndep","1pctCO2Ndep-bgc")){
    assign(paste("CNveg",model,experiment,sep="_"),
           get(paste("cVeg",model,experiment,sep="_"))/
             get(paste("nVeg",model,experiment,sep="_")))
  }}

# Average files
average_files_function<-function(variable){
  assign(paste(variable,"1pctCO2",sep="_"),
         rbind(
           get(paste(variable,"MIROC-ES2L","1pctCO2",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_bgc",sep="_"),
         rbind(
           get(paste(variable,"MIROC-ES2L","1pctCO2-bgc",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-bgc",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-bgc",sep="_"))),
  envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_rad",sep="_"),
         rbind(
           get(paste(variable,"MIROC-ES2L","1pctCO2-rad",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-rad",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-rad",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2Ndep",sep="_"),
         rbind(
           get(paste(variable,"MIROC-ES2L","1pctCO2Ndep",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2Ndep",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2Ndep",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2Ndep_bgc",sep="_"),
         rbind(
           get(paste(variable,"MIROC-ES2L","1pctCO2Ndep-bgc",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2Ndep-bgc",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2Ndep-bgc",sep="_"))),
         envir = .GlobalEnv)
}
average_files_function_N<-function(variable){
  assign(paste(variable,"1pctCO2_N",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2",sep="_")),
           get(paste(variable,"CESM2","1pctCO2",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_bgc_N",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CESM2","1pctCO2-bgc",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2-bgc",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2-bgc",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-bgc",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-bgc",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_rad_N",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2-rad",sep="_")),
           get(paste(variable,"CESM2","1pctCO2-rad",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2-rad",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2-rad",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-rad",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-rad",sep="_"))),
         envir = .GlobalEnv)
}
average_files_function_C<-function(variable){
  assign(paste(variable,"1pctCO2_C",sep="_"),
         rbind(
           get(paste(variable,"BCC-CSM2-MR","1pctCO2",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_bgc_C",sep="_"),
         rbind(
           get(paste(variable,"BCC-CSM2-MR","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2-bgc",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2-bgc",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2-bgc",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_rad_C",sep="_"),
         rbind(
           get(paste(variable,"BCC-CSM2-MR","1pctCO2-rad",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2-rad",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2-rad",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2-rad",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2-rad",sep="_"))),
         envir = .GlobalEnv)
}
average_files_function_all<-function(variable){
  assign(paste(variable,"1pctCO2_all",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2",sep="_")),
           get(paste(variable,"BCC-CSM2-MR","1pctCO2",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2",sep="_")),
           get(paste(variable,"CESM2","1pctCO2",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_bgc_all",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2-bgc",sep="_")),
           get(paste(variable,"BCC-CSM2-MR","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CESM2","1pctCO2-bgc",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2-bgc",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2-bgc",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2-bgc",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-bgc",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2-bgc",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2-bgc",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-bgc",sep="_"))),
         envir = .GlobalEnv)
  assign(paste(variable,"1pctCO2_rad_all",sep="_"),
         rbind(
           get(paste(variable,"ACCESS-ESM1-5","1pctCO2-rad",sep="_")),
           get(paste(variable,"BCC-CSM2-MR","1pctCO2-rad",sep="_")),
           get(paste(variable,"CanESM5","1pctCO2-rad",sep="_")),
           get(paste(variable,"CESM2","1pctCO2-rad",sep="_")),
           get(paste(variable,"CNRM-ESM2-1","1pctCO2-rad",sep="_")),
           get(paste(variable,"IPSL-CM6A-LR","1pctCO2-rad",sep="_")),
           get(paste(variable,"MIROC-ES2L","1pctCO2-rad",sep="_")),
           get(paste(variable,"MPI-ESM1-2-LR","1pctCO2-rad",sep="_")),
           get(paste(variable,"GFDL-ESM4","1pctCO2-rad",sep="_")),
           get(paste(variable,"NorESM2-LM","1pctCO2-rad",sep="_")),
           get(paste(variable,"UKESM1-0-LL","1pctCO2-rad",sep="_"))),
         envir = .GlobalEnv)
}

average_files_function("netAtmosLandCO2Flux")
average_files_function("fgco2")
average_files_function("tas")
average_files_function("cVeg")
average_files_function("nVeg")
average_files_function("CNveg")
average_files_function("fBNF")
average_files_function("fNnetmin")
average_files_function("fN2O")
average_files_function("fgn2o")

##### Calculate feedback #####

calculate_gamma<-function(model,domain){
  v1<-get(paste(domain,model,"1pctCO2",sep="_"))
  v2<-get(paste(domain,model,"1pctCO2-bgc",sep="_"))
  v3<-get(paste("tas",model,"1pctCO2",sep="_"))
  lm <- lm(v3~poly(1:140,4))
  v3_smoothed <- predict(lm,data.frame(x=1:140),interval='confidence',level=0.99)
  
  feedback<-rep(NA,140)
  for(i in 2:140){
    feedback[i]<-(sum(v1[1:i])-sum(v2[1:i]))/(sum(diff(v3_smoothed[1:i,1])))
  }

  return(list(feedback,(sum(v1[1:140])-sum(v2[1:140]))))
  
}
calculate_beta<-function(model,domain){
  v<-get(paste(domain,model,"1pctCO2-bgc",sep="_"))
  feedback<-rep(NA,140)
  for(i in 1:140){
    feedback[i]<-(sum(v[1:i]))/(sum(diff(CO2[1:i,2])))
  }
  
  return(feedback)
  
}
calculate_epsilon<-function(model,domain,bgc=FALSE){
  if(model%in%c("MIROC-ES2L",
                "MPI-ESM1-2-LR",
                "UKESM1-0-LL")){
  if(bgc==FALSE){
    v1<-get(paste(domain,model,"1pctCO2Ndep",sep="_"))
    v2<-get(paste(domain,model,"1pctCO2",sep="_"))
  }
  if(bgc==TRUE){
    v1<-get(paste(domain,model,"1pctCO2Ndep-bgc",sep="_"))
    v2<-get(paste(domain,model,"1pctCO2-bgc",sep="_"))
  }
  feedback<-rep(NA,140)
  for(i in 1:140){
    feedback[i]<-(sum(v1[1:i])-sum(v2[1:i]))/(sum(Ndep_total[1:i]))
  }}
  else{feedback<-rep(NA,140)}
  
  return(feedback)
  
}

calculate_feedbacks<-function(model,modelname,epsilon=FALSE){

  assign(paste("gamma",modelname,"land",sep="_"),calculate_gamma(model,"netAtmosLandCO2Flux")[[1]], envir = .GlobalEnv)
  assign(paste("beta",modelname,"land",sep="_"),calculate_beta(model,"netAtmosLandCO2Flux"), envir = .GlobalEnv)
  assign(paste("gamma",modelname,"ocean",sep="_"),calculate_gamma(model,"fgco2")[[1]], envir = .GlobalEnv)
  assign(paste("beta",modelname,"ocean",sep="_"),calculate_beta(model,"fgco2"), envir = .GlobalEnv)
  assign(paste("gamma",modelname,"land_PgC",sep="_"),calculate_gamma(model,"netAtmosLandCO2Flux")[[2]], envir = .GlobalEnv)

  if(epsilon){
    assign(paste("epsilon",modelname,"land",sep="_"),calculate_epsilon(model,"netAtmosLandCO2Flux"), envir = .GlobalEnv)
    assign(paste("epsilon",modelname,"ocean",sep="_"),calculate_epsilon(model,"fgco2"), envir = .GlobalEnv)
  }

}

calculate_feedbacks("ACCESS-ESM1-5","ACCESSESM15")
calculate_feedbacks("BCC-CSM2-MR","BCCCSM2MR")
calculate_feedbacks("CanESM5","CanESM5")
calculate_feedbacks("CESM2","CESM2")
calculate_feedbacks("CNRM-ESM2-1","CNRMESM21")
calculate_feedbacks("IPSL-CM6A-LR","IPSLCM6ALR")
calculate_feedbacks("MIROC-ES2L","MIROCES2L",epsilon=TRUE)
calculate_feedbacks("MPI-ESM1-2-LR","MPIESM12LR",epsilon=TRUE)
calculate_feedbacks("GFDL-ESM4","GFDLESM4")
calculate_feedbacks("NorESM2-LM","NorESM2LM")
calculate_feedbacks("UKESM1-0-LL","UKESM10LL",epsilon=TRUE)

epsilon_MIROCES2L_land_bgc<-calculate_epsilon("MIROC-ES2L","netAtmosLandCO2Flux",bgc=TRUE)
epsilon_UKESM10LL_land_bgc<-calculate_epsilon("UKESM1-0-LL","netAtmosLandCO2Flux",bgc=TRUE)
epsilon_MPIESM12LR_land_bgc<-calculate_epsilon("MPI-ESM1-2-LR","netAtmosLandCO2Flux",bgc=TRUE)
epsilon_MIROCES2L_ocean_bgc<-calculate_epsilon("MIROC-ES2L","fgco2",bgc=TRUE)
epsilon_UKESM10LL_ocean_bgc<-calculate_epsilon("UKESM1-0-LL","fgco2",bgc=TRUE)
epsilon_MPIESM12LR_ocean_bgc<-calculate_epsilon("MPI-ESM1-2-LR","fgco2",bgc=TRUE)

epsilon_N2O_MIROCES2L_land<-calculate_epsilon("MIROC-ES2L","fN2O")
epsilon_N2O_MPIESM12LR_land<-calculate_epsilon("MPI-ESM1-2-LR","fN2O")
epsilon_N2O_MIROCES2L_ocean<-calculate_epsilon("MIROC-ES2L","fgn2o")

df_feedbacks<-matrix(NA,nrow=11,ncol=12)
i=1
for(model in c("ACCESS-ESM1-5",
               "BCC-CSM2-MR",
               "CanESM5",
               "CESM2",
               "CNRM-ESM2-1",
               "IPSL-CM6A-LR",
               "MIROC-ES2L",
               "MPI-ESM1-2-LR",
               "GFDL-ESM4",
               "NorESM2-LM",
               "UKESM1-0-LL")){
  print(model)
  df_feedbacks[i,]<-c(calculate_gamma(model,"netAtmosLandCO2Flux")[[1]][140],
               calculate_beta(model,"netAtmosLandCO2Flux")[140],
               calculate_epsilon(model,"netAtmosLandCO2Flux")[140],
               calculate_gamma(model,"fgco2")[[1]][140],
               calculate_beta(model,"fgco2")[140],
               calculate_epsilon(model,"fgco2")[140],
               calculate_gamma(model,"netAtmosLandCO2Flux")[[2]],
               calculate_beta(model,"netAtmosLandCO2Flux")[140]*(sum(diff(CO2[1:140,2]))),
               calculate_epsilon(model,"netAtmosLandCO2Flux")[140]*(sum(Ndep_total[1:140])),
               calculate_gamma(model,"fgco2")[[2]],
               calculate_beta(model,"fgco2")[140]*(sum(diff(CO2[1:140,2]))),
               calculate_epsilon(model,"fgco2")[140]*(sum(Ndep_total[1:140])))
  i<-i+1
}

models<-as.data.frame(t(matrix(
  c("ACCESS-ESM1.5",TRUE,
    "BCC-CSM2-MR",FALSE,
    "CanESM5",FALSE,
    "CESM2",TRUE,
    "CNRM-ESM2-1",FALSE,
    "IPSL-CM6A-LR",FALSE,
    "MIROC-ES2L",TRUE,
    "MPI-ESM1.2-LR",TRUE,
    "NOAA-GFDL-ESM4",FALSE,
    "NorESM2-LM",TRUE,
    "UKESM1-0-LL",TRUE),
  ncol=11,nrow=2)))
colnames(models)=c("name","Ncycle")
models$Ncycle<-as.logical(models$Ncycle)

df_feedbacks<-data.frame(models,df_feedbacks)
colnames(df_feedbacks)<-c("name","Ncycle",
                   "gamma_L","beta_L","epsilon_L",
                   "gamma_O","beta_O","epsilon_O",
                   "gamma_L_PgC","beta_L_PgC","epsilon_L_PgC",
                   "gamma_O_PgC","beta_O_PgC","epsilon_O_PgC")
df_feedbacks$gamma_O[2]<-NA
df_feedbacks$beta_O[2]<-NA
df_feedbacks$gamma_O_PgC[2]<-NA
df_feedbacks$beta_O_PgC[2]<-NA

t.test(df_feedbacks$gamma_L[df_feedbacks$Ncycle==TRUE],df_feedbacks$gamma_L[df_feedbacks$Ncycle==FALSE])
t.test(df_feedbacks$beta_L[df_feedbacks$Ncycle==TRUE],df_feedbacks$beta_L[df_feedbacks$Ncycle==FALSE])

t.test(df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==TRUE],df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==FALSE])
t.test(df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==TRUE],df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==FALSE])

##### Calculate feedback by latitude #####

calculate_lat<-function(model,modelname,epsilon=FALSE){
  
  lat_land<-ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"lat")
  lat_ocean<-ncvar_get(nc_open(paste("fgco2",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"lat")
  
  temp_lat<-ncvar_get(nc_open(paste("tas",model,"1pctCO2_zonmean_timdiff.nc",sep="_")),"tas")
  assign(paste("temp",modelname,"lat",sep="_"),
         aggregate(temp_lat,by=list(cut(lat_land,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
         envir = .GlobalEnv)
  
  
  gamma_land_lat<-(ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"netAtmosLandCO2Flux")-
                ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2-bgc_zonsum_timsum.nc",sep="_")),"netAtmosLandCO2Flux"))*60*60*24*365.25/10^12
    #/ncvar_get(nc_open(paste("tas",model,"1pctCO2_zonmean_timdiff.nc",sep="_")),"tas")
  assign(paste("gamma_land",modelname,"lat",sep="_"),
         aggregate(gamma_land_lat,by=list(cut(lat_land,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
         envir = .GlobalEnv)
  
  gamma_ocean_lat<-(ncvar_get(nc_open(paste("fgco2",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"fgco2")-
                     ncvar_get(nc_open(paste("fgco2",model,"1pctCO2-bgc_zonsum_timsum.nc",sep="_")),"fgco2"))*60*60*24*365.25/10^12
  assign(paste("gamma_ocean",modelname,"lat",sep="_"),
         aggregate(gamma_ocean_lat,by=list(cut(lat_ocean,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
         envir = .GlobalEnv)
  
  beta_land_lat<-(ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2-bgc_zonsum_timsum.nc",sep="_")),"netAtmosLandCO2Flux"))*60*60*24*365.25/10^12
    #/(sum(diff(CO2[1:140,2])))
  assign(paste("beta_land",modelname,"lat",sep="_"),
         aggregate(beta_land_lat,by=list(cut(lat_land,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
         envir = .GlobalEnv)
  
  beta_ocean_lat<-(ncvar_get(nc_open(paste("fgco2",model,"1pctCO2-bgc_zonsum_timsum.nc",sep="_")),"fgco2"))*60*60*24*365.25/10^12
  assign(paste("beta_ocean",modelname,"lat",sep="_"),
         aggregate(beta_ocean_lat,by=list(cut(lat_ocean,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
         envir = .GlobalEnv)
  
  if(epsilon){
    epsilon_land_lat<-((ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2Ndep_zonsum_timsum.nc",sep="_")),"netAtmosLandCO2Flux")-
                ncvar_get(nc_open(paste("netAtmosLandCO2Flux",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"netAtmosLandCO2Flux")))*60*60*24*365.25/10^12
    assign(paste("epsilon_land",modelname,"lat",sep="_"),
           aggregate(epsilon_land_lat,by=list(cut(lat_land,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
           envir = .GlobalEnv)
    epsilon_ocean_lat<-((ncvar_get(nc_open(paste("fgco2",model,"1pctCO2Ndep_zonsum_timsum.nc",sep="_")),"fgco2")-
                      ncvar_get(nc_open(paste("fgco2",model,"1pctCO2_zonsum_timsum.nc",sep="_")),"fgco2")))*60*60*24*365.25/10^12
    assign(paste("epsilon_ocean",modelname,"lat",sep="_"),
           aggregate(epsilon_ocean_lat,by=list(cut(lat_ocean,seq(-90,90,5))),mean,na.rm=TRUE)[,2],
           envir = .GlobalEnv)
  }
  
}

calculate_lat("BCC-CSM2-MR","BCCCSM2MR")
calculate_lat("ACCESS-ESM1-5","ACCESSESM15")
calculate_lat("CanESM5","CanESM5")
calculate_lat("CESM2","CESM2")
calculate_lat("CNRM-ESM2-1","CNRMESM21")
calculate_lat("IPSL-CM6A-LR","IPSLCM6ALR")
calculate_lat("MIROC-ES2L","MIROCES2L",epsilon=TRUE)
calculate_lat("MPI-ESM1-2-LR","MPIESM12LR",epsilon=TRUE)
calculate_lat("GFDL-ESM4","GFDLESM4")
calculate_lat("NorESM2-LM","NorESM2LM")
calculate_lat("UKESM1-0-LL","UKESM10LL",epsilon=TRUE)

temp_lat_all_df<-matrix(NA,nrow=11,ncol=36)
gamma_land_lat_all_df<-matrix(NA,nrow=11,ncol=36)
beta_land_lat_all_df<-matrix(NA,nrow=11,ncol=36)
epsilon_land_lat_all_df<-matrix(NA,nrow=11,ncol=36)
gamma_ocean_lat_all_df<-matrix(NA,nrow=11,ncol=36)
beta_ocean_lat_all_df<-matrix(NA,nrow=11,ncol=36)
epsilon_ocean_lat_all_df<-matrix(NA,nrow=11,ncol=36)
i=1
for(model in c("ACCESSESM15",
               "BCCCSM2MR",
               "CanESM5",
               "CESM2",
               "CNRMESM21",
               "IPSLCM6ALR",
               "MIROCES2L",
               "MPIESM12LR",
               "GFDLESM4",
               "NorESM2LM",
               "UKESM10LL")){
  print(model)
  temp_lat_all_df[i,]<-get(paste("temp",model,"lat",sep="_"))
  gamma_land_lat_all_df[i,]<-get(paste("gamma_land",model,"lat",sep="_"))
  beta_land_lat_all_df[i,]<-get(paste("beta_land",model,"lat",sep="_"))
  if(exists(paste("epsilon_land",model,"lat",sep="_"))){
    epsilon_land_lat_all_df[i,]<-get(paste("epsilon_land",model,"lat",sep="_"))
  }
  gamma_ocean_lat_all_df[i,]<-get(paste("gamma_ocean",model,"lat",sep="_"))
  beta_ocean_lat_all_df[i,]<-get(paste("beta_ocean",model,"lat",sep="_"))
  if(exists(paste("epsilon_ocean",model,"lat",sep="_"))){
    epsilon_ocean_lat_all_df[i,]<-get(paste("epsilon_ocean",model,"lat",sep="_"))
  }
  
  i<-i+1
}

gamma_ocean_lat_all_df[2,]<-NA
beta_ocean_lat_all_df[2,]<-NA

##### Analysis #####

beta_land_vec<-c(beta_MIROCES2L_land[140],beta_UKESM10LL_land[140],beta_MPIESM12LR_land[140])
beta_ocean_vec<-c(beta_MIROCES2L_ocean[140],beta_UKESM10LL_ocean[140],beta_MPIESM12LR_ocean[140])
gamma_land_vec<-c(gamma_MIROCES2L_land[140],gamma_UKESM10LL_land[140],gamma_MPIESM12LR_land[140])
gamma_ocean_vec<-c(gamma_MIROCES2L_ocean[140],gamma_UKESM10LL_ocean[140],gamma_MPIESM12LR_ocean[140])
epsilon_land_vec<-c(epsilon_MIROCES2L_land[140],epsilon_UKESM10LL_land[140],epsilon_MPIESM12LR_land[140])
epsilon_ocean_vec<-c(epsilon_MIROCES2L_ocean[140],epsilon_UKESM10LL_ocean[140],epsilon_MPIESM12LR_ocean[140])

setwd(figures_directory)
##### Figure 2 #####

pdf("Fig2.pdf",width=10.5,height=7)

par(mfrow=c(1,1),mar=c(5.1,5.1,4.1,2.1))
layout(matrix(c(1, 2, 3, 4,4,4,4,4,4), nrow = 3, ncol = 3))

plot(0,0,
     col="white",
     yaxt="n",ylab="",
     cex.lab=1.5,#cex.axis=1.5,
     xlab=TeX("$beta$ (Pg C $ppm^{-1}$)"),
     xlim=c(-2,2),
     ylim=c(0,2))

v<-df_feedbacks$beta_L[df_feedbacks$Ncycle==FALSE]
rect(mean(v)-sd(v), 1.25, mean(v)+sd(v), 1.75, col=rgb(224/256,224/256,224/256,0.5),border="black",lwd=2)
segments(mean(v),1.25,mean(v),1.75,col="black",lwd=2)
points(df_feedbacks$beta_L[df_feedbacks$Ncycle==FALSE],jitter(rep(1.5,5),factor=5),
       col="black",pch=20)

v<-df_feedbacks$beta_L[df_feedbacks$Ncycle==TRUE]
rect(mean(v)-sd(v), 1.15, mean(v)+sd(v), 1.85, col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(mean(v),1.15,mean(v),1.85,col="darkorange",lwd=2)
points(df_feedbacks$beta_L[df_feedbacks$Ncycle==TRUE],jitter(rep(1.5,6),factor=5),col="darkorange",pch=20)

v<-na.omit(df_feedbacks$beta_O)
rect(mean(v)-sd(v), 0.25, mean(v)+sd(v), 0.75, col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(mean(v),0.25,mean(v),0.75,col="deepskyblue",lwd=2)
points(df_feedbacks$beta_O,jitter(rep(0.5,11),factor=5),col="deepskyblue",pch=20)

abline(v=0)
title("a) Carbon-concentration feedback",adj=0,cex.main=1.5)
axis(side=2,at=c(0.5,1.5),labels=c("ocean","land"),las=2,cex.axis=1.5)

plot(0,0,col="white",
     yaxt="n",
     xlab=TeX("$gamma$ (Pg C $\\degree C^{-1}$)"),
     ylab="",cex.lab=1.5,#cex.axis=1.5,
     xlim=c(-200,200),#xlim=c(-200,50),
     ylim=c(0,2))

v<-df_feedbacks$gamma_L[df_feedbacks$Ncycle==FALSE]
rect(mean(v)-sd(v), 1.25, mean(v)+sd(v), 1.75, col=rgb(224/256,224/256,224/256,0.5),border="black",lwd=2)
segments(mean(v),1.25,mean(v),1.75,col="black",lwd=2)
points(df_feedbacks$gamma_L[df_feedbacks$Ncycle==FALSE],jitter(rep(1.5,5),factor=5),
       col="black",pch=20)

v<-df_feedbacks$gamma_L[df_feedbacks$Ncycle==TRUE]
rect(mean(v)-sd(v), 1.15, mean(v)+sd(v), 1.85, col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(mean(v),1.15,mean(v),1.85,col="darkorange",lwd=2)
points(df_feedbacks$gamma_L[df_feedbacks$Ncycle==TRUE],jitter(rep(1.5,6),factor=5),col="darkorange",pch=20)

v<-na.omit(df_feedbacks$gamma_O)
rect(mean(v)-sd(v), 0.25, mean(v)+sd(v), 0.75, col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(mean(v),0.25,mean(v),0.75,col="deepskyblue",lwd=2)
points(df_feedbacks$gamma_O,jitter(rep(0.5,11),factor=5),col="deepskyblue",pch=20)

abline(v=0)
title("b) Carbon-climate feedback",adj=0,cex.main=1.5)
axis(side=2,at=c(0.5,1.5),labels=c("ocean","land"),las=2,cex.axis=1.5)

plot(0,0,col="white",
     yaxt="n",ylab="",cex.lab=1.5,#cex.axis=1.5,
     xlab=TeX("$epsilon$ (Pg C Tg $N^{-1}$)"),
     xlim=c(-0.01,0.01),
     ylim=c(0,2))

v<-c(epsilon_MIROCES2L_land[140],epsilon_UKESM10LL_land[140],epsilon_MPIESM12LR_land[140])
rect(mean(v)-sd(v), 1.25, mean(v)+sd(v), 1.75, col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(mean(v),1.25,mean(v),1.75,col="darkorange",lwd=2)
points(v,jitter(rep(1.5,3),factor=5),pch=20,col="darkorange")

v<-c(epsilon_MIROCES2L_ocean[140],epsilon_UKESM10LL_ocean[140],epsilon_MPIESM12LR_ocean[140])
rect(mean(v)-sd(v), 0.25, mean(v)+sd(v), 0.75, col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(mean(v),0.25,mean(v),0.75,col="deepskyblue",lwd=2)
points(v,jitter(rep(0.5,3),factor=5),pch=20,col="deepskyblue")

abline(v=0)
title("c) Carbon-nitrogen feedback",adj=0,cex.main=1.5)
axis(side=2,at=c(0.5,1.5),labels=c("ocean","land"),las=2,cex.axis=1.5)

plot(0,0,col="white",
     axes=FALSE,xlab="",
     ylab=TeX("Cumulative net $CO_2$ uptake (Pg C)"),ylim=c(-1000,2000),cex.lab=1.5,cex.axis=1.5,
     xlim=c(0,9.5))
v<-df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==FALSE]
rect(0,mean(v)-sd(v), 1, mean(v)+sd(v), col=rgb(224/256,224/256,224/256,0.1), border="black",lwd=2)
segments(0,mean(v),1,mean(v),col="black",lwd=2)
v<-df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==TRUE]
rect(1,mean(v)-sd(v), 2, mean(v)+sd(v), col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(1,mean(v), 2, mean(v),col="darkorange",lwd=2)
v<-df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==FALSE]
rect(2.25,mean(v)-sd(v), 3.25, mean(v)+sd(v), col=rgb(224/256,224/256,224/256,0.1), border="black",lwd=2)
segments(2.25,mean(v), 3.25, mean(v),col="black",lwd=2)
v<-df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==TRUE]
rect(3.25,mean(v)-sd(v), 4.25, mean(v)+sd(v), col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(3.25,mean(v), 4.25, mean(v),col="darkorange",lwd=2)
v<-na.omit(df_feedbacks$epsilon_L_PgC)
rect(4.5,mean(v)-sd(v), 5.5, mean(v)+sd(v), col=rgb(255/256,140/256,0,0.1), border="darkorange",lwd=2)
segments(4.5,mean(v), 5.5, mean(v),col="darkorange",lwd=2)
v<-na.omit(df_feedbacks$beta_O_PgC)
rect(6,mean(v)-sd(v), 7, mean(v)+sd(v), col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(6,mean(v), 7, mean(v),col="deepskyblue",lwd=2)
v<-na.omit(df_feedbacks$gamma_O_PgC)
rect(7.25,mean(v)-sd(v), 8.25, mean(v)+sd(v), col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(7.25,mean(v), 8.25, mean(v),col="deepskyblue",lwd=2)
v<-na.omit(df_feedbacks$epsilon_O_PgC)
rect(8.5,mean(v)-sd(v), 9.5, mean(v)+sd(v), col=rgb(0,191/256,255/256,0.1), border="deepskyblue",lwd=2)
segments(8.5,mean(v), 9.5, mean(v),col="deepskyblue",lwd=2)

axis(side=2,at=seq(-1000,2000,500),cex.axis=1.5)
abline(v=5.75,lty=2)
abline(h=0)
points(rep(0.5,5),jitter(df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==FALSE]),pch=20)
points(rep(1.5,6),jitter(df_feedbacks$beta_L_PgC[df_feedbacks$Ncycle==TRUE]),pch=20,col="darkorange")
points(rep(2.75,5),jitter(df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==FALSE]),pch=20)
points(rep(3.75,6),jitter(df_feedbacks$gamma_L_PgC[df_feedbacks$Ncycle==TRUE]),pch=20,col="darkorange")
points(rep(5,3),jitter(na.omit(df_feedbacks$epsilon_L_PgC)),pch=20,col="darkorange")
points(rep(6.5,10),jitter(na.omit(df_feedbacks$beta_O_PgC)),pch=20,col="deepskyblue")
points(rep(7.75,10),jitter(na.omit(df_feedbacks$gamma_O_PgC)),pch=20,col="deepskyblue")
points(rep(9,3),jitter(na.omit(df_feedbacks$epsilon_O_PgC)),pch=20,col="deepskyblue")
text(c(1,3.25,5,6.5,7.75,9),-1000,
     c("carbon-\nconcentration\nfeedback",
       "carbon-\nclimate\nfeedback",
       "carbon-\nnitrogen \nfeedback",
       "carbon-\nconcentration\nfeedback",
       "carbon-\nclimate\nfeedback",
       "carbon-\nnitrogen\nfeedback"),
     xpd=NA,cex=1.25) #srt=30,adj=1,
text(c(1,3.25,5,6.5,7.75,9),-1200,
     c(TeX("($beta_{L}$)"),
       TeX("($gamma_{L}$)"),
       TeX("($epsilon_{L}$)"),
       TeX("($beta_{O}$)"),
       TeX("($gamma_{O}$)"),
       TeX("($epsilon_{O}$)")),
     xpd=NA,cex=1.25)
text(c(3,7.75),-1400,c("land","ocean"),xpd=NA,cex=1.5)

legend("topleft",
       c("land C cycle (5 CMIP6 ESMs)","land C-N cycle (6 CMIP6 ESMs)","ocean (11 CMIP6 ESMs)"),
       fill=c("black","darkorange","deepskyblue"),border=NA,bty="n", cex=1.5)

title("d)",adj=0,cex.main=1.5)

dev.off()

##### Figure 3 #####

pdf("Fig3.pdf",width=10.5,height=10.5)

par(mfrow=c(2,3),mar=c(6.1,5.1,4.1,2.1))

plot(0,0,col="white",
     xlim=c(0,1000),xlab=TeX("$Delta CO_2 \\ (ppm)$"),
     ylim=c(-90,90),ylab="Latitude",cex.axis=1.5,cex.lab=1.5)
abline(v=(sum(diff(CO2[1:140,2]))),lwd=2)
abline(v=0)
title(TeX("\\textbf{a) $Delta CO_2$}"),adj=0,cex.main=1.5)

plot(0,0,
     col="white",
     xlim=c(0,25),xlab=TeX("$Delta T \\ (\\degree C)$"),
     ylim=c(-90,90),ylab="",cex.axis=1.5,cex.lab=1.5)
polygon(c(rev(apply(temp_lat_all_df,2,quantile,na.rm=TRUE,0.975)),
          apply(temp_lat_all_df,2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
lines(colMeans(temp_lat_all_df,na.rm=TRUE),
      seq(-87.5,87.5,5),
      typ="l",lwd=2)
abline(v=0)
title(TeX("\\textbf{b) $Delta T$}"),adj=0,cex.main=1.5)

plot(Ndep_lat*60*60*24*365.25/10^9,latitudes_Ndep,
     typ="l",lwd=2,
     xlim=c(0,300),xlab=TeX("$Delta N_{dep}$ (Tg N)"),
     ylim=c(-90,90),ylab="",cex.axis=1.5,cex.lab=1.5)
abline(v=0)
title(TeX("\\textbf{c) $Delta N_{dep}$}"),adj=0,cex.main=1.5)

plot(0,0,col="white",
     xlim=c(-10,80),xlab=TeX("Cumulative net $CO_2$ uptake (Pg C)"),
     ylim=c(-90,90),ylab="Latitude",cex.axis=1.5,cex.lab=1.5)
polygon(c(rev(apply(beta_land_lat_all_df[!models$Ncycle,],2,quantile,na.rm=TRUE,0.975)),
          apply(beta_land_lat_all_df[!models$Ncycle,],2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(apply(beta_land_lat_all_df[models$Ncycle,],2,quantile,na.rm=TRUE,0.975)),
          apply(beta_land_lat_all_df[models$Ncycle,],2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(255/256,140/256,0,0.1), border = NA)
polygon(c(rev(apply(beta_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.975)),
          apply(beta_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(0,191/256,255/256,0.1), border = NA)
lines(colMeans(beta_land_lat_all_df[!models$Ncycle,],na.rm=TRUE),seq(-87.5,87.5,5),lwd=2)
lines(colMeans(beta_land_lat_all_df[models$Ncycle,],na.rm=TRUE),seq(-87.5,87.5,5),lwd=2,col="darkorange")
lines(colMeans(beta_ocean_lat_all_df,na.rm=TRUE),seq(-87.5,87.5,5),lwd=2,col="deepskyblue")
abline(v=0)
title("d) Carbon-concentration feedback",adj=0,cex.main=1.5)
legend("topright",c("land (C cycle)","land (C-N cycle)","ocean"),lty=1,lwd=2,col=c("black","darkorange","deepskyblue"),bty="n",cex=1.5)

plot(0,0,col="white",
     xlim=c(-20,30),xlab=TeX("Cumulative net $CO_2$ uptake (Pg C)"),
     ylim=c(-90,90),ylab="",cex.axis=1.5,cex.lab=1.5)
polygon(c(rev(apply(gamma_land_lat_all_df[!models$Ncycle,],2,quantile,na.rm=TRUE,0.975)),
          apply(gamma_land_lat_all_df[!models$Ncycle,],2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(apply(gamma_land_lat_all_df[models$Ncycle,],2,quantile,na.rm=TRUE,0.975)),
          apply(gamma_land_lat_all_df[models$Ncycle,],2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(255/256,140/256,0,0.1), border = NA)
polygon(c(rev(apply(gamma_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.975)),
          apply(gamma_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(0,191/256,255/256,0.1), border = NA)
lines(colMeans(gamma_land_lat_all_df[!models$Ncycle,],na.rm=TRUE),seq(-87.5,87.5,5),lwd=2)
lines(colMeans(gamma_land_lat_all_df[models$Ncycle,],na.rm=TRUE),seq(-87.5,87.5,5),lwd=2,col="darkorange")
lines(colMeans(gamma_ocean_lat_all_df,na.rm=TRUE),seq(-87.5,87.5,5),lwd=2,col="deepskyblue")
abline(v=0)
title("e) Carbon-climate feedback",adj=0,cex.main=1.5)

plot(0,0,col="white",
     xlim=c(-1,3),xlab=TeX("Cumulative net $CO_2$ uptake (Pg C)"),
     ylim=c(-90,90),ylab="",cex.axis=1.5,cex.lab=1.5)
polygon(c(rev(apply(epsilon_land_lat_all_df,2,quantile,na.rm=TRUE,0.975)),
          apply(epsilon_land_lat_all_df,2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(255/256,140/256,0,0.1), border = NA)
polygon(c(rev(apply(epsilon_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.975)),
          apply(epsilon_ocean_lat_all_df,2,quantile,na.rm=TRUE,0.025)),
        c(rev(seq(-87.5,87.5,5)), seq(-87.5,87.5,5)),
        col = rgb(0,191/256,255/256,0.1), border = NA)
lines(colMeans(epsilon_land_lat_all_df,na.rm=TRUE),seq(-87.5,87.5,by=5),col="darkorange",typ="l",lwd=2)
lines(colMeans(epsilon_ocean_lat_all_df,na.rm=TRUE),seq(-87.5,87.5,by=5),typ="l",lwd=2,col="deepskyblue")
abline(v=0)
title("f) Carbon-nitrogen feedback",adj=0,cex.main=1.5)

dev.off()

##### Figure 4 #####

pdf("Fig4.pdf",width=10.5,height=7)

par(mfrow=c(1,1),mar=c(5.1,5.1,2.1,5.1))
layout(matrix(c(1, 2, 1, 2, 3, 3), nrow = 2, ncol = 3))

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Cumulative net $CO_2$ uptake (Pg C)"),
     xlim=c(0,140),
     ylim=c(0,800),cex.axis=1.5,cex.lab=1.5)
polygon(c(rev(1:140), 1:140), c(rev(apply(apply(netAtmosLandCO2Flux_1pctCO2,1,cumsum),1,quantile,na.rm=TRUE,0.975)),
                                apply(apply(netAtmosLandCO2Flux_1pctCO2,1,cumsum),1,quantile,na.rm=TRUE,0.025)),
        col = rgb(255/256,140/256,0/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(apply(netAtmosLandCO2Flux_1pctCO2Ndep,1,cumsum),1,quantile,na.rm=TRUE,0.975)),
                                apply(apply(netAtmosLandCO2Flux_1pctCO2Ndep,1,cumsum),1,quantile,na.rm=TRUE,0.025)),
        col = rgb(255/256,140/256,0/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(apply(fgco2_1pctCO2,1,cumsum),1,quantile,na.rm=TRUE,0.975)),
                                apply(apply(fgco2_1pctCO2,1,cumsum),1,quantile,na.rm=TRUE,0.025)),
        col = rgb(0/256,191/256,255/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(apply(fgco2_1pctCO2Ndep,1,cumsum),1,quantile,na.rm=TRUE,0.975)),
                                apply(apply(fgco2_1pctCO2Ndep,1,cumsum),1,quantile,na.rm=TRUE,0.025)),
        col = rgb(0/256,191/256,255/256,0.1), border = NA)
lines(1:140,cumsum(colMeans(netAtmosLandCO2Flux_1pctCO2)),col="darkorange",lwd=2)
lines(1:140,cumsum(colMeans(fgco2_1pctCO2)),col="deepskyblue",lwd=2)
lines(1:140,cumsum(colMeans(netAtmosLandCO2Flux_1pctCO2Ndep)),col="darkorange",lwd=2,lty=2)
lines(1:140,cumsum(colMeans(fgco2_1pctCO2Ndep)),col="deepskyblue",lwd=2,lty=2)
legend("topleft",c("land 1pctCO2","ocean 1pctCO2","land 1pctCO2Ndep","ocean 1pctCO2Ndep"),
       col=c("darkorange","deepskyblue"),lty=c(1,1,2,2),lwd=2,bty="n",
       cex=1.5)
title("a",adj=0,cex.main=1.5)

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Cumulative net $N_2O$ emissions (Tg N)"),
     xlim=c(0,140),
     ylim=c(0,2000),cex.axis=1.5,cex.lab=1.5)
df<-rbind(cumsum(get("fN2O_MIROC-ES2L_1pctCO2")),
          cumsum(get("fN2O_MPI-ESM1-2-LR_1pctCO2")))
polygon(c(rev(1:140), 1:140), c(rev(apply(df,2,quantile,na.rm=TRUE,0.975)),apply(df,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(255/256,140/256,0/256,0.1), border = NA,xpd=T)
df<-rbind(cumsum(get("fN2O_MIROC-ES2L_1pctCO2Ndep")),
          cumsum(get("fN2O_MPI-ESM1-2-LR_1pctCO2Ndep")))
polygon(c(rev(1:140), 1:140), c(rev(apply(df,2,quantile,na.rm=TRUE,0.975)),apply(df,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(255/256,140/256,0/256,0.1), border = NA,xpd=T)
lines(1:140,cumsum(colMeans(fN2O_1pctCO2,na.rm=TRUE)),col="darkorange",lwd=2)
lines(1:140,cumsum(colMeans(fN2O_1pctCO2Ndep,na.rm=TRUE)),col="darkorange",lwd=2,lty=2)
lines(1:140,cumsum(fgn2o_1pctCO2[1,]),col="deepskyblue",lwd=2)
lines(1:140,cumsum(fgn2o_1pctCO2Ndep[1,]),col="deepskyblue",lwd=2,lty=2)
axis(side=4,at=seq(0,200,50)/(1/1000*273*44.01/(14.01*2)*12.01/44.01),seq(0,200,50),
     col="gray",col.axis="gray",xpd=TRUE,cex.axis=1.5)
#mtext(TeX("Cumulative net land-atmosphere $N_2O$ flux (Pg C)"),side=4,line=3,col="gray")
mtext(TeX("(Pg C equivalent)"),side=4,line=3,col="gray")
title("b",adj=0,cex.main=1.5)

CO2uptake_land<-rowSums(netAtmosLandCO2Flux_1pctCO2Ndep)-rowSums(netAtmosLandCO2Flux_1pctCO2)
N2Oemissions_land_C<-(rowSums(fN2O_1pctCO2Ndep)-rowSums(fN2O_1pctCO2))/1000*273*44.01/(14.01*2)*12.01/44.01
CO2uptake_ocean<-rowSums(fgco2_1pctCO2Ndep)-rowSums(fgco2_1pctCO2)
N2Oemissions_ocean_C<-(rowSums(fgn2o_1pctCO2Ndep)-rowSums(fgn2o_1pctCO2))/1000*273*44.01/(14.01*2)*12.01/44.01

bottoms <- c(0,mean(CO2uptake_land,na.rm=T)-mean(N2Oemissions_land_C,na.rm=T),0,
             0,mean(CO2uptake_ocean,na.rm=T)-mean(N2Oemissions_ocean_C,na.rm=T),0)
tops <- c(mean(CO2uptake_land),mean(CO2uptake_land),mean(CO2uptake_land,na.rm=T)-mean(N2Oemissions_land_C,na.rm=T),
          mean(CO2uptake_ocean),mean(CO2uptake_ocean),mean(CO2uptake_ocean,na.rm=T)-mean(N2Oemissions_ocean_C,na.rm=T))
xvals<-c(1:3,5:7)
colours<-rep(c("darkolivegreen3","orangered","gray"),2)

par(mar=c(2.1,6.1,4.1,2.1))
plot(0, 0, type = "n", xlim=c(0,8),ylim = c(-5, 40),
     xlab = "", ylab = "", xaxt = 'n', bty="n",cex.axis=1.5,cex.lab=1.5)
title(ylab=TeX("Cumulative net $CO_{2}$ uptake or $N_{2}O$ emissions driven by N deposition"),
      line=4,cex.lab=1.5)
title(ylab=TeX("(Pg C equivalent)"),line=2,cex.lab=1.5)
for (i in 1:length(bottoms)) {
  rect(xleft = xvals[i] - 0.4, ybottom = bottoms[i], xright = xvals[i] + 0.4, ytop = tops[i],
       col = colours[i], border=NA)
}
segments(1,quantile(CO2uptake_land,0.025),1,quantile(CO2uptake_land,0.975))
segments(2,mean(CO2uptake_land)-quantile(N2Oemissions_land_C,0.025,na.rm=T),
         2,mean(CO2uptake_land)-quantile(N2Oemissions_land_C,0.975,na.rm=T))
segments(5,quantile(CO2uptake_ocean,0.025),5,quantile(CO2uptake_ocean,0.975))

abline(h=0)
lines(c(4,4),c(0,40),lty=2,col="gray")
lines(c(1-0.4,2+0.4),rep(mean(CO2uptake_land),2))
lines(c(2-0.4,3+0.4),rep(mean(CO2uptake_land,na.rm=T)-mean(N2Oemissions_land_C,na.rm=T),2))
lines(c(5-0.4,6+0.4),rep(mean(CO2uptake_ocean),2))
lines(c(6-0.4,7+0.4),rep(mean(CO2uptake_ocean,na.rm=T)-mean(N2Oemissions_ocean_C,na.rm=T),2))
text(1,-2.5,TeX("$CO_{2}$"),srt=90,adj=1,xpd=T)
text(2,-2.5,TeX("$N_{2}O$"),srt=90,adj=1,xpd=T)
text(3,-2.5,TeX("$CO_{2}-N_{2}O$"),srt=90,adj=1,xpd=T)
text(5,-2.5,TeX("$CO_{2}$"),srt=90,adj=1,xpd=T)
text(6,-2.5,TeX("$N_{2}O$"),srt=90,adj=1,xpd=T)
text(7,-2.5,TeX("$CO_{2}-N_{2}O$"),srt=90,adj=1,xpd=T)
text(2,40,"Land",xpd=T,cex=1.5)
text(6,40,"Ocean",xpd=T,cex=1.5)
title("c",adj=0,cex.main=1.5)

dev.off()

##### Figure S1 ######

pdf("FigS1.pdf",width=10.5,height=7)

par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1))

plot(1:140,CO2[1:140,2],typ="l",
     xlab="Year",
     ylim=c(0,1200),ylab=TeX("$CO_2$ (ppm)"))
title("a",adj=0)

plot(1:140,Ndep_total,typ="l",
     xlab="Year",
     ylim=c(0,120),ylab=TeX("N deposition (Tg N $yr^{-1}$)"))
title("b",adj=0)

dev.off()

##### Figure S2 #####

pdf("FigS2.pdf",width=7,height=7)

par(mfrow=c(1,1),mar=c(5.1,5.1,4.1,2.1))

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Net N mineralisation rate (Tg N $yr^{-1}$)"),
     xlim=c(0,140),
     ylim=c(0,1000))
polygon(c(rev(1:140), 1:140), c(rev(apply(fNnetmin_1pctCO2,2,quantile,na.rm=TRUE,0.975)),apply(fNnetmin_1pctCO2,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fNnetmin_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.975)),apply(fNnetmin_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fNnetmin_1pctCO2_rad,2,quantile,na.rm=TRUE,0.975)),apply(fNnetmin_1pctCO2_rad,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(1,0,0,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fNnetmin_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.975)),apply(fNnetmin_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fNnetmin_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.975)),apply(fNnetmin_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)

lines(1:140,colMeans(fNnetmin_1pctCO2),col="black",lwd=2)
lines(1:140,colMeans(fNnetmin_1pctCO2_bgc),col="mediumblue",lwd=2)
lines(1:140,colMeans(fNnetmin_1pctCO2_rad),col="red",lwd=2)
lines(1:140,colMeans(fNnetmin_1pctCO2Ndep),col="black",lwd=2,lty=2)
lines(1:140,colMeans(fNnetmin_1pctCO2Ndep_bgc),col="mediumblue",lwd=2,lty=2)

legend("topleft",c("fully coupled","biogeochemically coupled","radiatively coupled",
                   "fully coupled + N deposition","biogeochemically coupled + N deposition"),
       col=c("black","mediumblue","red","black","mediumblue"),lty=c(1,1,1,2,2),lwd=2,
       bty="n")

dev.off()
##### Figure S3 #####

pdf("FigS3.pdf",width=7,height=7)

world_data <- ne_coastline(scale = "medium", returnclass = "sf")
world_vect <- vect(world_data)

par(mfrow=c(1,1))
plot(terra::rotate(mean(Ndep), left=TRUE), mar=c(3.1, 3.1, 2.1, 6.1), 
     axes=F,
     main=TeX("N deposition (kg N $m^{-2}$ $yr^{-1}$)"))
plot(world_vect,add=TRUE)

dev.off()

##### Figure S4 #####

pdf("FigS4.pdf",width=10.5,height=7)

par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1))

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Biological N fixation (Tg N $yr^{-1}$)"),
     xlim=c(0,140),
     ylim=c(0,300))
polygon(c(rev(1:140), 1:140), c(rev(apply(fBNF_1pctCO2,2,quantile,na.rm=TRUE,0.975)),apply(fBNF_1pctCO2,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fBNF_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.975)),apply(fBNF_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fBNF_1pctCO2_rad,2,quantile,na.rm=TRUE,0.975)),apply(fBNF_1pctCO2_rad,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(1,0,0,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fBNF_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.975)),apply(fBNF_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(fBNF_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.975)),apply(fBNF_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)

lines(1:140,colMeans(fBNF_1pctCO2),col="black",lwd=2)
lines(1:140,colMeans(fBNF_1pctCO2_bgc),col="mediumblue",lwd=2)
lines(1:140,colMeans(fBNF_1pctCO2_rad),col="red",lwd=2)
lines(1:140,colMeans(fBNF_1pctCO2Ndep),col="black",lwd=2,lty=2)
lines(1:140,colMeans(fBNF_1pctCO2Ndep_bgc),col="mediumblue",lwd=2,lty=2)

mean(fBNF_1pctCO2Ndep[1:3,1:20])
sd(fBNF_1pctCO2Ndep[1:3,1:20])
mean(fBNF_1pctCO2Ndep[1:3,120:140])
sd(fBNF_1pctCO2Ndep[1:3,120:140])

legend("topleft",c("fully coupled","biogeochemically coupled","radiatively coupled",
                   "fully coupled + N deposition","biogeochemically coupled + N deposition"),
       col=c("black","mediumblue","red","black","mediumblue"),lty=c(1,1,1,2,2),lwd=2,
       bty="n")

title("a",adj=0)

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Vegetation C:N ratio (kg C kg $N^{-1}$)"),
     xlim=c(0,140),
     ylim=c(0,300))
polygon(c(rev(1:140), 1:140), c(rev(apply(CNveg_1pctCO2,2,quantile,na.rm=TRUE,0.975)),apply(CNveg_1pctCO2,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(CNveg_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.975)),apply(CNveg_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(CNveg_1pctCO2_rad,2,quantile,na.rm=TRUE,0.975)),apply(CNveg_1pctCO2_rad,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(1,0,0,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(CNveg_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.975)),apply(CNveg_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(CNveg_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.975)),apply(CNveg_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)

lines(1:140,colMeans(CNveg_1pctCO2),col="black",lwd=2)
lines(1:140,colMeans(CNveg_1pctCO2_bgc),col="mediumblue",lwd=2)
lines(1:140,colMeans(CNveg_1pctCO2_rad),col="red",lwd=2)
lines(1:140,colMeans(CNveg_1pctCO2Ndep),col="black",lwd=2,lty=2)
lines(1:140,colMeans(CNveg_1pctCO2Ndep_bgc),col="mediumblue",lwd=2,lty=2)

mean(CNveg_1pctCO2Ndep[1:3,1:20])
sd(CNveg_1pctCO2Ndep[1:3,1:20])
mean(CNveg_1pctCO2Ndep[1:3,120:140])
sd(CNveg_1pctCO2Ndep[1:3,120:140])

title("b",adj=0)

dev.off()
##### Figure S5 #####

pdf("FigS5.pdf",width=7,height=7)

par(mfrow=c(1,1),mar=c(5.1,5.1,4.1,2.1))

plot(0,0,lwd=2,bty="n",col="white",
     xlab="Year",
     ylab=TeX("Temperature ($\\degree$C)"),
     xlim=c(0,140),
     ylim=c(10,20))
polygon(c(rev(1:140), 1:140), c(rev(apply(tas_1pctCO2,2,quantile,na.rm=TRUE,0.975)),apply(tas_1pctCO2,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(tas_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.975)),apply(tas_1pctCO2_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(tas_1pctCO2_rad,2,quantile,na.rm=TRUE,0.975)),apply(tas_1pctCO2_rad,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(1,0,0,0.1), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(tas_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.975)),apply(tas_1pctCO2Ndep,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(224/256,224/256,224/256,0.5), border = NA)
polygon(c(rev(1:140), 1:140), c(rev(apply(tas_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.975)),apply(tas_1pctCO2Ndep_bgc,2,quantile,na.rm=TRUE,0.025)),
        col = rgb(0,0,205/256,0.1), border = NA)

lines(1:140,colMeans(tas_1pctCO2),col="black",lwd=2)
lines(1:140,colMeans(tas_1pctCO2_bgc),col="mediumblue",lwd=2)
lines(1:140,colMeans(tas_1pctCO2_rad),col="red",lwd=2)
lines(1:140,colMeans(tas_1pctCO2Ndep),col="black",lwd=2,lty=2)
lines(1:140,colMeans(tas_1pctCO2Ndep_bgc),col="mediumblue",lwd=2,lty=2)

legend("topleft",c("fully coupled","biogeochemically coupled","radiatively coupled",
                   "fully coupled + N deposition","biogeochemically coupled + N deposition"),
       col=c("black","mediumblue","red","black","mediumblue"),lty=c(1,1,1,2,2),lwd=2,
       bty="n")
dev.off()

##### Figure S6 #####

pdf("FigS6.pdf",width=7,height=7)

# Plot feedbacks:

par(mfrow=c(3,2))

plot(gamma_MIROCES2L_land[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$gamma_L$ (Pg C $\\degree C^{-1}$)"),
     xlim=c(-100,0),ylim=c(-1,1))
points(gamma_UKESM10LL_land[140],0,pch=2)
points(gamma_MPIESM12LR_land[140],0,pch=3)
legend("topleft",c("MIROC-ES2L","UKESM1-0-LL","MPI-ESM1-2-LR"),pch=1:3,bty="n")
title("a",adj=0)

plot(gamma_MIROCES2L_ocean[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$gamma_O$ (Pg C $\\degree C^{-1}$)"),
     xlim=c(-100,0),ylim=c(-1,1))
points(gamma_UKESM10LL_ocean[140],0,pch=2)
points(gamma_MPIESM12LR_ocean[140],0,pch=3)
title("b",adj=0)

plot(beta_MIROCES2L_land[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$beta_L$ (Pg C $ppm^{-1}$)"),
     xlim=c(0,2),ylim=c(-1,1))
points(beta_UKESM10LL_land[140],0,pch=2)
points(beta_MPIESM12LR_land[140],0,pch=3)
title("c",adj=0)

plot(beta_MIROCES2L_ocean[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$beta_O$ (Pg C $ppm^{-1}$)"),
     xlim=c(0,2),ylim=c(-1,1))
points(beta_UKESM10LL_ocean[140],0,pch=2)
points(beta_MPIESM12LR_ocean[140],0,pch=3)
title("d",adj=0)

plot(epsilon_MIROCES2L_land[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$epsilon_L$ (Pg C Tg $N^{-1}$)"),
     xlim=c(-0.01,0.01),
     ylim=c(-1,1))
points(epsilon_UKESM10LL_land[140],0,pch=2)
points(epsilon_MPIESM12LR_land[140],0,pch=3)
points(epsilon_MIROCES2L_land_bgc[140],0,col="red")
points(epsilon_UKESM10LL_land_bgc[140],0,pch=2,col="red")
points(epsilon_MPIESM12LR_land_bgc[140],0,pch=3,col="red")
legend("topleft",col=c("black","red"),pch=20,c("Fully coupled method","Biogeochemically coupled method"),bty="n")
title("e",adj=0)

plot(epsilon_MIROCES2L_ocean[140],0,
     yaxt="n",ylab="",
     xlab=TeX("$epsilon_O$ (Pg C Tg $N^{-1}$)"),
     xlim=c(-0.01,0.01),
     ylim=c(-1,1))
points(epsilon_UKESM10LL_ocean[140],0,pch=2)
points(epsilon_MPIESM12LR_ocean[140],0,pch=3)
points(epsilon_MIROCES2L_ocean_bgc[140],0,col="red")
points(epsilon_UKESM10LL_ocean_bgc[140],0,pch=2,col="red")
points(epsilon_MPIESM12LR_ocean_bgc[140],0,pch=3,col="red")
title("f",adj=0)

dev.off()

##### Figure S7 #####

pdf("FigS7.pdf",width=10.5,height=7)

par(mfrow=c(1,2),mar=c(5.1,5.1,4.1,2.1))

plot(1:140,get("fNnetmin_MIROC-ES2L_1pctCO2"),col="purple",typ="l",
     ylab=TeX("Net N mineralisation (Tg N $yr^{-1}$)"),
     main="",xlab="Year",ylim=c(0,1000))
lines(1:140,get("fNnetmin_UKESM1-0-LL_1pctCO2"),col="pink")
lines(1:140,get("fNnetmin_MPI-ESM1-2-LR_1pctCO2"),col="blue")
lines(1:140,get("fNnetmin_MIROC-ES2L_1pctCO2Ndep"),col="purple",lty=2)
lines(1:140,get("fNnetmin_UKESM1-0-LL_1pctCO2Ndep"),col="pink",lty=2)
lines(1:140,get("fNnetmin_MPI-ESM1-2-LR_1pctCO2Ndep"),col="blue",lty=2)
legend("bottomleft",lty=c(1,1,1,1,2),col=c("purple","blue","pink","black","black"),
       c("MIROC-ES2L","MPI-ESM1-2-LR","UKESM1-0-LL","1pctCO2","1pctCO2Ndep"),bty="n")
title("a Fully coupled",adj=0)

plot(1:140,get("fNnetmin_MIROC-ES2L_1pctCO2-bgc"),col="purple",typ="l",
     ylab=TeX("Net N mineralisation (Tg N $yr^{-1}$)"),
     main="",xlab="Year",ylim=c(0,1000))
lines(1:140,get("fNnetmin_UKESM1-0-LL_1pctCO2-bgc"),col="pink")
lines(1:140,get("fNnetmin_MPI-ESM1-2-LR_1pctCO2-bgc"),col="blue")
lines(1:140,get("fNnetmin_MIROC-ES2L_1pctCO2Ndep-bgc"),col="purple",lty=2)
lines(1:140,get("fNnetmin_UKESM1-0-LL_1pctCO2Ndep-bgc"),col="pink",lty=2)
lines(1:140,get("fNnetmin_MPI-ESM1-2-LR_1pctCO2Ndep-bgc"),col="blue",lty=2)
title("b Biogeochemically coupled",adj=0)

dev.off()
