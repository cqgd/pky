library(data.table)
#' @import data.table
library(ggplot2)
library(gamlss)
require(gamlss.tr)
gen.trun(0,family="NBI",type="left",name=".0tr") #THIS HAS TO BE IN THE GLOBAL ENVIRONMENT, CANNOT BE WITHIN THE FIT FUNCTION!
NBI.0tr <<- trun(0,family="NBI",type="left",name=".0tr", local=FALSE) #Just local == FALSE is not enough to make it global. Just locals (which gen.trun generates in in peaky:::) somehow aren't enough although all of these are referenced explcitly with peaky::: below

#environment(NBI.0tr) <- globalenv() #asNamespace('peaky')
library(devtools)
install_github("cqgd/R2BGLiMS")
library(R2BGLiMS)
library(plyr)
library(parallel)
library(fastcluster)
library(WGCNA) #IN THIS ORDER!
allowWGCNAThreads()

note = function(logfile=NA,logtime=T,...){
  if(logtime==TRUE){time=paste0(format(Sys.time(),"%d-%m-%Y %H:%M:%OS"),"\n")}else{time=""}
  msg = paste0(...,"\n"); cat(paste0(time,msg))
  if(!is.na(logfile)){write(paste0(time,msg),file=logfile,append=TRUE,sep="")}
}




######################################################
### BIN.start
######################################################
#' Bins putative interactions by distance
#'
#' Places putative interactions into equally-sized bins based on the distances they span, and stores these on disk.
#'
#' @param interactions_file Path to (tsv or csv) file containing putative interactions. Its first line must state the columns names: baitID, preyID, N. Subsequent lines report one putative interaction each: a bait fragment ID, prey fragment ID and readcount.
#' @param fragments_file Path to (bed) file containing fragment information. Its first line must state the column names: chrom, chromStart, chromEnd, ID. Each subsequent line reports the chromosome, starting coordinate, ending coordinate and ID of a fragment.
#' @param output_dir Directory where the generated bins will be stored. Will be created if it does not exist.
#' @param bins Number of bins to place the interactions into.
#' @return List containing the output directory and an overview of bins.
#'
#' @examples
#' base = "/local/data/public/hic/demo"
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' \dontrun{bin_interactions(interactions_file, fragments_file, output_dir=bins_dir)}
#'
#' @export

bin_interactions  = function(interactions_file, fragments_file, output_dir, bins=50){
  L = paste0(output_dir,"/log_bins.txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("BINNING\n",file=L,append=FALSE,sep="")

  note(L,T, "Reading interactions from ",interactions_file)
  D = fread(interactions_file)
  if(ncol(D)<3 | all.equal(colnames(D)[1:3],c("baitID","preyID","N"))!=TRUE){note(L,T,"First three columns of interactions file should be named: baitID preyID N"); stop()}

  note(L,T, "Reading fragment information from ",fragments_file)
  fragments = fread(fragments_file)
  if(ncol(fragments)<4 | all.equal(colnames(fragments)[1:4],c("chrom","chromStart","chromEnd","ID"))!=TRUE){note(L,T,"First four columns of fragments file should be named: chrom chromStart chromEnd ID"); stop()}
  setnames(fragments,c("chr","start","end","ID"))

  note(L,T,"Calculating fragment characteristics...")
  print(fragments)
  print(is.data.table(fragments))
  fragments = fragments[,.(chr, mid=round((end+start)/2), length=end-start, ID)]

  note(L,T,"Adding fragment characteristics for baits...")
  D = merge(D,fragments[fragments$ID %in% D$baitID],by.x="baitID",by.y="ID",all.x=TRUE)
  setnames(D,paste0(ifelse(names(D) %in% names(fragments),"b.",""),names(D)))

  note(L,T,"Adding fragment characteristics for preys...")
  D = merge(D,fragments[fragments$ID %in% D$preyID],by.x="preyID",by.y="ID",all.x=TRUE)
  setnames(D,paste0(ifelse(names(D) %in% names(fragments),"p.",""),names(D)))

  #Ditch all trans here?

  note(L,T,"Calculating interaction distances...")
  D[b.chr==p.chr,dist:=(p.mid-b.mid)] #no clue, but this changes the number of rows somehow, with is for the whole dataset apparently, not the subset
  #D[b.chr==p.chr,dist:=(D[b.chr==p.chr,p.mid]-D[b.chr==p.chr,b.mid])]

  note(L,T,"Calculating total trans-chromosomal read counts for each bait...")
  trans = D[,.(b.trans=sum((b.chr!=p.chr)*N)),by=.(baitID, b.chr)]

  note(L,T,"Modelling those as a function of bait chromosome...")
  trans_model = gamlss(b.trans ~ as.factor(b.chr), data=trans, family=NBI)
  trans[,b.trans_res:=trans_model$residuals]
  D = merge(D,trans[,.(baitID, b.trans, b.trans_res)],by="baitID",all.x=TRUE)

  note(L,T,"Assigning ",bins," distance bins...")
  D[p.chr==b.chr, dist.bin:=as.factor(as.integer(cut_number(D[p.chr==b.chr,abs(dist)],n=bins)))]

  save_bin = function(Dbin,output_dir,bins){
    filename = paste0(output_dir,"/bin_",as.character(Dbin$dist.bin)[1],".rds")
    note(L,F,filename)
    saveRDS(Dbin, file=filename)
  }

  note(L,T,"Saving bin details to ",output_dir,"/bins.txt")
  overview = D[,.(dist.abs.min=min(abs(dist)),dist.abs.max=max(abs(dist)),interactions=.N),by=dist.bin][order(dist.abs.min),]
  fwrite(overview,file=paste0(output_dir,"/bins.txt"),sep=" ",na=0)

  note(L,T,"Saving binned interactions:")
  by(D, D$dist.bin, FUN=save_bin, output_dir=output_dir, bins=bins) #Extra 'custom' columns from the input file are saved as well.

  note(L,T,"Done.")
  return(list(output_dir=output_dir,overview=overview))
}

######################################################
### BIN.end
######################################################


######################################################
### FIT.start
######################################################
#' Obtain readcounts adjusted for several forms of technical and biological background noise
#'
#' Parametrizes a negative binomial null model for the readcounts in a given distance bin. Covariates are the distance between an interaction's bait and prey fragments, their lengths, and transchromosomal bait activity. Randomized quantile residuals computed from this model are taken as noise-adjusted readcounts and stored on disk.
#'
#' @param bins_dir Directory containing putative interactions that are binned by distance.
#' @param bin_index Index of the bin to be processed.
#' @param output_dir Directory where the parametrized model and adjusted readcounts will be stored. Will be created if it does not exist.
#' @param subsample_size Number of interactions based on which the null-model is parametrized.
#' @param gamlss_cycles GAMLSS maximum number of cycles for convergence (see gamlss::gamlss.control).
#' @param gamlss_crit GAMLSS convergence criterion (see gamlss::gamlss.control).
#' @return List containing the path to the output directory, the null model, and the adjusted read counts.
#'
#' @examples
#' # fits_dir = paste0(base,"/fits")
#'
#' for(bin_index in 1:8){
#'   \dontrun{model_bin(bins_dir,bin_index,output_dir=fits_dir)}
#' }
#' @export


model_bin = function(bins_dir,bin_index,output_dir,subsample_size=NA,gamlss_cycles=200,gamlss_crit=0.1,log_file=TRUE){
  if(log_file==TRUE){L = paste0(output_dir,"/log_fit_",bin_index,".txt")}else{L=NA}
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  if(!is.na(log_file)){write("FITTING\n",file=L,append=FALSE,sep="")}

  b=bin_index
  output_path_fit = paste0(output_dir,"/fit_",b,".rds")
  output_path_residuals = paste0(output_dir,"/residuals_",b,".rds")


  bin_path = paste0(bins_dir,"/bin_",b,".rds")
  note(L,F,"Loading bin from ",bin_path)
  bin = readRDS(bin_path)

  if(is.na(subsample_size)){
    subsample_size=nrow(bin)
    note(L,T,"No subsampling size provided, using all ",subsample_size," interactions for the null model regression...")
  }else{
    note(L,T,"Subsampling ",subsample_size,"/",nrow(bin)," interactions for the null model regression...")
  }
  subset = sample(1:nrow(bin),subsample_size)

  note(L,T,"Fitting with a maximum of ",gamlss_cycles," iterations...")
  control = gamlss.control(c.crit=gamlss_crit, n.cyc=gamlss_cycles)

  fit = gamlss(N ~ log(abs(dist)) + b.trans_res  + sqrt(b.length) + sqrt(p.length),
               data=bin[subset,], family=peaky:::NBI.0tr, sigma.formula = ~log(abs(dist)), control=control)

  if(is.gamlss(fit)){
    note(L,T,"Converged: ",ifelse(fit$converged,"YES","NO"),"\nIterations: ",fit$iter,"\n\nCoefficients:")
    note(L,F,paste0(names(coef(fit)),"\t",coef(fit)))
    note(L,F,"Saving fit to ",output_path_fit,"...")
    saveRDS(fit, output_path_fit)

    note(L,T,"Obtaining normalized randomized quantile residuals for the full dataset...")
    residuals_all = gamlss:::rqres(pfun="peaky:::pNBI.0tr", type="Discrete", ymin=1,
                                   y=bin$N,
                                   mu=predict(fit,what="mu",newdata=bin,type="response",data=bin[subset,]),
                                   sigma=predict(fit,what="sigma",newdata=bin,type="response",data=bin[subset,]))
    note(L,T,"Saving all residuals to ",output_path_residuals,"...")
    saveRDS(residuals_all, output_path_residuals)
  }else{
    note(L,T,"No fit object obtained.")
  }

  note(L,T,"Done.")
  return(list(output_dir=output_dir, fit=fit, residuals=residuals_all))
}

######################################################
### FIT.end
######################################################


######################################################
### SPLIT_BAITS.start
######################################################
#' Regroups putative interactions by bait
#'
#' Regroups putative interactions by bait, rather than by distance, to prepare them for parallel RJMCMC processing. Generates a separate file for each bait. Paths are stored in baitlist.txt, which serves as a to-do list for peaky().
#'
#' @param bins_dir Directory containing putative interactions that are binned by distance.
#' @param residuals_dir Directory where the adjusted read counts from each distance bin are stored.
#' @param output_dir Directory where all putative interactive will be stored, one file per bait. Will be created if it does not exist.
#' @param indices The array of distance bins of which the baits to be processed. These must all have had null models fitted.
#' @param plots Whether adjusted readcounts are to be plotted aganst distance and stored for each bait.
#' @return The output directory.
#'
#' @examples
#' base = "/local/data/public/hic/demo"
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' \dontrun{bin_interactions(interactions_file, fragments_file, output_dir=bins_dir)}
#'
#' @export


split_baits = function(bins_dir, residuals_dir, indices, output_dir, plots=TRUE, log_file=TRUE){
  if(log_file==TRUE){L = paste0(output_dir,"/log_split.txt")}else{L=NA}
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  if(!is.na(log_file)){write("SPLITTING BAITS\n",file=L,append=FALSE,sep="")}

  bin_files = paste0(bins_dir,"/bin_",indices,".rds")
  residual_files = paste0(residuals_dir,"/residuals_",indices,".rds")

  note(L,T,"Loading bins..."); note(L,F,bin_files)
  bins = sapply(bin_files,function(file){readRDS(file)},simplify = FALSE)

  note(L,T,"Loading residuals..."); note(L,F,residual_files)
  residuals = sapply(residual_files,function(file){readRDS(file)},simplify = FALSE)

  pvalues = function(residuals_set) {
    p.res = pnorm(abs(residuals_set),lower.tail=FALSE) * 2
    fdr.res = p.adjust(p.res,method="fdr")
    p.res.onesided = pnorm(residuals_set,lower.tail=FALSE)
    fdr.res.onesided = p.adjust(p.res.onesided,method="fdr")
    data.table(p.res,fdr.res,p.res.onesided,fdr.res.onesided)
  }

  note(L,T,"Calculating p-values and performing bin-wise FDR-adjustments...")
  ps = lapply(residuals,pvalues)

  D = cbind(rbindlist(bins),residual=unlist(residuals),rbindlist(ps))

  save_bait = function(baitset,output_dir){
    filename = paste0(output_dir,"/bait_",as.character(baitset$baitID)[1],".rds")
    note(L,F,filename)
    saveRDS(baitset, file=filename)
  }

  plot_bait = function(baitset, output_dir){
    filename = paste0(output_dir,"/bait_",as.character(baitset$baitID)[1],".pdf")
    plot = ggplot(data=baitset,aes(x=dist, y=residual)) +
      geom_hline(yintercept=0, color="grey") + geom_vline(xintercept=0, color="steelblue") +
      geom_point() +
      xlab("Distance from bait (bp)") + ylab("Residuals") + ggtitle(paste0("Bait ",as.character(baitset$baitID)[1]))
    ggsave(plot, filename = filename,width=10,height=4)
  }

  note(L,T,"Saving baits...")
  by(D,as.factor(D$baitID), FUN = save_bait, output_dir=output_dir)

  bait_paths = paste0(output_dir,"/bait_",unique(D$baitID),".rds")
  write(bait_paths,paste0(output_dir,"/baitlist.txt"),sep="\n")

  if(plots==TRUE){
    note(L,T,"Plotting baits...")
    by(D,as.factor(D$baitID), FUN = plot_bait, output_dir=output_dir)
  }

  note(L,T,"Done.")
  return(list(output_dir=output_dir,bait_paths=bait_paths))
}


######################################################
### SPLIT_BAITS.end
######################################################



######################################################
### RJMCMC.start
######################################################

# library(optparse)
# opt_list = list(
#   make_option(c("-b", "--baitsfile"), type="character", metavar="STR", default=NULL,
#               help="REQUIRED: File containg bait paths on separate line."),
#   make_option(c("-i", "--index"), type="integer", metavar="INT", default=NULL,
#               help="REQUIRED: Line number to read from the file."),
#   make_option(c("-r", "--rjmcmcdir"), type="character", metavar="DIR", default=NULL,
#               help="REQUIRED: Absolute path to the directory where rjMCMC results will be stored. No trailing slash.")
# )
#
# opt = parse_args(OptionParser(option_list=opt_list))
# print(opt)
#
# baitsfile = opt$baitsfile
# output_dir = opt$r

#' Regroups putative interactions by bait
#'
#' Regroups putative interactions by bait, rather than by distance, to prepare them for RJMCMC processing.
#'
#' @param baitlist Path to the list of baits to be processed.
#' @param index Which bait on the list to process, with 1 corresponding to the first one on the list.
#' @param output_dir Directory where all RJMCMC results will be stored. Will be created if it does not exist.
#' @param omega_power Expected decay of adjusted read counts around a truly interacting prey. See details.
#' @param min_interactions Whether adjusted readcounts are to be plotted aganst distance and stored for each bait.
#' @param log_file Whether or not a log-file is required.
#' @param min_interactions Minimum requirement for the number of prey fragments (and thus counts) associated with a bait, baits with fewer are skipped.
#' @return The output directory.
#'
#' Details.
#' The steepness of the function to be fitted to putative peaks, \eqn{\beta \exp{- \abs{\omega * d}}}.
#'
#' @examples
#' base = "/local/data/public/hic/demo"
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' \dontrun{bin_interactions(interactions_file, fragments_file, output_dir=bins_dir)}
#'
#' @export
#'
peaky = function(baitlist, index, output_dir, omega_power, min_interactions=20,log_file=TRUE,iterations=1e6){
  bait_path = scan(file=baitlist,what="character",skip=index-1,nline=1,sep="\n")
  P = readRDS(file=bait_path)

  if(!all(c("baitID","dist","residual")%in%names(P))|!(nrow(P)>=min_interactions)){
    stop("Column names and/or minimum interaction number requirements not met.")
  }

  P = P[order(P$dist),]
  genome = P$dist
  signal = P$residual
  baitID = unique(P$baitID)

  if(log_file==TRUE){L = paste0(output_dir,"/log_rjmcmc_",baitID,".txt")}else{L=NA}
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  if(!is.na(log_file)){write("PEAKY\n",file=L,append=FALSE,sep="")}

  note(L,T,"Loaded bait ",baitID," from path: ",bait_path,"...")

  note(L,T,"Constructing distance matrix with omega=10^-", omega_power,"...")

  omega = 10^omega_power

  dist_exp = function(offset, strength, omega){strength * exp(-(abs(offset*omega)))}

  distance_matrix = mapply(function(position, strength, genome, omega){dist_exp(genome-position, strength=strength, omega=omega)},
                           position=genome, strength=1,
                           MoreArgs=list(genome=genome, omega=omega),
                           SIMPLIFY=TRUE)

  extra.arguments <- list(
    "AlphaPriorMu" = 0,
    "AlphaPriorSd" = 0,
    "Alpha_Initial_Value" = 0,
    "GaussianResidual_Initial_Value" = 1,
    "GaussianResidualPrior_UnifArg1" = 0.99,
    "GaussianResidualPrior_UnifArg2" = 1.01,
    "GaussianResidualPriorFamily" = 1,
    "Add_Move_Probability"= 0.1,
    "Delete_Move_Probability" = 0.1,
    "Swap_Move_Probability" = 0.01
  )

  model.space.prior=list(list("a"=1, "b"=length(unique(P$preyID)), "Variables"=paste0("P",1:length(genome))))

  data=data.frame(cbind(distance_matrix,signal))
  colnames(data) = c(paste0("P",1:length(genome)),"S")
  results = R2BGLiMS(
    likelihood="Gaussian",
    data=data,
    outcome.var="S",
    model.space.priors = model.space.prior,
    extra.arguments = extra.arguments,
    seed=sample(1:1e6,1),
    n.iter=iterations
  )

  results_path = paste0(output_dir,"/rjmcmc_",baitID,".rds")
  note(L,T,"Saving RJMCMC results to ",results_path)
  saveRDS(results,file=results_path)
  note(L,T,"Done.")
  return(list(output_dir=output_dir,rjmcmc=results))
}

######################################################
### RJMCMC.end
######################################################


######################################################
### INTERPRET RJMCMC.start
######################################################

mppi = function(thin,zero_min,zero_max,given_present=FALSE){
  condition = thin<zero_min | thin>zero_max
  if(given_present==FALSE){
    mppi = colMeans(condition)
  }else{
    present = abs(thin)>1e-5
    mppi = colSums(condition & present) / colSums(present)
  }
  return(mppi)
}

get_modules_fast = function(thin){ #chris
  S = WGCNA::cor(thin[,which(thin[,lapply(.SD,sd)]!=0),with=FALSE])
  d = as.dist(1-abs(S))
  h = fastcluster::hclust(d)
  hcut = cutree(h,h=0.99)
  return(data.table(pos=names(hcut),module=hcut)) #slight modification
}

crset_threshold_fast = function(module,thin){
  modelspace_roi = thin[,lapply(.SD,function(col){abs(col)>1e-5}),.SDcols=module]
  modelspace_roi = modelspace_roi[rowSums(modelspace_roi)>=1] #for frequencies given 1+ interaction
  models = data.table(count(modelspace_roi))
  models = models[rev(order(freq)),][,pp:=freq/nrow(modelspace_roi)][,cumpp_lower:=shift(cumsum(pp),n=1,fill=0,type="lag")]
  thresholds = models[,lapply(.SD,function(col,cumpp){cumpp[which.max(col)]},cumpp=cumpp_lower),.SDcols=module]
  return(as.numeric(thresholds))
}

module_pp_fast = function(module, thin){
  modelspace_roi = thin[,lapply(.SD,function(col){abs(col)>1e-5}),.SDcols=module]
  pp = mean(as.logical(rowSums(modelspace_roi)))
  return(pp)
}

rjmcmc_list = function(rjmcmc_dir, filename=NULL){
  note(NA,T,"Looking for RJMCMC results in ",rjmcmc_dir)
  rjmcmc_dir_files = list.files(rjmcmc_dir)
  note(NA,F,"Found ",length(rjmcmc_dir_files)," files in total.")
  rjmcmc_files = sort(rjmcmc_dir_files[grepl("^rjmcmc_.*rds$",rjmcmc_dir_files)])
  rjmcmc_paths = paste0(rjmcmc_dir,"/",rjmcmc_files)
  note(NA,F,"Found ",length(rjmcmc_files)," files containing RJMCMC results.")
  if(is.null(filename)){filename=paste0(rjmcmc_dir,"/rjmcmclist.txt")}
  write(rjmcmc_paths,filename,sep="\n")
  return(list(filename=filename, rjmcmc_paths=rjmcmc_paths))
}


#' Interprets RJMCMC results
#'
#' Extracts the marginal posterior probabilities of contact and other relevant statistics from the RJMCMC results and merges these with bait-associated information.
#'
#' @param rjmcmclist Path to the list of RJMCMC results to be interpreted.
#' @param index Which RJMCMC result on the list to process, with 1 corresponding to the first one on the list.
#' @param output_dir Directory where the merged RJMCMC results and bait information will be stored, one file per bait. Will be created if it does not exist.
#' @param omega_power The same value as used when running peaky, i.e. the expected decay of adjusted read counts around a truly interacting prey.
#' @return The output directory.
#'
#' @examples
#' base = "/local/data/public/hic/demo"
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' \dontrun{bin_interactions(interactions_file, fragments_file, output_dir=bins_dir)}
#'
#' @export
#'
interpret_peaky = function(rjmcmclist, index, baits_dir, output_dir, omega_power){
  rjmcmc_path = scan(file=rjmcmclist,what="character",skip=index-1,nline=1)
  path_split = strsplit(rjmcmc_path,"_")[[1]]
  id = as.numeric(gsub(".rds","",path_split[length(path_split)]))
  bait_path = paste0(baits_dir,"/bait_",id,".rds")

  L=paste0(output_dir,"/log_baits_rjmcmc_",id,".txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("INTERPRETING RJMCMC RESULTS\n",file=L,append=FALSE,sep="")

  note(L,F,"Loading RJMCMC result and bait:\n",rjmcmc_path,"\n",bait_path)

  rjmcmc = readRDS(rjmcmc_path)
  D = readRDS(bait_path)
  setorder(D,dist)

  thin = as.matrix(rjmcmc@mcmc.output)
  loci = grep("^P[0-9]+$",colnames(thin),value=TRUE)
  thin = thin[,loci]

  ### ALPHAS

  note(L,T,"Calculating alphas...")
  alphas = rjmcmc@mcmc.output[,"alpha"]; rm(rjmcmc)

  ### POSTERIORS

  note(L,T,"Calculating posteriors...")
  zero=1e-5
  pars = data.table(zero_min=c(-zero,-Inf, -Inf, -zero, -zero),
                    zero_max=c(zero, zero, zero, Inf, Inf),
                    present=c(FALSE,FALSE,TRUE,FALSE,TRUE))
  posteriors = mapply(mppi,list(thin),pars$zero_min, pars$zero_max, pars$present)
  colnames(posteriors) = c("mppi","mppi_pos","mppi_pos_present","mppi_neg","mppi_neg_present")

  ### BETAS

  note(L,T,"Calculating betas means and medians...")

  b_mean = colMeans(thin); thin = data.table(thin)
  b_mean_present = as.numeric(thin[,lapply(.SD,function(col){mean(col[col>abs(1e-5)])})])
  b_median = as.numeric(thin[,lapply(.SD,median)])
  b_median_present = as.numeric(thin[,lapply(.SD,function(col){median(col[col>abs(1e-5)])})])

  ### PREDICTED

  gamma_power = -4.69897 #WHY IS THIS HARDCODED?
  genome = D$dist

  note(L,T,"Constructing distance matrix with omega=10^-", omega_power,"...")

  omega = 10^omega_power

  dist_exp = function(offset, strength, omega){strength * exp(-(abs(offset*omega)))}

  distance_matrix = mapply(function(position, strength, genome, omega){dist_exp(genome-position, strength=strength, omega=omega)},
                           position=genome, strength=1,
                           MoreArgs=list(genome=genome, omega=omega),
                           SIMPLIFY=TRUE)

  note(L,F,"Finding RJMCMC predictions based on mean betas...")
  predicted = distance_matrix%*%b_mean
  predicted_present = distance_matrix%*%as.numeric(b_mean_present)

  ### MODULES

  note(L,T,"Finding modules...")

  modules = get_modules_fast(thin)
  modules[,crset:=crset_threshold_fast(pos,thin),by=module]
  modules[,module_pp:=module_pp_fast(pos,thin),by=module]

  ### OUTPUT

  note(L,T,"Wrapping up...")

  res = data.table(posteriors,
                   b_mean, b_mean_present, b_median, b_median_present,
                   predicted=predicted, predicted_present=predicted_present)
  colnames(res) = c("rjmcmc","rjmcmc_pos","rjmcmc_pos_present","rjmcmc_neg","rjmcmc_neg_present",
                    "beta_mean","beta_mean_present","beta","beta_present",
                    "predicted","predicted_present")
  res[,pos:=loci]
  res = merge(res,modules,by="pos",all.x=TRUE)
  res = res[order(as.numeric(gsub("P","",pos))),]
  D = cbind(D,res)

  output_path = paste0(output_dir,"/bait_rjmcmc_",id,".rds")
  note(L,T,"Saving the bait with its RJMCMC results to",output_path)
  saveRDS(D,file = output_path)
  return(list(output_path=output_path,output=D))
}


######################################################
### INTERPRET RJMCMC.end
######################################################





######################################################
### PACKAGE BUILD
######################################################

#Save first....
if(exists("GO")){
   rm(GO)
   setwd("C:/Users/cq/Documents/peaky")
   document()
   install("../peaky"); library(peaky)
}

######################################################
### WALKTHROUGH
######################################################
#system.file("extdata", "fragments.bed", package = "peaky")
#system.file("extdata", "counts.tsv", package = "peaky")
# base = system.file("extdata",package="peaky")
#
# #BINNING
# base = "/local/data/public/hic/demo"
# interactions_file = paste0(base,"/counts.tsv")
# bins_dir = paste0(base,"/bins")
# fragments_file = paste0(base,"/fragments.bed")
#
# bin_interactions(interactions_file, fragments_file, output_dir=bins_dir)
#
# #FITTING
# fits_dir = paste0(base,"/fits")
#
# for(bin_index in 1:8){
#   model_bin(bins_dir,bin_index,output_dir=fits_dir)
# }
#
# #SPLITTING
# baits_dir = paste0(base,"/baits")
#
# split_baits(bins_dir,residuals_dir = fits_dir, indices=1:8, output_dir = baits_dir)
#
# #RJMCMC
# baitlist = paste0(baits_dir,"/baitlist.txt")
# rjmcmc_dir = paste0(base,"/rjmcmc")
# omega_power = 4.7
#
# for(i in 1:49){
#   peaky(baitlist,i,output_dir=rjmcmc_dir,omega_power=omega_power)
# }
#
# rjmcmc_list(rjmcmc_dir)
#
# #MPPI/MPPC
# rjmcmclist = paste0(rjmcmc_dir,"/rjmcmclist.txt")
# baits_rjmcmc_dir = paste0(base,"/baits_rjmcmc")
# interpret_peaky(rjmcmclist,1,omega_power,baits_dir,baits_rjmcmc_dir)
#
# ##QUESTIONS
# #remove p-values, modules? what about .rds's? what to return? merging chains? relabel bait ids because no support for ids > 10e3







