######################################################
### PACKAGE BUILD
######################################################

#Save first....
if(exists("GO")){
  library(devtools)
  library(roxygen2)
  rm(GO)
  setwd("~/Work/pky")
  document()
  install("../pky"); library(peaky)
}

#######################################################

#' @importFrom ggplot2 cut_number ggplot geom_hline geom_vline geom_point xlab ylab ggtitle aes ggsave
#' @import data.table R2BGLiMS gamlss gamlss.tr gamlss.dist
#Import for R2BGLiMS: Error in path.package("R2BGLiMS") : none of the packages are loaded

library(gamlss) #Complains about not finding NBI() if not loaded via the library
gamlss.tr::gen.trun(0,family="NBI",type="left",name=".0tr") #THIS HAS TO BE IN THE GLOBAL ENVIRONMENT, CANNOT BE WITHIN THE FIT FUNCTION!
NBI.0tr <<- gamlss.tr::trun(0,family="NBI",type="left",name=".0tr", local=FALSE) #Just local == FALSE is not enough to make it global. Just locals (which gen.trun generates in in peaky:::) somehow aren't enough although all of these are referenced explcitly with peaky::: below

note = function(logfile=NA,logtime=T,...){
  if(logtime==TRUE){time=paste0(format(Sys.time(),"%d-%m-%Y %H:%M:%OS"),"\n")}else{time=""}
  msg = paste0(...,"\n"); cat(paste0(time,msg))
  if(!is.na(logfile)){write(paste0(time,msg),file=logfile,append=TRUE,sep="")}
  return(msg)
}

######################################################
### BIN.start
######################################################

#' Bins putative interactions by distance
#'
#' Places putative interactions into equally-sized bins based on the distances they span.
#'
#' @param interactions Data table containing putative interactions. Columns called baitID, preyID, N, storing the bait fragment ID, prey fragment ID and readcount, respectively.
#' @param fragments Data table containing fragment information. Columns called chrom, chromStart, chromEnd, ID, storing the chromosome, starting coordinate, ending coordinate and ID of a fragment, respectively.
#' @param bins Number of bins to place the interactions into.
#' @param log_file Path to a log file.
#' @return List containing the binned interactions ($interactions) and an overview of bins ($bins).
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#'
#' interactions_file = paste0(base,"/counts.tsv")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' interactions = data.table::fread(interactions_file)
#' fragments = data.table::fread(fragments_file)
#'
#' \donttest{BI = bin_interactions(interactions, fragments, bins=5)}
#' print(BI)
#'
#' @export

bin_interactions = function(interactions, fragments, bins=5, log_file=NA){
  D = interactions; rm(interactions)
  L = log_file; rm(log_file)

  if(ncol(D)<3 | all.equal(colnames(D)[1:3],c("baitID","preyID","N"))!=TRUE){note(L,T,"First three columns of interactions file should be named: baitID preyID N"); stop()}
  if(ncol(fragments)<4 | all.equal(colnames(fragments)[1:4],c("chrom","chromStart","chromEnd","ID"))!=TRUE){note(L,T,"First four columns of fragments file should be named: chrom chromStart chromEnd ID"); stop()}
  setnames(fragments,c("chr","start","end","ID"))

  note(L,T,"Calculating fragment characteristics...")
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
  trans_model = gamlss(b.trans ~ as.factor(b.chr), data=trans, family=gamlss.dist::NBI)
  trans[,b.trans_res:=trans_model$residuals]
  D = merge(D,trans[,.(baitID, b.trans, b.trans_res)],by="baitID",all.x=TRUE)

  note(L,T,"Assigning ",bins," distance bins...")
  D[p.chr==b.chr, dist.bin:=as.factor(as.integer(cut_number(D[p.chr==b.chr,abs(dist)],n=bins)))]

  overview = D[,.(dist.abs.min=min(abs(dist)),dist.abs.max=max(abs(dist)),interactions=.N),by=dist.bin][order(dist.abs.min),]

  note(L,T,"Done.")
  return(list(interactions=D,bins=overview))
}

#' Wrapper for bin_interactions() that interacts with the filesystem
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
#' base = system.file("extdata",package="peaky")
#'
#' interactions_file = paste0(base,"/counts.tsv")
#' fragments_file = paste0(base,"/fragments.bed")
#' bins_dir = paste0(base,"/bins")
#' \donttest{
#' BI = bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir)
#' print(BI)
#' }
#' @export

bin_interactions_fs  = function(interactions_file, fragments_file, output_dir, bins=5){
  L = paste0(output_dir,"/log_bins.txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("BINNING\n",file=L,append=FALSE,sep="")

  note(L,T, "Reading interactions from ",interactions_file)
  D = fread(interactions_file)

  note(L,T, "Reading fragment information from ",fragments_file)
  fragments = fread(fragments_file)

  BI = bin_interactions(D, fragments, bins, L)

  save_bin = function(Dbin,output_dir,bins){
    filename = paste0(output_dir,"/bin_",as.character(Dbin$dist.bin)[1],".rds")
    note(L,F,filename)
    saveRDS(Dbin, file=filename)
  }

  note(L,T,"Saving binned interactions:")
  by(BI$interactions, BI$interactions$dist.bin, FUN=save_bin, output_dir=output_dir, bins=bins) #Extra 'custom' columns from the input file are saved as well.

  note(L,T,"Saving bin details to ",output_dir,"/bins.txt")
  fwrite(BI$bins,file=paste0(output_dir,"/bins.txt"),sep=" ",na=0)

  return(list(output_dir=output_dir, interactions=BI$interactions, bins=BI$bins))
}


######################################################
### BIN.end
######################################################


######################################################
### FIT.start
######################################################
#' Obtain readcounts adjusted for several forms of technical and biological background noise
#'
#' Parametrizes a negative binomial null model for the readcounts in a given distance bin. Covariates are the distance between an interaction's bait and prey fragments, their lengths, and transchromosomal bait activity. Randomized quantile residuals computed from this model are taken as noise-adjusted readcounts.
#'
#' @param bin Data table containing putative interactions in the same distance bin.
#' @param subsample_size Number of interactions based on which the null-model is parametrized. By default, all are used.
#' @param gamlss_cycles GAMLSS maximum number of cycles for convergence (see gamlss::gamlss.control).
#' @param gamlss_crit GAMLSS convergence criterion (see gamlss::gamlss.control).
#' @param log_file Path to a log file.
#' @return List containing the fitted null model ($fit) and the adjusted readcounts ($residuals).
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#'
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' interactions = data.table::fread(interactions_file)
#' fragments = data.table::fread(fragments_file)
#'
#' \donttest{
#' BI = bin_interactions(interactions, fragments, bins=5)
#' BM = model_bin(BI$interactions[dist.bin==2,])
#'
#' print(BM)
#' plot(BM$fit)
#' }
#'
#' @export

model_bin = function(bin, subsample_size=NA, gamlss_cycles=200, gamlss_crit=0.1, log_file=NA){
  L = log_file; rm(log_file)

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
               data=bin[subset,], family=NBI.0tr, sigma.formula = ~log(abs(dist)), control=control)

  if(is.gamlss(fit)){
    note(L,T,"Converged: ",ifelse(fit$converged,"YES","NO"),"\nIterations: ",fit$iter,"\n\nCoefficients:")
    note(L,F,paste0(names(coef(fit)),"\t",coef(fit)))

    note(L,T,"Obtaining normalized randomized quantile residuals for the full dataset...")
    residuals_all = gamlss:::rqres(pfun="peaky:::pNBI.0tr", type="Discrete", ymin=1,
                                   y=bin$N,
                                   mu=predict(fit,what="mu",newdata=bin,type="response",data=bin[subset,]),
                                   sigma=predict(fit,what="sigma",newdata=bin,type="response",data=bin[subset,]))
    return(list(fit=fit, residuals=residuals_all))
  }else{
    note(L,T,"No fit object obtained.")
  }
}

#' Wrapper for model_bin() that interacts with the filesystem
#'
#' Parametrizes a negative binomial null model for the readcounts in a given distance bin. Covariates are the distance between an interaction's bait and prey fragments, their lengths, and transchromosomal bait activity. Randomized quantile residuals computed from this model are taken as noise-adjusted readcounts and stored on disk.
#'
#' @param bins_dir Directory containing putative interactions that are binned by distance.
#' @param bin_index Index of the bin to be processed.
#' @param output_dir Directory where the parametrized model and adjusted readcounts will be stored. Will be created if it does not exist.
#' @param subsample_size Number of interactions based on which the null-model is parametrized. By default, all are used.
#' @param gamlss_cycles GAMLSS maximum number of cycles for convergence (see gamlss::gamlss.control).
#' @param gamlss_crit GAMLSS convergence criterion (see gamlss::gamlss.control).
#' @return List containing the path to the output directory, the null model, and the adjusted read counts.
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' \donttest{bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir)}
#'
#' fits_dir = paste0(base,"/fits")
#'
#' for(bin_index in 1:8){
#'   \donttest{model_bin_fs(bins_dir,bin_index,output_dir=fits_dir)}
#' }
#' @export

model_bin_fs = function(bins_dir,bin_index,output_dir,subsample_size=NA,gamlss_cycles=200,gamlss_crit=0.1){
  L = paste0(output_dir,"/log_fit_",bin_index,".txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("FITTING\n",file=L,append=FALSE,sep="")

  b=bin_index
  output_path_fit = paste0(output_dir,"/fit_",b,".rds")
  output_path_residuals = paste0(output_dir,"/residuals_",b,".rds")

  bin_path = paste0(bins_dir,"/bin_",b,".rds")
  note(L,F,"Loading bin from ",bin_path)
  bin = readRDS(bin_path)

  BM = model_bin(bin, subsample_size=subsample_size, gamlss_cycles=gamlss_cycles, gamlss_crit=gamlss_crit, log_file=L)

  if(is.gamlss(BM$fit)){
    note(L,F,"Saving fit to ",output_path_fit,"...")
    saveRDS(BM$fit, output_path_fit)

    note(L,T,"Saving all residuals to ",output_path_residuals,"...")
    saveRDS(BM$residuals, output_path_residuals)
  }

  return(list(output_dir=output_dir, fit=BM$fit, residuals=BM$residuals))
}

######################################################
### FIT.end
######################################################


######################################################
### SPLIT_BAITS.start
######################################################
#' Regroups putative interactions by bait and calculates p-values under the negative binomial model
#'
#' Merges bin-wise data, calculates the p-values for putative interactions under the negative binomial model. Regroups putative interactions by bait rather than by distance bin.
#'
#' @param bins List of data tables containing putative interactions that are binned by distance.
#' @param residuals List containing residuals (adjusted read counts) for each bin (matching the interaction order).
#' @param log_file Path to a log file.
#' @return Data table with p-values, residuals and putative interaction attributes that can be sorted or split by bait ID.
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' fragments_file = paste0(base,"/fragments.bed")
#' interactions = data.table::fread(interactions_file)
#' fragments = data.table::fread(fragments_file)
#'
#' BI = bin_interactions(interactions, fragments, bins=5)
#' models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=1000)
#' residuals = lapply(models, "[[", "residuals")
#' bins = split(BI$interactions, BI$interactions$dist.bin)
#'
#' split_baits(bins, residuals)
#'
#' @export


split_baits = function(bins, residuals, log_file=NA){
  bins = rbindlist(bins)
  L = log_file; rm(log_file)

  pvalues = function(residuals_set) {
    p.res = pnorm(abs(residuals_set),lower.tail=FALSE) * 2
    fdr.res = p.adjust(p.res,method="fdr")
    p.res.onesided = pnorm(residuals_set,lower.tail=FALSE)
    fdr.res.onesided = p.adjust(p.res.onesided,method="fdr")
    data.table(p.res,fdr.res,p.res.onesided,fdr.res.onesided)
  }

  note(L,T,"Calculating p-values and performing bin-wise FDR-adjustments...")
  ps = rbindlist(lapply(residuals,pvalues))

  residuals = unlist(residuals)
  D = cbind(bins,residual=residuals,ps)
  return(D)

}

#' Wrapper for split_baits() that interacts with the filesystem
#'
#' Calculates p-values under the negative binomial model. Subsequently regroups putative interactions by bait, rather than by distance, to prepare them for parallel RJMCMC processing. Generates a separate file for each bait. Paths are stored in baitlist.txt, which serves as a to-do list for peaky().
#'
#' @param bins_dir Directory containing putative interactions that are binned by distance.
#' @param residuals_dir Directory where the adjusted read counts from each distance bin are stored.
#' @param output_dir Directory where all putative interactive will be stored, one file per bait. Will be created if it does not exist.
#' @param indices Indices of distance bins whose baits are processed. These must all have had null models fitted.
#' @param plots Whether adjusted readcounts are to be plotted aganst distance and stored for each bait.
#' @return The output directory.
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir)
#'
#' fits_dir = paste0(base,"/fits")
#'
#' for(bin_index in 1:5){
#'   model_bin_fs(bins_dir,bin_index,output_dir=fits_dir,subsample_size=1000)
#' }
#'
#' baits_dir = paste0(base,"/baits")
#'
#' split_baits_fs(bins_dir,residuals_dir = fits_dir, indices=1:5, output_dir = baits_dir)
#'
#' @export

split_baits_fs = function(bins_dir, residuals_dir, indices, output_dir, plots=TRUE){
  L = paste0(output_dir,"/log_split.txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("SPLITTING BAITS\n",file=L,append=FALSE,sep="")

  bin_files = paste0(bins_dir,"/bin_",indices,".rds")
  residual_files = paste0(residuals_dir,"/residuals_",indices,".rds")

  note(L,T,"Loading bins..."); note(L,F,bin_files)
  bins = lapply(bin_files,readRDS)
  #bins = rbindlist(bins)

  note(L,T,"Loading residuals..."); note(L,F,residual_files)
  residuals = lapply(residual_files,readRDS)
  #residuals = unlist(residuals)

  D = split_baits(bins, residuals, log_file=L)

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
#' Obtains peak locations and heights based on adjusted readcounts
#'
#' Fits additive models of peaks of varying strengths in various locations to the adjusted readcounts via RJMCMC, and stores these models on disk.
#'
#' @param bait Data table containing the putative interactions of a bait, having the columns 'baitID', 'dist', and 'residual'. These report the bait ID, its distance to putative preys, and the adjusted readcounts for its interactions with them, respectively.
#' @param omega_power Expected decay of adjusted read counts around a truly interacting prey. See details.
#' @param iterations Number of models to parametrize. Greated numbers should lead to increased reproducibility.
#' @param log_file Path to a log file.
#' @return The output directory.
#'
#' Details.
#' The steepness of the function to be fitted to putative peaks is determined by \eqn{\omega} according to \eqn{\beta \exp{- \abs{\omega * d}}}, where \beta represents peak height and \eqn{d} the distance from the center of the peak in bp.
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' fragments_file = paste0(base,"/fragments.bed")
#' interactions = data.table::fread(interactions_file)
#' fragments = data.table::fread(fragments_file)
#'
#' BI = bin_interactions(interactions, fragments, bins=5)
#' models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=1000)
#' residuals = lapply(models, "[[", "residuals")
#' bins = split(BI$interactions, BI$interactions$dist.bin)
#'
#' BTS = split_baits(bins, residuals)
#'
#' peaky(BTS[baitID==618421], omega_power=4.7, iterations=10e3)
#' @export

peaky = function(bait, omega_power, iterations=1e6, min_interactions=20, log_file=NA){
  L=log_file; rm(log_file)
  P=bait; rm(bait)

  if(!all(c("baitID","dist","residual")%in%names(P))){
    stop("Column names ought to include baitID, dist and residual")
  }else if(!(nrow(P)>=min_interactions)){
    stop("Minimum number of putative interactions (i.e. peak locations) not present.")
  }

  P = P[order(P$dist),]
  genome = P$dist
  signal = P$residual
  baitID = unique(P$baitID)

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

  return(results)
}

#' Wrapper for peaky() that interacts with the filesystem
#'
#' Fits additive models of peaks of varying strengths in various locations to the adjusted readcounts via RJMCMC, and stores these models on disk.
#'
#' @param baitlist Path to the list of baits to be processed.
#' @param index Which bait on the list to process, with 1 corresponding to the first one on the list.
#' @param output_dir Directory where all RJMCMC results will be stored. Will be created if it does not exist.
#' @param omega_power Expected decay of adjusted read counts around a truly interacting prey. See details.
#' @param min_interactions Minimum requirement for the number of prey fragments (and thus counts) associated with a bait, baits with fewer are skipped.
#' @param iterations Number of models to parametrize. Greated numbers should lead to increased reproducibility.
#' @return List containing the model output directory and the models themselves.
#'
#' Details.
#' The steepness of the function to be fitted to putative peaks is determined by \eqn{\omega} according to \eqn{\beta \exp{- \abs{\omega * d}}}, where \beta represents peak height and \eqn{d} the distance from the center of the peak in bp.
#'
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir)
#'
#' fits_dir = paste0(base,"/fits")
#'
#' for(bin_index in 1:5){
#'   model_bin_fs(bins_dir,bin_index,output_dir=fits_dir,subsample_size=1000)
#' }
#'
#' baits_dir = paste0(base,"/baits")
#'
#' split_baits_fs(bins_dir,residuals_dir = fits_dir, indices=1:5, output_dir = baits_dir)
#'
#' baitlist = paste0(baits_dir,"/baitlist.txt")
#' rjmcmc_dir = paste0(base,"/rjmcmc")
#' omega_power = 4.7
#' \donttest{
#' for(i in 1:3){
#'  peaky_fs(baitlist,i,output_dir=rjmcmc_dir,omega_power=omega_power)
#' }
#' }
#'
#' @export

peaky_fs = function(baitlist, index, output_dir, omega_power, iterations=1e6, min_interactions=20){
  bait_path = scan(file=baitlist,what="character",skip=index-1,nline=1,sep="\n")
  P = readRDS(file=bait_path)
  baitID = unique(P$baitID)

  L = paste0(output_dir,"/log_rjmcmc_",baitID,".txt")
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive=TRUE)}
  write("PEAKY\n",file=L,append=FALSE,sep="")


  note(L,T,"Loaded bait ",baitID," from path: ",bait_path,"...")

  PKY = peaky(bait=P, omega_power=omega_power, iterations=iterations, min_interactions=min_interactions, log_file=L)

  results_path = paste0(output_dir,"/rjmcmc_",baitID,".rds")
  note(L,T,"Saving RJMCMC results to ",results_path)
  saveRDS(PKY,file=results_path)
  note(L,T,"Done.")
  return(list(output_dir=output_dir,rjmcmc=PKY))
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

#' Create a list of peaky model files
#'
#' Reports and lists files that contain peaky models in preparation for their interpretation with interpret_peaky_fs()
#'
#' @param Path to a folder containing the results produced by peaky_fs()
#' @return List containing the path to a report of all results found ($filename), and the paths to the results themselves ($rjmcmc_paths.
#'
#' @export

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


#' Interprets models fitted with peaky
#'
#' Extracts the marginal posterior probabilities of contact and other relevant statistics from the RJMCMC results and merges these with bait-associated information.
#'
#' @param bait Data table containing the putative interactions of a bait, having the columns 'baitID', 'dist', and 'residual'. These report the bait ID, its distance to putative preys, and the adjusted readcounts for its interactions with them, respectively.
#' @param peaks The models fitted by peaky.
#' @param omega_power The same value as used when running peaky, i.e. the expected decay of adjusted read counts around a truly interacting prey. See details.
#' @return A data table containing bait-associated information, posterior probabilities of contact and other statistics.
#'
#' Details.
#' The steepness of the function to be fitted to putative peaks is determined by \eqn{\omega} according to \eqn{\beta \exp{- \abs{\omega * d}}}, where \eqn{\beta} represents peak height and \eqn{d} the distance from the center of the peak in bp.
#
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' fragments_file = paste0(base,"/fragments.bed")
#' interactions = data.table::fread(interactions_file)
#' fragments = data.table::fread(fragments_file)
#'
#' BI = bin_interactions(interactions, fragments, bins=5)
#' models = by(BI$interactions, BI$interactions$dist.bin, model_bin, subsample_size=1000)
#' residuals = lapply(models, "[[", "residuals")
#' bins = split(BI$interactions, BI$interactions$dist.bin)
#'
#' BTS = split_baits(bins, residuals)
#'
#' relevant_bait = BTS[baitID==618421]
#' omega_power=4.7
#' PKS = peaky(relevant_bait, omega_power, iterations=1e6)
#'
#' interpret_peaky(relevant_bait, PKS, omega_power)
#' @export

interpret_peaky = function(bait, peaks, omega_power, log_file=NA){
  L = log_file; rm(log_file)
  D = bait; rm(bait)
  rjmcmc = peaks; rm(peaks)

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

  note(L,T,"Wrapping up...")

  res = data.table(posteriors,
                   b_mean, b_mean_present, b_median, b_median_present,
                   predicted=predicted, predicted_present=predicted_present)
  colnames(res) = c("rjmcmc","rjmcmc_pos","rjmcmc_pos_present","rjmcmc_neg","rjmcmc_neg_present",
                    "beta_mean","beta_mean_present","beta","beta_present",
                    "predicted","predicted_present")
  res[,pos:=loci]
  #res = merge(res,modules,by="pos",all.x=TRUE)
  res = res[order(as.numeric(gsub("P","",pos))),]
  D = cbind(D,res)

}

#' Wrapper for interpret_peaky() that interacts with the filesystem
#'
#' Extracts the marginal posterior probabilities of contact and other relevant statistics from the RJMCMC results and merges these with bait-associated information.
#'
#' @param rjmcmclist Path to the list of RJMCMC results to be interpreted.
#' @param index Which RJMCMC result on the list to process, with 1 corresponding to the first one on the list.
#' @param output_dir Directory where the merged RJMCMC results and bait information will be stored, one file per bait. Will be created if it does not exist.
#' @param omega_power The same value as used when running peaky, i.e. the expected decay of adjusted read counts around a truly interacting prey. See details.
#' @return The output directory.
#'
#' Details.
#' The steepness of the function to be fitted to putative peaks is determined by \eqn{\omega} according to \eqn{\beta \exp{- \abs{\omega * d}}}, where \eqn{\beta} represents peak height and \eqn{d} the distance from the center of the peak in bp.
#
#' @examples
#' base = system.file("extdata",package="peaky")
#' interactions_file = paste0(base,"/counts.tsv")
#' bins_dir = paste0(base,"/bins")
#' fragments_file = paste0(base,"/fragments.bed")
#'
#' bin_interactions_fs(interactions_file, fragments_file, output_dir=bins_dir)
#'
#' fits_dir = paste0(base,"/fits")
#'
#' for(bin_index in 1:5){
#'   model_bin_fs(bins_dir,bin_index,output_dir=fits_dir,subsample_size=1000)
#' }
#'
#' baits_dir = paste0(base,"/baits")
#'
#' split_baits_fs(bins_dir,residuals_dir = fits_dir, indices=1:5, output_dir = baits_dir)
#'
#' baitlist = paste0(baits_dir,"/baitlist.txt")
#' rjmcmc_dir = paste0(base,"/rjmcmc")
#' omega_power = 4.7
#' \donttest{
#' for(i in 1:3){
#'  peaky_fs(baitlist,i,output_dir=rjmcmc_dir,omega_power=omega_power)
#' }
#'
#' rjmcmc_list(rjmcmc_dir)
#'
#' rjmcmclist = paste0(rjmcmc_dir,"/rjmcmclist.txt")
#' baits_rjmcmc_dir = paste0(base,"/baits_rjmcmc")
#' interpret_peaky_fs(rjmcmclist,1,omega_power,baits_dir,baits_rjmcmc_dir)
#' }
#' @export

interpret_peaky_fs = function(rjmcmclist, index, baits_dir, output_dir, omega_power){
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

  ### OUTPUT

  D = interpret_peaky(bait=D, peaks=rjmcmc, omega_power=omega_power, log_file=L)

  output_path = paste0(output_dir,"/bait_rjmcmc_",id,".rds")
  note(L,T,"Saving the bait with its RJMCMC results to",output_path)
  saveRDS(D,file = output_path)
  return(list(output_path=output_path,output=D))
}


######################################################
### INTERPRET RJMCMC.end
######################################################








