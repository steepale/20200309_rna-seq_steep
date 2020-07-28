#!/bin/R
# Author: Nicole Gay
# Feb 19 2020

# parallelize process by gene instead, which should always be maxing out a CPU
# saves one RData object per gene
impulse_optimization = foreach (i = 1:length(keep), .verbose=T, .combine = rbind) %dopar% {
        maximize_single_gene(keep[i], logfc_filt, outdir, .N=100, .maxiter=1000)
        #.N: number of different initializations
        #.maxiter: max number of iterations to optimize SSE
}



# from Sara Mostafavi
sigmoid = function(x, beta, t){
        return(1/(1+exp(-beta*(x-t))))
}

# from Sara Mostafavi (adjusted)
impulse_model_objective = function(theta, x, y){
        
        h0 = theta[1]
        h1 = theta[2]
        h2 = theta[3]
        t1 = theta[4]
        t2 = theta[5]
        beta1 = theta[6]
        beta2 = theta[7]
        
        # estimated values at given x:
        fits = get_impulse_vals(h0, h1, h2, t1, t2, beta1, beta2, x)
        f1 = fits$f
        S1 = fits$S1
        S2 = fits$S2
        s1 = fits$s1
        s2 = fits$s2
        
        # sum of squared error:
        # I don't think this should actually be multiplied by 2?
        obj = (1/2)*sum((f1-y)^2)
        
        # gradient of the impulse function:
        gh0 = ((1/h1)*s2*(1-S1))
        gh1 = (-(1/h1^2)*s1*s2 + (1/h1)*s1*S2 + (1/h1)*s2*S1)
        gh2 = ((1/h1)*s1*(1-S2))
        gt1 = (-(1/h1)*s2*((h1-h0)*S1*(1-S1)*beta1))
        gt2 = (-(1/h1)*s1*((h1-h2)*S2*(1-S2)*beta2))
        gb1 = -((1/h1)*s2*(h1-h0)*S1*(1-S1)*(t1-x) )
        gb2 = -(1/h1)*s1*(h1-h2)*S2*(1-S2)*(t2-x)
        
        grad = colSums(t(repmat((f1-y),7,1)) * matrix(c(gh0, gh1, gh2, gt1, gt2, gb1, gb2),
                                                      nrow=length(x),ncol=7,byrow = FALSE))
        
        return(c(obj,grad))
}

# from Sara Mostafavi (adjusted)
get_impulse_vals = function(h0, h1, h2, t1, t2, beta1, beta2, x){
        
        S1 = sigmoid(x, beta1, t1)
        S2 = sigmoid(x, beta2, t2)
        s1 = h0 + (h1 - h0)*S1
        s2 = h2 + (h1 - h2)*S2
        
        f = (1/h1)*s1*s2
        
        return(list(f=f, S1=S1, S2=S2, s1=s1, s2=s2))
}

estimate_starting_point = function(y,t){
        
        timecourse = y
        names(timecourse) = t
        
        h0_prior = 0 # Consider generating a more indicitive value (ACS)
        h1_prior = y[[which.max(abs(y))]]
        h2_prior = y[length(y)]
        t4h2_prior = t[length(y)]
        
        tmax = as.numeric(names(timecourse[which.max(abs(timecourse))]))
        t0 = -1 # Consider generating a more indicative value (ACS)
        
        t1_prior = (tmax - t0)/2
        # Adjust t2_prior to be more robust estimate
        t2_prior = (t4h2_prior-tmax)/2 # this is pretty inaccurate
        beta1_prior = (h1_prior - h0_prior)/(tmax - t0)
        
        # pick one of two starting points for beta2:
        # get closest time between tmax and 48h
        t_down = (t4h2_prior-tmax)/2
        # Consider taking the median for midval (ACS)
        # midval literally take the counts at middle time, first in group (ACS)
        midval = timecourse[[as.character(t[which.min(abs(t-t_down))])]]
        # Conditional statement determines if there is peak or trough
        if(which.min(abs(midval - c(h0_prior, h1_prior/2)))==1){
                beta2_prior = (h2_prior - h1_prior)/(t4h2_prior - tmax)
        }else{
                beta2_prior = -beta1_prior
        }
        
        # these don't have to be perfect, just not nonsensical:
        X = c(h0_prior, h1_prior, h2_prior, t1_prior, t2_prior, beta1_prior, beta2_prior)
        
        return(X)
}

random_starting_point = function(){
        h0 = runif(1, -5, 5)
        h1 = runif(1, -10, 10) # between -10 and 10
        h2 = runif(1, -10, 10) # between -10 and 10
        t1 = runif(1, -10, 20) # between -10 and 50
        t2 = runif(1, 1, 50) # between -10 and 50
        beta1 = runif(1, -5, 5) # between -10 and 10
        beta2 = runif(1, -5, 5) # between -10 and 10
        
        if(beta1*beta2 > 0){
                beta2 = -beta2
        }
        
        if(t2 < t1){
                ttmp = t2
                t2 = t1
                t1 = ttmp
        }
        
        X = c(h0, h1, h2, t1, t2, beta1, beta2)
        return(X)
}

#i <- 1
single_optimization = function(i, t, y, lim=500){
        if(i == 1){
                # always start with informed estimates since this seems to do well
                # you can tell which result corresponds to this one because h0_init = 0
                X = estimate_starting_point(y,t)
        }else{
                X = random_starting_point()
        }
        
        # from Sara Mustafavi:
        result = darch::minimize(X,impulse_model_objective,lim,t,y)
        all_res = c(X,
                    result[[1]],
                    result[[2]][length(result[[2]])], 
                    result[[3]])
        names(all_res) = c('h0_init','h1_init','h2_init','t1_init','t2_init','beta1_init','beta2_init',
                           'h0','h1','h2','t1','t2','beta1','beta2',
                           'sse',
                           'n_iter')
        return(all_res)
}



#by_gene_df <- by_gene_df_bk
# by_gene_df <- by_gene_df %>%
#         filter(SYMBOL_RAT == 'Arntl')
.gene <- 'Pdk4'
by_gene_df <- by_gene_df %>%
        filter(SYMBOL_RAT == 'Pdk4')
colnames(by_gene_df$data[[1]])
by_gene_df$data[[1]]$count
by_gene_df$data[[1]]$specimen.collection.t_exercise_hour

by_gene_df$data[[1]] %>%
        ggplot(aes(x = specimen.collection.t_exercise_hour, y = count)) +
        geom_point()
#outdir <- paste0(WD,'/test_dir')
maximize_single_gene = function(.gene, logfc_table, outdir, .N=100, .maxiter=1000){
        
        if(file.exists(sprintf('%s/%s_%s_%s.RData',outdir, .gene, .N, .maxiter))){
                 load(sprintf('%s/%s_%s_%s.RData',outdir, .gene, .N, .maxiter))
                 if(!grepl('gene',colnames(df))){
                         df$gene = .gene
                 }
                 return(df)
         }
        
        # # get logfc vals
        # log_vals = unlist(logfc_table[gene==.gene, .(log2FoldChange_0h, log2FoldChange_0.5h, log2FoldChange_1h,
        #                                              log2FoldChange_4h, log2FoldChange_7h, log2FoldChange_24h,
        #                                              log2FoldChange_48h)]) # expected column names in logfc_table
        
        #t = c(0, 0.5, 1, 1.5, 4.5, 7.5, 24.5, 48.5) # timepoints (0h post-exercise is actually 0.5 hours; make 0 the baseline)
        # Arrange the variables for the impulse model
        # by_gene_df
        
        # Create a dataframe that takes time and count data, then arrange by time
        count_df <- data.frame(TPE = by_gene_df$data[[1]]$specimen.collection.t_exercise_hour,
                   COUNTS = by_gene_df$data[[1]]$count) %>%
                arrange(TPE)
        
        t = count_df$TPE
        # The actual normalized gene counts
        
        #y = c(0, unname(log_vals)) # logFC vals (relative to 7h control) # first value is 0 because logFC between baseline and baseline is 0
        y = count_df$COUNTS
        # handle missing values
        # names(y) = t
        # y = y[!is.na(y)]
        # t = as.numeric(names(y))
        # y = unname(y)
        
        if(length(t) < 7){
                # can't fit. not enough time points (i.e. missing values)
                df = data.frame(matrix(nrow=0,ncol=17))
                colnames(df) = c("h0_init","h1_init","h2_init","t1_init","t2_init","beta1_init","beta2_init","h0","h1","h2","t1","t2","beta1","beta2","sse","n_iter","gene")
                return(df)
        }
        
        # get fit with good starting estimates
        prelim_optim = single_optimization(1, t, y, lim=.maxiter)
        if(prelim_optim['sse'] > 1){
                # stop looking
                df = data.frame(t(data.frame(prelim_optim)),row.names=1)
                df$gene = .gene
                save(df, file=sprintf('%s/%s_%s_%s.RData',outdir, .gene, .N, .maxiter))
                return(df)
        }
        
        # keep going
        prelim_optim2 = foreach (i = 2:.N) %do% {
                single_optimization(i, t, y, lim=.maxiter)
        }
        df = data.frame(do.call(rbind, prelim_optim2))
        df = data.frame(rbind(data.frame(t(data.frame(prelim_optim)),row.names=1), df))
        df$gene = .gene
        df = df[order(df$sse, decreasing = F),]
        save(df, file=sprintf('%s/%s_%s_%s.RData',outdir, .gene, .N, .maxiter))
        return(df)
}

##################################################################################################################################################################
## PLOTTING
##################################################################################################################################################################

# gene: Ensembl gene name
# impulse params: list of optimized parameters in this order: 'h0','h1','h2','t1','t2','beta1','beta2'
# log_vals: list of log2FCs in this order: 0h, 0.5h, 1h, 4h, 7h, 24h, 48h (post-exercise groups relative to controls)
plot_mostafavi_impulse = function(gene, impulse_params, log_vals, xlim = 48.5, err = NA, err2 = NA, errors = NA, .annot=T){
        
        # format logFC
        orig = data.frame(y=c(0,unname(log_vals)), x=c(0,0.5,1,1.5,4.5,7.5,24.5,48.5), point=1)
        if(!is.na(errors)){
                orig$se = c(0, unname(unlist(errors)))
        }
        
        # get fitted values from equations
        names(impulse_params) = c('h0','h1','h2','t1','t2','beta1','beta2')
        
        t = linspace(0,48,1000)
        y_fitted = get_impulse_vals(impulse_params[['h0']],
                                    impulse_params[['h1']],
                                    impulse_params[['h2']],
                                    impulse_params[['t1']],
                                    impulse_params[['t2']],
                                    impulse_params[['beta1']],
                                    impulse_params[['beta2']],
                                    t)$f
        
        fitted_vals = data.frame(x=t,y=y_fitted,point=0)
        if(!is.na(errors)){
                fitted_vals$se = NA_real_
        }
        fitted_vals = data.frame(rbind(fitted_vals, orig))
        
        # get the gene symbol
        if(!exists('.map')){
                .map = data.table(bitr(geneID = rownames(counts_round), fromType = 'ENSEMBL', toType = 'SYMBOL', OrgDb = org.Rn.eg.db, drop = TRUE))
        }
        if(gene %in% .map[,ENSEMBL]){
                title=sprintf('%s %s (%s %s)',gene, .map[ENSEMBL == gene, SYMBOL], TISSUE, SEX)
        }else{
                title=sprintf('%s (%s %s)',gene, TISSUE, SEX)
        }
        
        if(!is.na(err) & !is.na(err2)){
                ann = sprintf('beta1 = %s\nbeta2 = %s\nsse = %s\nscaled sse=%s',
                              round(impulse_params[['beta1']],2),
                              round(impulse_params[['beta2']],2),
                              signif(err, 3),
                              signif(err2, 3))
        }else if(!is.na(err) & is.na(err2)){
                ann = sprintf('beta1 = %s\nbeta2 = %s\nsse = %s',
                              round(impulse_params[['beta1']],2),
                              round(impulse_params[['beta2']],2),
                              signif(err, 3))
        }else{
                ann = sprintf('beta1 = %s\nbeta2 = %s',
                              round(impulse_params[['beta1']],2),
                              round(impulse_params[['beta2']],2))
        }
        
        ann_y = (max(fitted_vals$y,na.rm=T) + min(fitted_vals$y,na.rm=T))/2
        
        if(.annot){
                
                g = ggplot()+
                        geom_line(data=fitted_vals[fitted_vals$point==0,], aes(x=x,y=y), colour='gray40') +
                        geom_point(data=fitted_vals[fitted_vals$point==1,], aes(x=x,y=y)) +
                        theme_classic() +
                        geom_point() + # add original logFC values
                        annotate('text',x=xlim-15, y=ann_y,
                                 label=ann,
                                 hjust=0) +
                        labs(title=title, x='Time', y='logFC')
                
        }else{
                g = ggplot()+
                        geom_line(data=fitted_vals[fitted_vals$point==0,], aes(x=x,y=y), colour='gray40') +
                        geom_point(data=fitted_vals[fitted_vals$point==1,], aes(x=x,y=y)) +
                        theme_classic() +
                        geom_point() + # add original logFC values
                        labs(title=title, x='Time', y='logFC')
        }
        
        if(!is.na(errors)){
                g = g + geom_errorbar(data = fitted_vals[fitted_vals$point==1,],
                                      aes(x=x, ymin=y-se, ymax=y+se),
                                      width=.5, colour=tissuecols[[TISSUE_LAB]]) +
                        geom_point(data=fitted_vals[fitted_vals$point==1,], aes(x=x,y=y))
        }
        
        return(g)
}
