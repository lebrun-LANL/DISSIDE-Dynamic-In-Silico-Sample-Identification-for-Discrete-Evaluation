# Written by Erick LeBrun elebrun@lanl.gov
# Additional development by Chien-Chi Lo, Paul Li, and Jessica Salguero
# Version 1.5
# Date: 8/9/23
# LANL O4711

# Load Dependencies
suppressMessages(require("vegan"))
suppressMessages(require("ggplot2"))
suppressMessages(require("plotly"))
suppressMessages(require("htmlwidgets"))
suppressMessages(require("RColorBrewer"))
suppressMessages(require("reticulate"))
path_to_python <- "/opt/anaconda3/envs/disside/bin/python"
reticulate::use_python(path_to_python)

# Import args
args <- commandArgs()
# print(args)

# Input arguments
workdir<-args[6]
inputfile<-args[7]
projname<-args[8]
distmetric<-args[9]
minSample_inGroup<-args[10]
maxGroup<-args[11]
permutations<-as.numeric(args[12])
numcores<-as.numeric(args[13])
linkage_method<-args[14]

# Import data
feature_raw<-read.csv(inputfile, header=T, row.names=1, sep = ",", check.names = FALSE)

# Set working directory
setwd(workdir)

print("DISSIDE with Adjacency Matrix.")

if (length(rownames(feature_raw))!=length(colnames(feature_raw))) {
  print("Adjacency matrix is not symetrical or is improperly formatted. Exiting")
} else {
  # Count Samples
  samples<-length(rownames(feature_raw))

  # Check for reasonable n
  if (samples/as.numeric(minSample_inGroup)<2) {
    print(paste("Not enough samples to support groups at provided n =", minSample_inGroup, "value. Please lower n."))
  } else {

    # Alternative Clustering (non-euclidean)
    print("Building hclust tree and checking distances.")
    exitcheck=FALSE
    while (!isTRUE(exitcheck)) {
      feature_raw_dist<-vegdist(feature_raw, method=distmetric)
      feature_raw_clust<-hclust(feature_raw_dist, method = linkage_method)
      feature_raw_dist_frame<-as.matrix(feature_raw_dist)
      feature_raw_dist_frame<-as.data.frame(feature_raw_dist_frame, drop=FALSE)
      removeout<-NULL
      for (i in seq(1,length(colnames(feature_raw_dist_frame)), by=1)){
        removeout[i]<-min(feature_raw_dist_frame[,i][feature_raw_dist_frame[,i]!=0])
      }
      removeout<-which(removeout > .9)
      if (length(removeout)>=1){
        print("Removing samples with no reasonable neighbors. Min dist > 0.9")
        print(paste(length(removeout), "samples removed for new Raw data."))
        feature_out<-feature_raw[removeout,-removeout]
        feature_raw<-feature_raw[-removeout,-removeout]
      } else {
        exitcheck=TRUE
      }
    }
    
    # Recount Samples
    samples<-length(rownames(feature_raw))
    
    # Check for sample overlap in features
    if (any(feature_raw_dist_frame==1)){
      print("***WARNING***")
      print("Distance issue. Samples exist in dataset with no feature overlap. Splitting data into individual, coherent data sets.")
      print("***WARNING***")
      datasets<-NULL
      elements<-1
      for (i in seq(1,length(colnames(feature_raw_dist_frame)), by=1)){
        resamp1<-which(feature_raw_dist_frame[i,]<1)
        temp_df<-feature_raw_dist_frame[resamp1,resamp1]
        remsamp2<-as.data.frame(which(temp_df=="1", arr.ind=TRUE))
        remrows<-sort(unique(remsamp2$row))
        if (length(remrows>0)){
          data<-rownames(temp_df[-remrows,])
        } else {
          data<-rownames(temp_df)
        }
        if (!(list(data) %in% datasets)) {
          datasets[[elements]]<-data
          elements<-elements+1
        }
      }
      setwd(paste(workdir, "/Data_Files/", sep = ""))
      for (i in seq(1,length(datasets), by=1)) {
        new_data<-as.data.frame(t(feature_raw[rownames(feature_raw) %in% datasets[[i]],colnames(feature_raw) %in% datasets[[i]]]), drop = FALSE)
        new_data<-new_data[rowSums(new_data)>0,colSums(new_data)>0]
        write.csv(new_data, file = sprintf("%s_cohesive_data_%i.csv", projname, i))
      }
      # Add report for stats on cohesive data sets
      ranks<-NULL
      for (i in seq(1, i, by = 1)) {
        data_file<-read.csv(paste0(projname,"_cohesive_data_",i,".csv"), header=T, row.names=1, check.names=FALSE, sep = ",")
        data_file<-data_file[rowSums(data_file)>0,]
        data_file<-data_file[rowSums(data_file/data_file, na.rm=T)>=2,]
        samples<-length(colnames(data_file))
        otus<-length(rownames(data_file))
        otusample<-NULL
        for (x in seq(1, length(colnames(data_file)), by=1)) {
          otusample[x]<-(length(which(data_file[,x]>0)))
        }
        mincov<-min(otusample)/otus
        maxcov<-max(otusample)/otus
        meancov<-mean(otusample)/otus
        scores<-cbind(i, samples, otus, mincov, maxcov, meancov)
        ranks<-rbind(ranks,scores)
      }
      colnames(ranks)<-cbind("Data", "Samples", "CommonOTUs", "MinCov", "MaxCov", "MeanCov")
      ranks<-as.data.frame(ranks)
      ranks$RankScore<-rank(-apply(cbind(rank(ranks$Samples)*sd(ranks$Samples)/mean(ranks$Samples),
                                         rank(ranks$CommonOTUs)*sd(ranks$CommonOTUs)/mean(ranks$CommonOTUs),
                                         rank(ranks$MeanCov)*sd(ranks$MeanCov)/mean(ranks$MeanCov)),MARGIN=1, FUN="mean"))
      write.csv(ranks, file = "Data_set_stats.csv", row.names = FALSE)
      
      setwd(workdir)
      print(sprintf("%i new data set/s succesfully generated. Please run DISSIDE on individual data sets of interest.", i))
    } else {
      # Save hierarchical cluster plot
      pdf(file = paste("./Figures/Score_Plots/", projname, "_hclust_tree.pdf", sep = ""), width = 10, height =5)
      plot(feature_raw_clust, xlab = "Samples", cex = .5)
      invisible(dev.off())

      # Calculate group differences
      print("Comparing a priori groupings for optimal group count.")
      groups<-2
      groupFs<-NULL
      if (maxGroup=="NULL") {
        # Look at all possible groupings while retaining > 50% samples
        while (sum(table(cutree(feature_raw_clust, k=groups))[which(table(cutree(feature_raw_clust, k=groups))<as.numeric(minSample_inGroup))])<samples/2)
        {
          # Only calculate for groups with at least 2 groups to retain (samples >= n)
          if (length(which(table(cutree(feature_raw_clust, k=groups))>=as.numeric(minSample_inGroup)))>1) {
            feature_raw_clust_cut<-as.data.frame(cutree(feature_raw_clust, k=groups))
            colnames(feature_raw_clust_cut)<-"Group"
            group_adonis<- adonis2(feature_raw_dist ~ as.factor(feature_raw_clust_cut$Group),
                                  permutations = permutations, parallel = numcores)
            group_anosim<- anosim(feature_raw_dist, grouping = as.factor(feature_raw_clust_cut$Group),
                                  permutations = permutations, parallel = numcores)
            # Only retain significant values
            if (group_adonis$`Pr(>F)`[1]<=0.05 && group_anosim$signif<=0.05) {
              groupFs<-rbind(groupFs,cbind(group_adonis$F[1],group_adonis$`Pr(>F)`[1], group_anosim$statistic, group_anosim$signif))
            } else {
              groupFs<-rbind(groupFs,rep(NA,4))
            }
          } else {
            groupFs<-rbind(groupFs,rep(NA,4))
          }
          groups<-groups+1
        }
      } else {
        while (groups<=as.numeric(maxGroup))
        {
          if (length(which(table(cutree(feature_raw_clust, k=groups))>=as.numeric(minSample_inGroup)))>1) {
            feature_raw_clust_cut<-as.data.frame(cutree(feature_raw_clust, k=groups))
            colnames(feature_raw_clust_cut)<-"Group"
            group_adonis<- adonis2(feature_raw_dist ~ as.factor(feature_raw_clust_cut$Group),
                                  permutations = permutations, parallel = numcores)
            group_anosim<- anosim(feature_raw_dist, grouping = as.factor(feature_raw_clust_cut$Group),
                                  permutations = permutations, parallel = numcores)
            if (group_adonis$`Pr(>F)`[1]<=0.05 && group_anosim$signif<=0.05) {
              groupFs<-rbind(groupFs,cbind(group_adonis$F[1],group_adonis$`Pr(>F)`[1], group_anosim$statistic, group_anosim$signif))
            } else {
              groupFs<-rbind(groupFs,rep(NA,4))
            }
          } else {
            groupFs<-rbind(groupFs,rep(NA,4))
          }
          groups<-groups+1
        }
      }
      groups<-groups-1
      groupFs<-as.data.frame(groupFs)
      colnames(groupFs)<-cbind("F","Pp","R","Ap")
      rownames(groupFs)<-2:groups
      
      # Check that n was not too high
      if (length(rownames(na.omit(groupFs)))==0) {
        if (as.numeric(minSample_inGroup)==2){
          print("n=2 already, no a priori disrete patterns are discernable in the data.")
        }else{
          print(paste("A priori groups don't support n =",minSample_inGroup,". Reduce n."))
        }
      } else {
        # Plot Grouping Stats
        pdf(file = paste("./Figures/Score_Plots/", projname, "_optimalgroup_plots.pdf", sep = ""), width = 10, height =5)
        par(mfrow=c(2,2))
        plot(groupFs$F~rownames(groupFs), type="b", main = "PERMANOVA Stat", xlab = "Groups", ylab= "Statistic", cex.axis = .6, xaxt ='n')
        axis(1,rownames(groupFs), cex.axis = .6, labels = rownames(groupFs))
        plot(groupFs$Pp~rownames(groupFs), type="b", main = "PERMANOVA Sig", xlab = "Groups", ylab= "P-value", cex.axis = .6, xaxt ='n')
        axis(1,rownames(groupFs), cex.axis = .6, labels = rownames(groupFs))
        plot(groupFs$R~rownames(groupFs), type="b", main = "ANOSIM Stat", xlab = "Groups", ylab= "Statistic", cex.axis = .6, xaxt ='n')
        axis(1,rownames(groupFs), cex.axis = .6, labels = rownames(groupFs))
        plot(groupFs$Ap~rownames(groupFs), type="b", main = "ANOSIM Sig", xlab = "Groups", ylab= "P-value", cex.axis = .6, xaxt ='n')
        axis(1,rownames(groupFs), cex.axis = .6, labels = rownames(groupFs))
        invisible(dev.off())
        p1<-plot_ly(x =as.numeric(rownames(groupFs)), y=groupFs$F,mode = 'lines+markers', type="scatter",
            marker = list(size = 8, color = 'rgba(255, 182, 193, .9)',line = list(color = 'rgba(152, 0, 0, .8)', width = 2)),
            line = list(color = 'rgba(152, 0, 0, .8)', width = 2),showlegend = FALSE
        ) %>% layout(annotations = list(text="PERMANOVA Stat",xref = "paper", yref = "paper",x = 0.9,y = 1,showarrow = FALSE),xaxis=list(title="Groups"),yaxis=list(title="Statistic"))

        p3<-plot_ly(x =as.numeric(rownames(groupFs)), y=groupFs$R,mode = 'lines+markers', type="scatter",
                             marker = list(size = 8, color = 'rgba(255, 182, 193, .9)',line = list(color = 'rgba(152, 0, 0, .8)', width = 2)),
                             line = list(color = 'rgba(152, 0, 0, .8)', width = 2),showlegend = FALSE
        ) %>% layout(annotations = list(text="ANOSIM Stat",xref = "paper", yref = "paper",x = 0.9,y = 1,showarrow = FALSE),xaxis=list(title="Groups"),yaxis=list(title="Statistic"))

        p2<-plot_ly(x =as.numeric(rownames(groupFs)),y=groupFs$Pp,mode = 'lines+markers', type="scatter",
            marker = list(size = 8, color = 'rgba(30,144,255, .9)',line = list(color = 'rgba(0,0,255, .8)', width = 2)),
            line = list(color = 'rgba(0,0,255, .8)', width = 2),showlegend = FALSE
        ) %>% layout(annotations = list(text="PERMANOVA Sig",xref = "paper", yref = "paper",x = 0.9,y = 1,showarrow = FALSE),xaxis=list(title="Groups"),yaxis=list(title="P-Value"))

        p4<-plot_ly(x =as.numeric(rownames(groupFs)),y=groupFs$Ap,mode = 'lines+markers', type="scatter",
                    marker = list(size = 8, color = 'rgba(30,144,255, .9)',line = list(color = 'rgba(0,0,255, .8)', width = 2)),
                    line = list(color = 'rgba(0,0,255, .8)', width = 2),showlegend = FALSE
        ) %>% layout(annotations = list(text="ANOSIM Sig",xref = "paper", yref = "paper",x = 0.9,y = 1,showarrow = FALSE),xaxis=list(title="Groups"),yaxis=list(title="P-Value"))

        grouping_stats_plotly<-subplot(p1,p2,p3,p4, nrows = 2, margin = 0.05, titleY=TRUE ,titleX=TRUE) %>% layout(margin=list(l=75,r=75,b=75))
        setwd(paste(workdir, "/Figures/Score_Plots/", sep = ""))
        saveWidget(grouping_stats_plotly, paste(projname, "_optimalgroup_plots.html", sep = "") , selfcontained = T)
        setwd(workdir)
      
        # Select groups
        groupFs<-na.omit(groupFs)
        for (i in as.numeric(rownames(groupFs))) {
          # groups no outlier  has weight 1, all sample retained
          if (length(which(table(cutree(feature_raw_clust, k=i))<as.numeric(minSample_inGroup)))==0) {
            groupFs$weight<-1
          } else {
            groupFs$weight[which(rownames(groupFs)==i)]<-
              ((sum(table(cutree(feature_raw_clust, k=i))[which(table(cutree(feature_raw_clust, k=i))>=as.numeric(minSample_inGroup))]))/samples)*
              ((length(which(table(cutree(feature_raw_clust, k=i))<as.numeric(minSample_inGroup))))/length(table(cutree(feature_raw_clust, k=i))))
              ## here length(table(cutree(feature_raw_clust, k=i))) is actually == i

          }
        }
        groupFs$Score<-(((groupFs$F-min(groupFs$F))/(max(groupFs$F)-min(groupFs$F)))*(
          (groupFs$R-min(groupFs$R))/(max(groupFs$R)-min(groupFs$R))))*groupFs$weight
  
        opt_groups<-as.numeric(rownames(groupFs)[which.max(groupFs$Score)])

        pdf(file = paste("./Figures/Score_Plots/", projname,"_scoreplot.pdf",sep=""), width = 10, height = 5)
        plot(groupFs$Score ~ rownames(groupFs), type = "b", main = "Difference Scores", xlab = "Groups", ylab = "Statistic", xaxt="n")
        axis(1,rownames(groupFs), cex.axis = .6, labels = rownames(groupFs))
        invisible(dev.off())  
        score_plotly<-plot_ly(x =as.numeric(rownames(groupFs)),y=groupFs$Score,mode = 'lines+markers', type="scatter",
            marker = list(size = 10, color = 'rgba(255, 182, 193, .9)',line = list(color = 'rgba(152, 0, 0, .8)', width = 2)),
            line = list(color = 'rgba(152, 0, 0, .8)', width = 2)
            ) %>% layout(title="Difference Scores",xaxis=list(title="Groups"),yaxis=list(title="Statistic"),margin=list(l=75,r=75,b=75,t=75))
        setwd(paste(workdir, "/Figures/Score_Plots/", sep = ""))
        saveWidget(score_plotly, paste(projname,"_scoreplot.html",sep="") , selfcontained = T)
        setwd(workdir)
      
        feature_raw_clust_cut<-as.data.frame(cutree(feature_raw_clust, k=opt_groups))
        colnames(feature_raw_clust_cut)<-"Group"
  
        # Force Groups as Factors
        feature_raw_clust_cut$Group<-paste("G",feature_raw_clust_cut$Group, sep = "")
        
        # Remove outlier samples
        print("Removing outlier samples.")
        groupcounts<-table(feature_raw_clust_cut$Group)
        keep<-which(groupcounts>=as.numeric(minSample_inGroup))
        remove<-which(groupcounts<as.numeric(minSample_inGroup))
        feature_no_out<-feature_raw[which((rownames(feature_raw) %in% rownames(feature_raw_clust_cut)) & (feature_raw_clust_cut$Group %in% names(keep))),
                                    which((colnames(feature_raw) %in% colnames(feature_raw_clust_cut)) & (feature_raw_clust_cut$Group %in% names(keep)))]
        
        if (isTRUE(exists("feature_out"))) {
          feature_out<-rbind(feature_out, feature_raw[which((rownames(feature_raw) %in% rownames(feature_raw_clust_cut)) & (feature_raw_clust_cut$Group %in% names(remove))),])
        } else {
          feature_out<-feature_raw[which((rownames(feature_raw) %in% rownames(feature_raw_clust_cut)) & (feature_raw_clust_cut$Group %in% names(remove))),]
        }
        groupskept<-length(names(keep))
        groupsremoved<-length(names(remove))
        feature_clean_clust_cut<-feature_raw_clust_cut[feature_raw_clust_cut$Group %in% names(keep), drop =FALSE,]
  
        write.csv(feature_raw_clust_cut, file = paste("./Data_Files/", projname, "_raw_group_meta.csv", sep = ""))
      
        # Build NMDS for Raw
        print("Building raw NMDS and visualization.")
        sink("/dev/null")
        mds1<-metaMDS(feature_raw_dist, k=3, autotransform=F, na.remove=TRUE)
        sink()
        stress_raw<-mds1$stress

        # Extract points for plotting
        mds1.pts<-mds1$points[, c("MDS1","MDS2","MDS3")]
        mds1.pts<-as.data.frame(mds1.pts, drop=FALSE)
      
        # Ensure data is ordered properly
        feature_raw_clust_cut<-feature_raw_clust_cut[order(rownames(feature_raw_clust_cut), rownames(mds1.pts)),, drop = FALSE]
        mds1.pts<-mds1.pts[order(rownames(mds1.pts), rownames(feature_raw_clust_cut)),]

        # Plot with groups
        getPalette = colorRampPalette(brewer.pal(n=9, name="Set1"))
        colors <- getPalette(length(unique(feature_raw_clust_cut$Group)))
      
        # 2D plots
        views<-c(1,2,3)
        df<-list(mds1.pts[,1:2],mds1.pts[,-2],mds1.pts[2:3])
        myplots = lapply(views, function(col)
          ggplot(as.data.frame(df[[col]]), aes(df[[col]][,1], df[[col]][,2]))+
            geom_point(size=3, aes(color=as.factor(feature_raw_clust_cut$Group)))+
            scale_color_manual(values=colors)+
            guides(color = guide_legend(title="Group"))+
            ggtitle("Raw Data NMDS Optimal Groupings") +
            theme_bw()+
            ylab(colnames(df[[col]])[2])+
            xlab(colnames(df[[col]])[1])+
            theme(legend.key = element_blank(),
                  legend.position = "bottom",
                  legend.direction = "horizontal",
                  legend.box = "horizontal",
                  legend.box.just = "centre",
                  plot.title = element_text(face="bold", size=20)))
      
        plotnames<-c("MDS1xMDS2","MDS1xMDS3","MDS2xMDS3")
        setwd(paste(workdir, "/Figures/2D_NMDS_Plots/", sep = ""))
        for (i in views){
          ggsave(paste(projname, "_2dNMDS_", plotnames[i],"_raw.pdf", sep = ""), myplots[[i]], width = 7,
                 height = 5+.2*(length(names(keep))/5), device = "pdf", units="in")
          options(warn=-1)
          myplots[[i]] <- myplots[[i]] +  theme(plot.margin = margin(2, 4, 2, 2,'cm'))
          NMDS.ggplotly<-ggplotly(myplots[[i]])
          saveWidget(NMDS.ggplotly, paste(projname, "_2dNMDS_", plotnames[i], "_raw.html", sep = ""), selfcontained = T)
          options(warn=0)
        }
      
        # 3D Plots
        setwd(paste(workdir, "/Figures/3D_NMDS_Plots/", sep = ""))
        intfig <- plot_ly(mds1.pts, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~as.factor(feature_raw_clust_cut$Group), colors = colors)
        intfig <- intfig %>% add_markers()
        intfig <- intfig %>% layout(title="Raw NMDS Optimal Groupings", legend=list(title=list(text='Group')))
        saveWidget(intfig, paste(projname, "_3DintNMDS_raw.html", sep = ""), selfcontained = T)
        kaleido(intfig, file = paste(projname, "_3dNMDS_raw.pdf", sep = ""), parallel_limit = numcores, safe = TRUE)
        setwd(workdir)
      
        print("Running analysis on a priori groups with raw data.")
        # Envfit model
        envvars<-envfit(mds1, feature_raw_clust_cut, na.rm=TRUE)
        raw_env_r<-as.numeric(envvars$factors$r)
        raw_env_p<-as.numeric(envvars$factors$pvals)

        # CAP Model
        group_cap<-capscale(feature_raw_dist ~ as.factor(feature_raw_clust_cut$Group), na.action = na.exclude)
        capconst<-sum(group_cap$CCA$eig)/group_cap$tot.chi
        capnova<-anova(group_cap)
        capnova_F<-capnova$F[1]
        capnova_p<-capnova$`Pr(>F)`[1]

        # PERMANOVA
        featureadonis <- adonis2(feature_raw_dist ~ as.factor(feature_raw_clust_cut$Group),
                            permutations = permutations, parallel = numcores)
        raw_perm_F<-featureadonis$F[1]
        raw_perm_p<-featureadonis$`Pr(>F)`[1]

        #ANOSIM
        featureanosim <- anosim(feature_raw_dist, grouping = as.factor(feature_raw_clust_cut$Group),
                            permutations = permutations, parallel = numcores)
        raw_anosim_R<-featureanosim$statistic
        raw_anosim_p<-featureanosim$signif

        # Store Raw Analysis Results
        header<-c("NMDS_Stress(k=3)","Samples","Groups","ENVFIT_r2","ENVFIT_p","CAP_constrained","CAP_ANOVA_F",
                  "CAP_ANOVA_p","PERMANOVA_F","PERMANOVA_p","ANOSIM_R","ANOSIM_p")
        raw_analysis<-c(stress_raw,samples,opt_groups,raw_env_r,raw_env_p,capconst,capnova_F,capnova_p,raw_perm_F,
                        raw_perm_p,raw_anosim_R,raw_anosim_p)
        
        # Calculate raw centroids
        #Find the centroids of the clusters: raw -> feature_raw_dist, group_num = opt_groups
        #cluster the distances based on optimal group numbers which were previously calculated
        clust_raw <- cutree(feature_raw_clust, k=opt_groups)
        clust_raw_distMatrix <-as.matrix(feature_raw_dist)
        clustered_items_raw <- split(names(clust_raw), clust_raw)
        
        #create a table for collecting the centroids names for each cluster
        centroid_df_raw <- data.frame(matrix(ncol = opt_groups, nrow = 1))
        df_labels_raw <- paste("G", labels(clustered_items_raw), sep="")  
        colnames(centroid_df_raw) <- df_labels_raw
        rownames(centroid_df_raw) <- c("Centroid")
        
        for (group in 1:opt_groups) {
          #the current cluster will be the samples in the current group # of the for loop
          curr_cluster_raw <- clustered_items_raw[[group]]
          
          #if there is only item in the cluster, it is automatically the centroid
          if (length(curr_cluster_raw) <= 1) {
            centroid_df_raw[group] <- curr_cluster_raw
          } else { 
            #centroid will have minimum distance to all of the items in the cluster
            sums <- rowSums(clust_raw_distMatrix[curr_cluster_raw, curr_cluster_raw])
            centroid_df_raw[group] <- names(which.min(sums))
          }
        }
        centroid_t_raw <- t(centroid_df_raw)
        write.csv(centroid_t_raw, file=paste("./Data_Files/", projname, "_centroids_raw.csv", sep = ""))
    
        if (length(rownames(feature_no_out))==length(rownames(feature_raw))) {
          print("No additional samples removed as outliers. Performing analysis on full data set only.")
          clean_analysis<-raw_analysis
          diffcalc<-clean_analysis-raw_analysis
          anoutfile<-rbind(raw_analysis,clean_analysis,diffcalc)
          colnames(anoutfile)<-header
          rownames(anoutfile)<-c("Raw","Cleaned","Difference")
          write.csv(anoutfile, file = paste(projname, "_analysis.csv", sep = ""))
        } else {
          feature_no_out_dist<-vegdist(feature_no_out, method=distmetric)
          # Cleaned NMDS
          print(paste(length(rownames(feature_out)),"samples removed as group outliers."))
          print("Building cleaned NMDS and visualization.")
          sink("/dev/null")
          mds2<-metaMDS(feature_no_out_dist, k=3, autotransform=F, na.remove=TRUE)
          sink()

          stress_clean<-mds2$stress

          # Extract points for plotting
          mds2.pts<-mds2$points[, c("MDS1","MDS2","MDS3")]
          mds2.pts<-as.data.frame(mds2.pts, drop=FALSE)
        
          # Ensure data is ordered properly
          feature_clean_clust_cut<-feature_clean_clust_cut[order(rownames(feature_clean_clust_cut), rownames(mds2.pts)),, drop = FALSE]
          mds2.pts<-mds2.pts[order(rownames(mds2.pts), rownames(feature_clean_clust_cut)),]
        
          # Plot with groups
          getPalette = colorRampPalette(brewer.pal(n=9, name="Set1"))
          colors <- getPalette(length(unique(feature_clean_clust_cut$Group)))
        
          #2D Plots
          views<-c(1,2,3)
          df<-list(mds2.pts[,1:2],mds2.pts[,-2],mds2.pts[2:3])
          myplots = lapply(views, function(col)
            ggplot(as.data.frame(df[[col]]), aes(df[[col]][,1], df[[col]][,2]))+
              geom_point(size=3, aes(color=as.factor(feature_clean_clust_cut$Group)))+
              scale_color_manual(values=colors)+
              guides(color = guide_legend(title="Group"))+
              ggtitle("Cleaned Data NMDS Optimal Groupings") +
              theme_bw()+
              ylab(colnames(df[[col]])[2])+
              xlab(colnames(df[[col]])[1])+
              theme(legend.key = element_blank(),
                    legend.position = "bottom", 
                    legend.direction = "horizontal",
                    legend.box = "horizontal",
                    legend.box.just = "centre",
                    plot.title = element_text(face="bold", size=20)))
        
          plotnames<-c("MDS1xMDS2","MDS1xMDS3","MDS2xMDS3")
          setwd(paste(workdir, "/Figures/2D_NMDS_Plots/", sep = ""))
          for (i in views){
            ggsave(paste(projname, "_2dNMDS_", plotnames[i],"_clean.pdf", sep = ""), myplots[[i]], width = 7,
                   height = 5+.2*(length(names(keep))/5), device = "pdf", units="in")
            options(warn=-1)
            myplots[[i]] <- myplots[[i]] +  theme(plot.margin = margin(2, 4, 2, 2,'cm'))
            NMDS.ggplotly<-ggplotly(myplots[[i]])
            saveWidget(NMDS.ggplotly, paste(projname, "_2dNMDS_", plotnames[i], "_clean.html", sep = ""), selfcontained = T)
            options(warn=0)
          }
        
          # 3D Plots
          setwd(paste(workdir, "/Figures/3D_NMDS_Plots/", sep = ""))
          intfig <- plot_ly(mds2.pts, x = ~MDS1, y = ~MDS2, z = ~MDS3, color = ~as.factor(feature_clean_clust_cut$Group), colors = colors)
          intfig <- intfig %>% add_markers()
          intfig <- intfig %>% layout(title="Cleaned NMDS Optimal Groupings", legend=list(title=list(text='Group')))
          saveWidget(intfig, paste(projname, "_3DintNMDS_clean.html", sep = "") , selfcontained = T)
          kaleido(intfig, file = paste(projname, "_3dNMDS_clean.pdf", sep = ""), parallel_limit = numcores, safe = TRUE)
          setwd(workdir)
        
          print("Running analysis on a priori groups with cleaned data.")
          # Envfit model
          envvars2<-envfit(mds2, feature_clean_clust_cut, na.rm=TRUE)
          clean_env_r<-as.numeric(envvars2$factors$r)
          clean_env_p<-as.numeric(envvars2$factors$pvals)

          # CAP Model
          group_cap_clean<-capscale(feature_no_out_dist ~ as.factor(feature_clean_clust_cut$Group), na.action = na.exclude)
          capconst_clean<-sum(group_cap_clean$CCA$eig)/group_cap_clean$tot.chi
          capnova_clean<-anova(group_cap_clean)
          capnova_F_clean<-capnova_clean$F[1]
          capnova_p_clean<-capnova_clean$`Pr(>F)`[1]

          # PERMANOVA
          featureadonis_c<-adonis2(feature_no_out_dist ~ as.factor(feature_clean_clust_cut$Group),
                              permutations = permutations, parallel = numcores)
          clean_perm_F<-featureadonis$F[1]
          clean_perm_p<-featureadonis$`Pr(>F)`[1]

          #ANOSIM
          featureanosim_c<-anosim(feature_no_out_dist, grouping = as.factor(feature_clean_clust_cut$Group),
                              permutations = permutations, parallel = numcores)
          clean_anosim_R<-featureanosim_c$statistic
          clean_anosim_p<-featureanosim_c$signif
      
          # Build output files
          print("Formatting and writing output files.")
          # Analysis details
          clean_analysis<-c(stress_clean,length(rownames(feature_no_out)),length(names(keep)),clean_env_r,clean_env_p,
                            capconst_clean,capnova_F_clean,capnova_p_clean,clean_perm_F,clean_perm_p,clean_anosim_R,clean_anosim_p)
          diffcalc<-clean_analysis-raw_analysis

          anoutfile<-rbind(raw_analysis,clean_analysis,diffcalc)
          colnames(anoutfile)<-header
          rownames(anoutfile)<-c("Raw","Cleaned","Difference")

          # Cleaned Samples
          feature_no_out<-t(feature_no_out)
          
          # Calculate clean centroids
          clust_clean <- clust_raw[][which(names(clust_raw) %in% rownames(feature_clean_clust_cut))]
          clust_clean_distMatrix <-as.matrix(feature_no_out_dist)
          clustered_items_clean <- split(names(clust_clean), clust_clean)
          
          #create a table for collecting the centroids names for each cluster
          centroid_df_clean <- data.frame(matrix(ncol = groupskept, nrow = 1))
          df_labels_clean <- paste("G", labels(clustered_items_clean), sep="")  
          colnames(centroid_df_clean) <- df_labels_clean
          rownames(centroid_df_clean) <- c("Centroid")
          
          for (group in 1:groupskept) {
            #the current cluster will be the samples that were placed into the current group # of the for loop
            curr_cluster_clean <- clustered_items_clean[[group]]
            
            #if there is only item in the cluster, it is automatically the centroid
            if (length(curr_cluster_clean) <= 1) {
              centroid_df_clean[group] <- curr_cluster_clean
            } else { 
              #centroid will have minimum distance to all of the items in the cluster
              sums <- rowSums(clust_clean_distMatrix[curr_cluster_clean, curr_cluster_clean])
              centroid_df_clean[group] <- names(which.min(sums))
            }
          }
          centroid_t_clean <- t(centroid_df_clean)
          write.csv(centroid_t_clean, file=paste("./Data_Files/", projname, "_centroids_cleaned.csv", sep = ""))
     
          # Write files out
          write.csv(anoutfile, file = paste(projname, "_analysis.csv", sep = ""))
          write.csv(feature_no_out, file = paste("./Data_Files/", projname, "_cleaned_sample_feature.csv", sep = ""))
          write.csv(feature_raw_clust_cut, file = paste("./Data_Files/", projname, "_raw_group_meta.csv", sep = ""))
          write.csv(feature_clean_clust_cut, file = paste("./Data_Files/", projname, "_cleaned_group_meta.csv", sep = ""))
        }
      }
      # Output raw data used
      feature_raw<t(feature_raw)
      write.csv(feature_raw, file = paste("./Data_Files/", projname, "_raw_sample_feature.csv", sep = ""))
    
      # Output removed samples file if exist
      if (isTRUE(exists("feature_out"))){
        feature_out<-t(feature_out)
        write.csv(feature_out, file = paste("./Data_Files/", projname, "_removed_sample_feature.csv", sep = ""))
      }
    }
  }
}