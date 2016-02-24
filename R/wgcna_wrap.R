#' soft-thresholding powers picking
#'
#' As suggested by WGCNA, the 1st number that hists
#' between 0.8-0.9 region should be picked as the soft-threshold
#'
#' @param x - adjacency matrix returned by phaseLockingMatrix, grangerTest or coherenceTest
#' @return  pdf file with R^2 plot and average connectivity
#' @examples
#' #example 1
#' #be aware that this example data is generating a random matirx, so the plot does not match the real situation.  
#' x = replicate(1000, qunif(sort(runif(1000,min=0,max=1)), min = 0, max = 1))
#' x[lower.tri(x)] = t(x)[lower.tri(x)]
#' wgcnaPower(x)
#'
#' @references
#' [1] Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008)
#'
#' @export
wgcnaPower <- function(x) {
  # The following setting is important, do not omit.
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
  
  if (min(x)<0){
    networkType="signed"
  }else {networkType="unsigned"}
  
  # Call the network topology analysis function
  sft = WGCNA::pickSoftThreshold(x, dataIsExpr = FALSE, networkType = networkType,powerVector = powers, verbose = 5)
  
  return(sft)
  # Plot the results:
  #sizeGrWindow(9, 5)

}


#' WGCNA modules
#'
#' As suggested by WGCNA, the 1st number that hists
#' between 0.8-0.9 region should be picked as the soft-threshold
#'
#' @param x - adjacency matrix returned by phaseLockingMatrix, grangerTest or coherenceTest
#' @param softPower - based on the plot of function wgcna_analysis_power, the 1st number that hists between 0.8-0.9 region should be picked as the soft-threshold. If none, use 1.
#' @return  3 plots.
#' @examples
#' #example 1
#' x = replicate(1000, rnorm(1000))
#' wgcna_result=wgcna_analysis(x,5)
#'
#' @references
#' [1] Langfelder, P. & Horvath, S. WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 9, 559 (2008)
#'
#' @export
wgcnaAnalysis <- function(x, softPower=1) {
  adjacency = WGCNA::adjacency(x, power = softPower);
  # Turn adjacency into topological overlap
  TOM = WGCNA::TOMsimilarity(adjacency);
  dissTOM = 1 - TOM
  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = "average");
  
  # par(mfrow = c(1,3));
#   # Plot the resulting clustering tree (dendrogram)
#   sizeGrWindow(12,9)
#   
#   plot(
#     geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based dissimilarity",labels = FALSE, hang = 0.04
#   );
  
  # We like large modules, so we set the minimum module size relatively high:
  minModuleSize = 30;
  # Module identification using dynamic tree cut:
  dynamicMods = dynamicTreeCut::cutreeDynamic(
    dendro = geneTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE,minClusterSize = minModuleSize
  );
  table(dynamicMods)
  # Convert numeric lables into colors
  dynamicColors = WGCNA::labels2colors(dynamicMods)
  table(dynamicColors)
  # Plot the dendrogram and colors underneath
  sizeGrWindow(8,6)
#   WGCNA::plotDendroAndColors(
#     geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors"
#   )
  
  # Calculate eigengenes
  MEList = WGCNA::moduleEigengenes(x, colors = dynamicColors)
  MEs = MEList$eigengenes
  # Calculate dissimilarity of module eigengenes
  MEDiss = 1 - cor(MEs);
  # Cluster module eigengenes
  METree = hclust(as.dist(MEDiss), method = "average");
  MEDissThres = 0.25
  # Call an automatic merging function
  merge = WGCNA::mergeCloseModules(x, dynamicColors, cutHeight = MEDissThres, verbose = 3)
  # The merged module colors
  mergedColors = merge$colors;
  # Eigengenes of the new merged modules:
  mergedMEs = merge$newMEs;
  # Rename to moduleColors
  moduleColors = mergedColors
  # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
  plotTOM = dissTOM ^ 7;
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA;
#   # Call the plot function
#   WGCNA::sizeGrWindow(9,9)
  #plot here
  # WGCNA::TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
  names(dynamicMods)=rownames(x)
  names(dynamicColors)=rownames(x)

  return(list("modules"=dynamicMods,"geneTree"=geneTree,"modulColors"=dynamicColors))
}


#' Plot modules
#'
#' To visualize modules in graph and in heatmap
#'
#' @param x - co-sychronization matrix returned by phaseLockingMatrix, grangerTest or coherenceTest
#' @param expr - expression profile
#' @param modulColors - color vector return by wgcnaAnalysis
#' @param thisColor - the color of the module to be plotted
#' @param cut - the threshold for the co-sychronization matrix, extacting submatrix that is bigger than this cutoff, and plot the graph.
#' 
#' @return  A graph and a heatmap
#' @examples
#'
#' @export
ModulePlots <-
  function(x,expr,mColors,thisColor = "blue",cut = 0.3) {
    #bug here
    #row names != col names of x
    #colnames(x)=rownames(x)

    
    #extract the adjacency matrix in the same module
    this_adj = x[names(mColors[mColors == thisColor]),
                 names(mColors[mColors == thisColor])]
    this_adj[this_adj<cut]=0
    diag(this_adj)=0
    
    #plot the vertices with edges incident to 
    all_names=unique(unlist(apply(this_adj,1,function(x)names(which(x>0)))))
    plot_adj=this_adj[all_names,all_names]
    
    net = igraph::graph.adjacency(
      as.matrix(plot_adj),mode = "undirected",weighted = TRUE,diag = FALSE
    )
    
    V(net)$size<-igraph::degree(net)
    
    #for phase gamma and coherence, we can set the edge weight this way. Not for granger p-value
    igraph::plot.igraph(
      net,vertex.label = V(net)$name,layout = layout.fruchterman.reingold,
      edge.color = "black",vertex.label.color='black',vertex.label.font=1,
      vertex.label.cex=0.5,edge.width = E(net)$weight,main = paste("Module", thisColor)
    )

    
    #heatmap, require library(gplot)
    this_expr = expr[names(mColors[mColors == thisColor]),]
    heatmap(this_expr,Rowv = TRUE, Colv=FALSE,labRow=FALSE,col=cm.colors(256),main = paste("Module", thisColor))
    
  }

#' Output the modules with GO annotation
#'
#' To export modules with GO.
#'
#' @param colors - numeric color vector return by wgcnaAnalysis.
#' @return  Table with the module color, gene names, GO annotation and statistics etc.
#' @examples
#'
#' @export
ModuleAnnotation <-
  function(colors) {
    geneList = names(colors)
    #Generate module annotation table
    #Step 3: take the results and output them to an excel file
    #warning:this uses GO.db and AnnotationDBI and the organism specific annotation
    #the inputs are colorlist and genename
    #colorlist <- lnsing$colors
    
    #the actual code converts the hgnc_symbol to entrez id, which is needed to use the GOenrichment analysis
    ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mapTab <-
      biomaRt::getBM(
        attributes = c("hgnc_symbol", "entrezgene"), filters = "hgnc_symbol", values = geneList, mart = ensembl, uniqueRows =
          FALSE
      )
    #the next three lines remove duplicate rows
    dupRows <-
      union(which(duplicated(mapTab[,1])), which(duplicated(mapTab[,2])))
    entrezIds <- mapTab[-dupRows, 2]
    names(entrezIds) <- mapTab[-dupRows, 1]
    #this matches the resulting entrezIds to their respective genes
    list <-
      lapply(geneList, function(x)
        if (length(which(x == names(entrezIds))) > 0)
          entrezIds[[which(x == names(entrezIds))]])
    #remove any null values and replace with 0
    list1 <- lapply(list, function(x)
      if (is.null(x))
        0
      else
        x);
    list2 = unlist(list1)
    moduleColors = WGCNA::labels2colors(colors)
    #non matched are returned as 0
    #this is the enrichment analysis, which only takes entrez id
    #it takes the 10 best fit, generally this enrichment is poorer than online ones and should only
    #be used as a preliminary test, not for publications
    GOenr = WGCNA::GOenrichmentAnalysis(moduleColors, list2, nBestP = 10);
    tab = GOenr$bestPTerms[[4]]$enrichment
    return(tab)
  }

#' Output the 2 modules overlaps.
#'
#' To export the number of genes shared between modules.
#'
#' @param modules1 - modules return by wgcnaAnalysis.
#' @param modules2 - modules return by wgcnaAnalysis.
#' @return  A figure and table with the number of overlapping genes, 
#'  correspondence (correspondence= number of overlap/ avergae module size),
#'  size of modules.
#' @examples
#'
#' @export
wgcnaOverlap <- function(modules1,modules2) {
  
  #   #begin the cross comparison
  #   # Isolate the module labels in the order they appear in ordered module eigengenes
  module1Colors=WGCNA::labels2colors(modules1)
  module2Colors=WGCNA::labels2colors(modules2)
  
  # Convert the numeric module labels to color labels
  color_set1= WGCNA::labels2colors(as.numeric(unique(modules1)))
  color_set2 = WGCNA::labels2colors(as.numeric(unique(modules2)))
  # Number of modules
  n_Mod1 = length(color_set1)
  n_Mod2 = length(color_set2)

    # Initialize tables of correspondence and of the corresponding counts
  correspd_table  = matrix(0, nrow = n_Mod1, ncol = n_Mod2);
  count_table = matrix(0, nrow = n_Mod1, ncol = n_Mod2);
  # Execute all pairwaise comparisons
  for (i in 1:n_Mod1)
    for (j in 1:n_Mod2)
    {
      members1 = (module1Colors == color_set1[i]);
      members2 = (module2Colors == color_set2[j]);
      if (sum(members1) > 400 || sum(members2) > 400) {
        correspd_table [i,j] = 0;
      }else{
        correspd_table [i,j] = sum(members1 + members2 == 2) / (0.5 * (sum(members1) +
                                                                        sum(members2)))
      };
      #pTable[i, j] = -log10(fisher.test(members1, members2, alternative = "greater")$p.value);
      count_table[i, j] = sum(module1Colors == color_set1[i] &
                            module2Colors ==
                            color_set2[j])
    }
  
  
  # Marginal counts (really module sizes)
  mod1_totals = apply(count_table, 1, sum)
  mod2_totals = apply(count_table, 2, sum)
  # Actual plotting
  sizeGrWindow(10,7 );
  par(mfrow=c(1,1));
  par(cex = 1.0);
  par(mar=c(8, 10.4, 2.7, 1)+0.3);
  # Use function labeledHeatmap to produce the color-coded table with all the trimmings
  WGCNA::labeledHeatmap(Matrix = correspd_table,
                 xLabels = paste("mod1", color_set1),
                 yLabels = paste("mod2", color_set2),
                 colorLabels = TRUE,
                 xSymbols = paste(color_set1, mod1_totals, sep=":"),
                 ySymbols = paste(color_set2, mod2_totals, sep=":"),
                 textMatrix = count_table,
                 colors = blueWhiteRed(100)[50:100],
                 main = "overlapping matrix colored by correspondence",
                 cex.text = 1.0, cex.lab = 1.0, setStdMargins = FALSE);
  
  return(list(correspondence=correspd_table,counts=count_table,
              mod1_totals=mod1_totals,mod2_totals=mod2_totals))
}

#' Output the 2 modules overlaps.
#'
#' To export the number of genes shared between modules.
#'
#' @param colors - numeric color vector return by wgcnaAnalysis.
#' @param selectColor - celcet the color of the module to perform Go enrichment analysis
#' @param identifier - default is "symbol". One in  c("Entrez", GenBank, Alias, Ensembl, GeneSymbol", "GeneName" and "UniGene"). Refer to topGO. 
#' @param species - "Human" or "Mouse"
#' @param pvalue - export the GO terms with p-value smaller than this number
#' @param fileTree - TRUE if you wish to export the GO tree
#' @param numSigTerm - a parameter past to topGO function printGraph, determines how many top GO terms will be plotted in rectangle.
#' @return  GO enrichment result, and the tree structure GO plot with significant GO terms highlighted in red and orange colors. The pdf file of the tree will be saved in the working directory.
#' @examples
#'
#' @export
topGOanalysis <- function(colors, selectColor="blue", identifier="symbol", species="Human", pvalue=0.05, fileTree=TRUE, numSigTerm=5) {
  geneColors=WGCNA::labels2colors(colors)

  if (species == "Human") {
    mapdb <- "org.Hs.eg.db"      
  } else if (species == "Mouse") {
    mapdb <- "org.Mm.eg.db"      
  }         
  res <- list()
  
  allgene <- names(colors)
  inputgene <- allgene[geneColors==selectColor]
  ##unique and remove 0
  allgene_unique<-setdiff(unique(allgene), 0)
  inputgene_unique<-setdiff(unique(inputgene), 0)
  
  geneList <- factor(as.integer(allgene_unique %in% inputgene_unique))
  names(geneList) <- allgene_unique
  
  #prepare the data
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.org, mapping = mapdb, ID = identifier)
  
  resultFisher <- topGO::runTest(GOdata, algorithm = "classic", statistic = "fisher")
  goscore<-sort(score(resultFisher))
  

  if(fileTree){
    #if (!exists(numSigTerm)){
    #  numSigTerm<-5
    #}

    topGO::showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = numSigTerm, useInfo = 'all')
  
    topGO::printGraph(GOdata, resultFisher, firstSigNodes = numSigTerm, fn.prefix = fileTree, useInfo = "all", pdfSW = TRUE)

    }

  goscore=goscore[goscore<=pvalue]
  res= matrix( rep(0, length(goscore)*3), ncol=3) 
  colnames(res)<-c("GOterm", "GODesc", "pValue")
  for (i in 1:length(goscore)){
    res[i,1]=names(goscore)[i][1]
    res[i,2]=Term(names(goscore)[i][1])
    res[i,3]=goscore[i]
  }
  return(res)
}


