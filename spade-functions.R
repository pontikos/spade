library(spade)

### The code has been modified to work with .RData files instead of .FCS files


###
applyTransforms <- function(x, transforms) {
    n <- 1
    transforms <- transforms[intersect(names(transforms),colnames(x))]
    apply(x, 2, function(y) {
             if (n > length(transforms)) return(y)
             y <- transforms[[colnames(x)[n]]](y)
             n <<- n+1
             return(y)
            })
}

###
SPADE.markerMedians2 <- function (files, num.clusters, cols = NULL, transforms = flowCore::arcsinhTransform(a = 0, b = 0.2), cluster_cols = NULL) {
    data <- c()
    print(files)
    for (f in files) {
        print(load(f))
        in_data <- fcs.data
        data <- rbind(data, in_data)
    }
    clst <- data[, "cluster"]
    print(head(data <- data[, colnames(data) != "cluster", drop = FALSE]))
    data <- applyTransforms(data, transforms)
    colnames(data) <- sapply(colnames(data), function(x) {
            if (x %in% cluster_cols) x <- paste(x, "clust", sep = "_")
            return(x)
    })
    ids <- 1:num.clusters
    if (any(is.na(match(unique(clst), ids)))) {
        stop("More clusters in FCS files than indicated")
    }
    count <- matrix(0, nrow = num.clusters, ncol = 1, dimnames = list(ids, "count"))
    medians <- matrix(NA, nrow = num.clusters, ncol = ncol(data), dimnames = list(ids, colnames(data)))
    cvs <- matrix(NA, nrow = num.clusters, ncol = ncol(data), dimnames = list(ids, colnames(data)))
    for (i in ids) {
        data_s <- subset(data, clst == i)
        count[i, 1] <- nrow(data_s)
        medians[i, ] <- apply(data_s, 2, median)
        cvs[i, ] <- apply(data_s, 2, function(d) 100 * sd(d)/abs(mean(d)) )
    }
    percenttotal <- matrix((count/sum(count)) * 100, nrow = num.clusters, ncol = 1, dimnames = list(ids, "percenttotal"))
    list(count = count, medians = medians, cvs = cvs, percenttotal = percenttotal)
}

###

###

### fcs.data
SPADE.cluster <- function(tbl, k) {
  if (nrow(tbl) > 60000) warning("Potentially too many observations for the clustering step",immediate=TRUE);
  # Transpose table before call into row major order
  clust <- .Call("SPADE_cluster",t(tbl),as.integer(k))
  # Invalid clusters have assgn == 0
  centers = c()
  is.na(clust$assgn) <- which(clust$assgn == 0)
  for (i in c(1:max(clust$assgn, na.rm=TRUE))) {  
       obs <- which(clust$assgn == i)
       if (length(obs) > 1) {
            centers <- rbind(centers,colMeans(tbl[obs,,drop=FALSE]))
            clust$assgn[obs] <- nrow(centers)
       } else {
            is.na(clust$assgn) <- obs
      }
    }
    return(list(centers=centers,assign=clust$assgn))
}

###
### fcs.data contains density column
### after this function, the boolean keep column
### determines whether to keep the data in the downsampling
SPADE.downsample <- function(fcs.data, exclude_pctile=0.01, target_pctile=0.05, desired_samples=NULL) {
    boundary <- quantile(fcs.data[, 'density'], c(exclude_pctile, target_pctile), names = FALSE)
    #out_data <- subset(fcs.data, fcs.data[,'density'] > boundary[1])
    fcs.data[,'keep'] <- FALSE
    #first step: include those which have a density higher than the lowest pct
    fcs.data[,'keep'] <- fcs.data[,'density'] > boundary[[1]]
    out_data <- fcs.data[fcs.data[,'keep'],]
    density <- out_data[,'density']
    if (is.null(desired_samples)) {
        boundary <- boundary[2]
        #out_data <- subset(out_data, boundary/density>runif(nrow(out_data)))
        fcs.data[fcs.data[,'keep'],][which(boundary/density<=runif(nrow(out_data))),'keep'] <- FALSE
    } else if (desired_samples < nrow(out_data)) {
        density_s <- sort(density)
        cdf <- rev(cumsum(1/rev(density_s)))
        boundary <- desired_samples/cdf[1]
        if (boundary > density_s[1]) {
            targets <- (desired_samples - 1:length(density_s))/cdf
            boundary <- targets[which.min(targets - density_s > 0)]
        }
        #out_data <- subset(out_data, boundary/density > runif(length(density)))
        fcs.data[fcs.data[,'keep'],][which(boundary/density<=runif(length(density))),'keep'] <- FALSE
    }
    print(prop.table(table(fcs.data[,'keep'])))
    return(fcs.data)
}

###
SPADE.plot.trees2 <- function (graph, individual=NULL, day=NULL, layout = SPADE.layout.arch, attr_pattern = "percent|medians|fold|cvs", scale = NULL, pctile_color=c(0.02, 0.98), normalize = "global", size_scale_factor=2, edge.color = "grey", bare = FALSE, palette = "bluered", DOSES=c('0U','01U','10U','1000U')) {
    print('mine')
    if (!is.igraph(graph)) stop("Not a graph object")
    if (!is.null(scale) && (!is.vector(scale) || length(scale) != 2)) stop("scale must be a two element vector")
    if (!is.vector(pctile_color) || length(pctile_color) != 2) stop("pctile_color must be a two element vector with values in [0,1]")
    if (length(files) == 1 && file.info(SPADE.strip.sep(files))$isdir) files <- dir(SPADE.strip.sep(files), full.names = TRUE, pattern = glob2rx(file_pattern))
    load_attr <- function(save_file) {
        anno <- NULL
        l <- load(save_file)
        stopifnot(l == "anno")
        return(anno)
    }
    boundaries <- NULL
    if (normalize == "global") {
        boundaries <- c()
        all_attrs <- c()
        for (f in files) {
            attrs <- load_attr(f)
            for (i in grep(attr_pattern, colnames(attrs))) {
                n <- colnames(attrs)[i]
                all_attrs[[n]] <- c(all_attrs[[n]], attrs[, i])
            }
        }
        for (i in seq_along(all_attrs)) {
            boundaries[[names(all_attrs)[i]]] <- quantile(all_attrs[[i]], probs = pctile_color, na.rm = TRUE)
        }
    }
    if (is.function(layout)) graph_l <- layout(graph)
    else graph_l <- layout
    if (palette == "jet") 
        palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    else if (palette == "bluered") 
        palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
    else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
    colorscale <- palette(100)
    print(files)
    #par(mfrow=c(2,2))
    #for (dose in c('0U','01U','10U','1000U')) {
    for (dose in DOSES) {
        #cluster file
        #print(f <- file.path( sprintf('%s.fcs.density.fcs.cluster.fcs',paste(individual, dose, day, sep='_')) ))
        #print(length(unique(clusters <- (flowCore::read.FCS(f))@exprs[,'cluster'])))
        print(f <- file.path( sprintf('%s.density.cluster',paste(individual, dose, day, sep='_')) ))
        print(load(f))
        print(length(unique(clusters <- fcs.data[,'cluster'])))
        #anno file
        f <- file.path( sprintf('%sdensity.cluster.anno.Rsave',paste(individual, dose, day, sep='_')) )
        attrs <- load_attr(f)
        vsize <- attrs$percenttotal
        vsize <- vsize/(max(vsize, na.rm = TRUE)^(1/size_scale_factor)) * 3 + 2
        #vsize <- vsize/(max(vsize, na.rm = TRUE)^(1/size_scale_factor))*3+1
        vsize[is.na(vsize) | (attrs$count == 0)] <- 1
        print(grep(attr_pattern, colnames(attrs), value=TRUE))
        #pdf(pdf.file)
        for (i in grep(attr_pattern, colnames(attrs))) {
            print(name <- colnames(attrs)[i])
            attr <- attrs[, i]
            if (!is.null(scale)) 
                boundary <- scale
            else if (normalize == "global") {
                boundary <- boundaries[[name]]
            }
            else boundary <- quantile(attr, probs = pctile_color, na.rm = TRUE)
            if (length(grep("^medians|percent|cvs", name))) boundary <- c(min(boundary), max(boundary))
            else boundary <- c(-max(abs(boundary)), max(abs(boundary)))
            boundary <- round(boundary, 2)
            if (boundary[1] == boundary[2]) {
                boundary <- c(boundary[1] - 1, boundary[2] + 1)
            }
            cat(boundary[1], boundary[2], '\n')
            grad <- seq(boundary[1], boundary[2], length.out = length(colorscale))
            color <- colorscale[findInterval(attr, grad, all.inside = TRUE)]
            color[is.na(attr) | (attrs$count == 0)] <- "grey"
            #if cluster in manual gate then color is gray
            #color[clusters[which(rowSums(CLR)>0)]] <- 'grey'
            #if (grepl("^percenttotalratiolog$", name)) color[is.na(attr) & attrs$count > 0] <- tail(colorscale, 1)
            fill_color <- color
            is.na(fill_color) <- is.na(attr)
            frame_color <- color
            #graph_aspect <- ((max(graph_l[, 2]) - min(graph_l[, 2]))/(max(graph_l[, 1]) - min(graph_l[, 1])))
            plot(graph, layout=graph_l, vertex.shape="circle", vertex.color=fill_color, vertex.frame.color=frame_color, edge.color=edge.color, vertex.size=vsize, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1) #, asp=graph_aspect)
            #plot(layout, pch=20, cex = vsize, col=color, yaxt='n', xaxt='n', ann=FALSE,frame.plot=FALSE,main=dose)
            #lines(G$x,G$y,lwd=2,col='purple')
            if (!bare) {
                if (length(grep("^medians", name))) name <- sub("medians", "Median of ", name)
                else if (length(grep("^fold", name))) name <- sub("fold", "Arcsinh diff. of ", name)
                else if (grepl("^percenttotal$", name)) name <- sub("percent", "Percent freq. of ", name)
                else if (grepl("^percenttotalratiolog$", name)) name <- "Log10 of Ratio of Percent Total of Cells in Each Cluster"
                else if (grepl("^cvs", name)) name <- sub("cvs", "Coeff. of Variation of ", name)
                if (grepl("_clust$", name)) name <- sub("_clust", "\n(Used for tree-building)", name)
                #title(main = paste(strsplit(basename(f), ".fcs")[[1]][1], sub = name, sep = "\n"))
                #legend
            }
        }
        #dev.off()
    }
    #subplot(image(grad, c(1), matrix(1:length(colorscale), ncol = 1), col = colorscale, xlab = ifelse(is.null(scale), paste("Range:", pctile_color[1], "to", pctile_color[2], "pctile"), ""), ylab = "", yaxt = "n", xaxp = c(boundary, 1)), x = "right,bottom", size = c(1, 0.2))
    return(data.frame(vsize=vsize, color=color))
}
tmpfun <- get("SPADE.plot.trees", envir = asNamespace("spade"))
environment(SPADE.plot.trees2) <- environment(tmpfun)
attributes(SPADE.plot.trees2) <- attributes(tmpfun)  # don't know if this is really needed
assignInNamespace('SPADE.plot.trees', SPADE.plot.trees2, ns='spade')
unlockBinding('SPADE.plot.trees', as.environment('package:spade'));
assignInNamespace('SPADE.plot.trees', SPADE.plot.trees2, envir=as.environment('package:spade'))
assign('SPADE.plot.trees', SPADE.plot.trees2, 'package:spade')
lockBinding('SPADE.plot.trees', as.environment('package:spade'))


###
SPADE.layout.arch2 <- function (mst_graph) {
    if (!is.igraph(mst_graph)) {
        stop("Input has to be igraph object")
    }
    if (!is.connected(mst_graph)) {
        stop("Cannot handle graph that has disjoint components")
    }
    if (girth(mst_graph)$girth > 0) {
        stop("Cannot handle graphs with cycles")
    }
    hops <- c()
    for (v in V(mst_graph)) {
        hops <- rbind(hops, unlist(lapply(igraph::get.shortest.paths(mst_graph, v)$vpath, length)) - 1)
    }
    v_pos <- array(0, c(vcount(mst_graph), 2))
    print(terminals <- which(hops == max(hops), arr.ind = TRUE)[1,] - 1)
    back_bone <- unlist(igraph::get.shortest.paths(mst_graph, from = terminals["row"], to = terminals["col"])$vpath)
    bb_span <- pi * 0.55
    bb_unit <- 50
    angles <- seq(pi/2 - bb_span/2, pi/2 + bb_span/2, length.out = length(back_bone))
    v_pos[back_bone, ] <- bb_unit * length(back_bone) * cbind(cos(angles), -sin(angles))
    for (v in back_bone) {
        n <- intersect(neighbors(mst_graph, v), back_bone)
        side_v <- sort(subcomponent(delete.edges(mst_graph, E(mst_graph, P = c(mapply(c, v, n)))), v))
        if (length(side_v) > 1) {
            side_h <- hops[v, side_v]
            side_g <- as.directed(subgraph(mst_graph, side_v), mode = "mutual")
            e <- get.edges(side_g, E(side_g))
            side_g <- delete.edges(side_g, subset(E(side_g), side_h[e[, 1]] > side_h[e[, 2]]))
            print(root <- which.min(side_h))
            print(side_g)
            layout <- layout.reingold.tilford(side_g, root = root)
            polar <- cbind(atan2(layout[, 2], layout[, 1]), sqrt(rowSums(layout^2)))
            polar[, 1] <- polar[, 1] + angles[back_bone == v] - pi/2
            layout <- bb_unit * polar[, 2] * cbind(cos(polar[, 1]), -sin(polar[, 1]))
            layout <- layout + matrix(v_pos[v, ] - layout[root, ], nrow = nrow(layout), ncol = 2, byrow = TRUE)
            v_pos[side_v + 1, ] <- layout
        }
    }
    v_pos
}

###
plot.mst <- function( fcs.data,
                      boundary=NULL,
                      layout.table=NULL,
                      param='PSTAT5',
                      pctile_color=c(0.02, 0.98),
                      size_scale_factor=2,
                      edge.color = "grey",
                      palette = "bluered" ) {
    graph <- read.graph('mst.gml',format='gml')
    if (palette == "jet") palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    else if (palette == "bluered") palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
    else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
    colorscale <- palette(100)
    fcs.data2 <- fcs.data
    print(load('clusters.RData'))
    print(nclust <- length(unique(fcs.data$cluster)))
    if (is.null(boundary)) print(boundary <- quantile(fcs.data[,param], probs=pctile_color))
    fcs.data <- fcs.data2
    print(length(unique(clusters <- fcs.data[,'cluster'])))
    # per cluster calculate percent total and median pstat5
    # if cluster does not exist in file then percentotal is zero and median pstat5 is na
    # color of node should be gray and size 1
    percentotal <- numeric(nclust)
    percentotal[sort(unique(clusters))] <- 100*as.numeric(prop.table(table(clusters)))
    #percentotal <- percentotal/sum(percentotal)
    marker <- rep(NA,length=nclust)
    marker[sort(unique(clusters))] <- tapply(fcs.data[,param], clusters, median)
    vsize <- percentotal
    #vsize <- vsize/(max(vsize)^(1/size_scale_factor)) * 3 + 2
    vsize <- vsize/(max(vsize)^(1/size_scale_factor)) * 3 + 1
    vsize[percentotal==0] <- 1
    grad <- seq(boundary[1], boundary[2], length.out = length(colorscale))
    color <- colorscale[findInterval(marker, grad, all.inside = TRUE)]
    color[is.na(color)] <- 'grey'
    fill_color <- color
    frame_color <- color
    #plot(graph, layout=SPADE.layout.arch(graph), vertex.shape="circle", vertex.color=fill_color, vertex.frame.color=rgb(0,0,0,alpha=0), edge.color=edge.color, vertex.size=vsize, vertex.label=NA, edge.arrow.size=0.25, edge.arrow.width=1)
    if (is.null(layout.table)) layout.table <- SPADE.layout.arch(graph)
    plot(layout.table, cex=vsize, col=color, pch=20, yaxt='n', xaxt='n', ann=FALSE,frame.plot=FALSE)
}


###
plot.mst.fold <- function( files,
                      boundaries,
                      layout.table=NULL,
                      param='PSTAT5',
                      size_scale_factor=2,
                      palette = "bluered",
                      DOSES=c('0U','01U','10U','1000U')) {
    if (is.null(layout.table)) {
        graph <- read.graph('mst.gml',format='gml')
        layout.table <- SPADE.layout.arch(graph)
    }
    if (palette == "jet") palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
    else if (palette == "bluered") palette <- colorRampPalette(c("blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red"))
    else stop("Please use a supported color palette.  Options are 'bluered' or 'jet'")
    colorscale <- palette(100)
    print(load('clusters.RData'))
    print(nclust <- length(unique(fcs.data$cluster))) 
    #
    load(files[[1]])
    print(length(unique(clusters <- fcs.data[,'cluster'])))
    marker.base <- rep(NA,length=nclust)
    marker.base[sort(unique(clusters))] <- tapply(fcs.data[,param], clusters, median)
    MARKERS <- matrix(0,nrow=nclust,ncol=length(files))
    MARKERS[,1] <- marker.base
    PERCENT <- matrix(0,nrow=nclust,ncol=length(files))
    #
    for (i in 1:length(files)) {
        print(load(files[[i]]))
        print(length(unique(clusters <- fcs.data[,'cluster'])))
        # per cluster calculate percent total and median pstat5
        # if cluster does not exist in file then percentotal is zero and median pstat5 is na
        # color of node should be gray and size 1
        PERCENT[sort(unique(clusters)),i] <- 100*as.numeric(prop.table(table(clusters)))
        MARKERS[sort(unique(clusters)),i] <- tapply(fcs.data[,param], clusters, median)
        MARKERS[,i] <- MARKERS[,i]-marker.base
    }
    MARKERS[is.na(MARKERS)] <- 0
    grad <- seq(boundaries[[1]], boundaries[[2]], length.out = length(colorscale))
    par(mfrow=c(2,2))
    figure.labels <- iter(paste(letters,')',sep=''))
    for (i in 1:length(files)) {
        vsize <- PERCENT[,i]
        vsize <- vsize/(max(vsize)^(1/size_scale_factor)) * 3 + 1
        vsize[PERCENT[,i]==0] <- 1
        color <- colorscale[findInterval(MARKERS[,i], grad, all.inside = TRUE)]
        color[is.na(color)] <- 'grey'
        plot(layout.table, cex=1, col=color, pch=20, yaxt='n', xaxt='n', ann=FALSE,frame.plot=FALSE)
        title(paste(nextElem(figure.labels), sprintf('pSTAT5 MFI at %s', DOSES[[i]]), sep='\t'), adj=0)
    }
}



###
SPADE.clustering <- function(files,
                             out_dir=".",
                             cluster_cols=NULL,
                             transforms=flowCore::arcsinhTransform(a=0, b=0.2),
                             #downsampling
                             downsampling_samples=20000,
                             downsampling_exclude_pctile=0.01,
                             downsampling_target_pctile=0.05,
                             #clustering
                             clustering_samples=50000,
                             k=300,
                             #density estimation
                             kernel_mult=5.0,
                             apprx_mult=1.5,
                             med_samples=2000) {
    # First step calculate density and downsample files
    for (f in files) {
        message('Computing density on ', f, ' ...')
        load(f)
        #data needs to be transformed for mvt density estimation and clustering
        print(dim(fcs.data <- data.frame(applyTransforms(fcs.data, transforms))))
        #this step can possibly be improved by ANN
        density <- SPADE.density(fcs.data[,cluster_cols], kernel_mult=kernel_mult, apprx_mult=apprx_mult, med_samples=med_samples)
        if (max(density) == 0) warning(paste(f, "has degenerate densities, possibly due to many identical observations", sep = " "))
        print(quantile(fcs.data[,'density'] <- density))
        fcs.data <- SPADE.downsample(fcs.data, exclude_pctile=downsampling_exclude_pctile, target_pctile=downsampling_target_pctile, desired_samples=downsampling_samples)
        #fcs.data now contains 'keep' column
        save(fcs.data, file=file.path(out_dir,basename(f)))
    }
    # Second step: pool files and do clustering
    # From clusters draw MST
    message("Pooling files...")
    #pool downsampled files
    pool.data <- c()
    i <- 1
    for (f in file.path(out_dir, basename(files))) {
        load(f)
        pool.data <- rbind(pool.data, cbind(fcs.data[fcs.data$keep,],file=i))
        i <- i+1
    }
    print(dim(pool.data))
    if (nrow(pool.data) > clustering_samples) pool.data <- pool.data[sample(1:nrow(pool.data), clustering_samples), ]
    else if (nrow(pool.data) == 0) stop("Number of observations to cluster is zero. Did the density/downsampling warn about data similarity?")
    message('Clustering pooled data...')
    clust <- SPADE.cluster(pool.data[,cluster_cols], k)
    print(dim(fcs.data <- subset(cbind(pool.data, cluster = clust$assign), !is.na(clust$assign))))
    cluster.data <- fcs.data
    #
    save(fcs.data, file=file.path(out_dir, "clusters.RData"))
    #
    write.table(clust$centers, file=file.path(out_dir, "clusters.table"), row.names=FALSE, col.names=TRUE)
    ### make MST from cluster centers
    adjacency  <- as.matrix(dist(clust$centers, method='manhattan'))
    full_graph <- graph.adjacency(adjacency,mode="undirected",weighted=TRUE)
    mst_graph  <- minimum.spanning.tree(full_graph)
    # write gml of MST
    write.graph(mst_graph, file.path(out_dir, "mst.gml"), format="gml")
    ### Third step: Upsampling
    ### Assign points in each file to closest cluster.
    ### After this step, each file will have columns density, keep and cluster.
    load(file.path(out_dir, "clusters.RData"))
    cluster.data <- fcs.data
    for (f in file.path(out_dir, basename(files))) {
        message("Upsampling file: ", f)
        print(load(f))
        assign <- .Call( "SPADE_assign",t(fcs.data[,cluster_cols]), t(cluster.data[,cluster_cols]), as.integer(cluster.data[, "cluster"]))
        fcs.data<-data.frame(fcs.data)
        fcs.data[,'cluster'] <- assign
        save(fcs.data, file=f)
    }
}

