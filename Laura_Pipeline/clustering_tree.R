# Created by Luke Zappia, githubID: lazappi
#
#' Get tree nodes
#'
#' Get the nodes required for a clustering tree
#'
#' @param clusterings data.frame describing clusterings at different
#' resolutions. Columns must be named "res.X" where X is numeric and indicates
#' the clustering resolution.
#' @param tSNE data.frame giving the t-SNE positions of each cell
#'
#' @return data.frame contining the nodes of a clustering tree
getTreeNodes <- function(clusterings, tSNE = NULL) {

    nodes <- lapply(colnames(clusterings), function(res) {
        clustering <- clusterings[, res]
        clusters <- sort(unique(clustering))

        node <- sapply(clusters, function(cluster) {
            is.cluster <- clustering == cluster
            size <- sum(is.cluster)
            res.clean <- stringr::str_replace(res, "res.", "")
            node.name <- paste0("R", res.clean, "C", cluster)

            node.data <- c(Node    = node.name,
                           Res     = res.clean,
                           Cluster = cluster,
                           Size    = size)

            if (!is.null(tSNE)) {
                node.data["tSNE1"] <- mean(tSNE[is.cluster, 1])
                node.data["tSNE2"] <- mean(tSNE[is.cluster, 2])
            }

            return(node.data)
        })

        return(t(node))
    })

    nodes <- do.call("rbind", nodes) %>%
        dplyr::as_tibble() %>%
        dplyr::mutate_at(c("Res", "Cluster", "Size"), as.numeric)

    if (!is.null(tSNE)) {
        nodes <- nodes %>%
            dplyr::mutate_at(c("tSNE1", "tSNE2"), as.numeric)
    }

    return(nodes)
}

#' Get tree edges
#'
#' Get the edges required for a clustering tree
#'
#' @param clusterings data.frame describing clusterings at different
#' resolutions. Columns must be named "res.X" where X is numeric and indicates
#' the clustering resolution.
#' @param nodes data.frame describing tree nodes
#' @param tSNE logical. Whether to add tSNE coordinates.
#'
#' @return data.frame containing the edges of a clustering tree
getTreeEdges <- function(clusterings, nodes, tSNE = FALSE) {

    res.values <- colnames(clusterings)
    res.clean <- stringr::str_replace(res.values, "res.", "")
    res.values <- res.values[order(as.numeric(res.clean))]

    edges <- lapply(1:(ncol(clusterings) - 1), function(i) {

        from.res <- res.values[i]
        to.res <- res.values[i + 1]

        from.clusters <- sort(unique(clusterings[, from.res]))
        to.clusters <- sort(unique(clusterings[, to.res]))

        trans.df <- expand.grid(FromClust = from.clusters,
                                ToClust = to.clusters, stringsAsFactors = FALSE)

        trans <- apply(trans.df, 1, function(x) {
            from.clust <- x[1]
            to.clust <- x[2]

            is.from <- clusterings[, from.res] == from.clust
            is.to <- clusterings[, to.res] == to.clust

            trans.count <- sum(is.from & is.to)

            from.size <- sum(is.from)
            to.size <- sum(is.to)
            trans.prop.from <- trans.count / from.size
            trans.prop.to <- trans.count / to.size

            return(c(trans.count, trans.prop.from, trans.prop.to))
        })

        trans.df$FromRes <- as.numeric(gsub("res.", "", from.res))
        trans.df$ToRes <- as.numeric(gsub("res.", "", to.res))
        trans.df$TransCount <- trans[1, ]
        trans.df$TransPropFrom <- trans[2, ]
        trans.df$TransPropTo <- trans[3, ]

        return(trans.df)
    })

    edges <- dplyr::bind_rows(edges)

    levs <- sort(unique(edges$ToClust))
    edges <- edges %>%
        dplyr::mutate(FromClust = factor(FromClust, levels = levs))  %>%
        dplyr::mutate(ToClust = factor(ToClust, levels = levs)) %>%
        dplyr::mutate(FromNode = paste0("R", FromRes, "C", FromClust)) %>%
        dplyr::mutate(ToNode = paste0("R", ToRes, "C", ToClust)) %>%
        dplyr::left_join(nodes, by = c(FromNode = "Node"))

    if (tSNE) {
        edges <- edges %>%
        dplyr::rename(tSNE1From = tSNE1, tSNE2From = tSNE2)
    }

    edges <- edges %>%
        dplyr::left_join(nodes, by = c(ToNode = "Node"))

    if (tSNE) {
        edges <- edges %>%
        dplyr::rename(tSNE1To = tSNE1, tSNE2To = tSNE2)
    }

    edges <- edges %>%
        dplyr::select(matches("From"), matches("To"), TransCount)

    return(edges)
}

#' Plot clustering tree
#'
#' Plot a tree showing the relationshop between clusterings at different
#' resolutions
#'
#' @param clusterings data.frame describing clusterings at different
#' resolutions. Columns must be named "res.X" where X is numeric and indicates
#' the clustering resolution.
#' @param count.filter Threshold for the number of cells associated with an edge
#' @param prop.filter Threhold for the to cluster proportion associated with an
#'        edge
#'
#' @return ggplot of the clustering tree
plotClusteringTree <- function(clusterings, count.filter = 0,
                               prop.filter = 0.1) {

    library("ggraph")

    nodes <- getTreeNodes(clusterings)
    edges <- getTreeEdges(clusterings, nodes)

    edges %>%
        dplyr::filter(TransCount > count.filter) %>%
        dplyr::filter(TransPropTo > prop.filter) %>%
        dplyr::select(FromNode, ToNode, everything()) %>%
        igraph::graph_from_data_frame(vertices = nodes) %>%
    ggraph(layout = "tree") +
        geom_edge_link(arrow = arrow(length = unit(1, 'mm')),
                       end_cap = circle(3.5, "mm"), edge_width = 1,
                       aes(colour = log(TransCount), alpha = TransPropTo)) +
        geom_node_point(aes(colour = factor(Res),
                            size = Size)) +
        geom_node_text(aes(label = Cluster), size = 3) +
        scale_size(range = c(4, 15)) +
        scale_edge_colour_gradientn(colours = viridis::viridis(100)) +
        guides(size = guide_legend(title = "Cluster Size",
                                   title.position = "top"),
               colour = guide_legend(title = "Clustering Resolution",
                                     title.position = "top"),
               edge_colour = guide_edge_colorbar(title = "Cell Count (log)",
                                                 title.position = "top"),
               edge_alpha = guide_legend(title = "Cluster Prop",
                                         title.position = "top", nrow = 2)) +
        cowplot::theme_nothing() +
        theme(legend.position = "bottom")
}

#' Plot clustering tree
#'
#' Plot a tree showing the relationshop between clusterings at different
#' resolutions
#'
#' @param seurat Seurat object that has been clustered at different resolutions
#' @param count.filter Threshold for the number of cells associated with an edge
#' @param prop.filter Threhold for the to cluster proportion associated with an
#'        edge
#' @param levels Resolutions to included. If NULL all are used.
#'
#' @return ggplot of the clustering tree
plotClusteringTreeSeurat <- function(seurat, count.filter = 0,
                                     prop.filter = 0.1, levels = NULL) {

    clusterings <- seurat@meta.data %>% dplyr::select(dplyr::contains("res."))

    if (!is.null(levels)) {
        clusterings <- clusterings[, paste0("res.", levels)]
    }

    plotClusteringTree(clusterings, count.filter = count.filter,
                       prop.filter = prop.filter)
}

#' Plot clustering tree-SNE
#'
#' Plot a clustering tree projected onto t-SNE coordinates
#'
#' @param clusterings data.frame describing clusterings at different
#' resolutions. Columns must be named "res.X" where X is numeric and indicates
#' the clustering resolution.
#' @param tSNE data.frame giving the t-SNE positions of each cell
#' @param count.filter Threshold for the number of cells associated with an edge
#' @param prop.filter Threhold for the to cluster proportion associated with an
#'        edge
#' @param plot.sides Whether to add plots from the point of view of the t-SNE
#'        dimensions
#'
#' @return ggplot of the clustering tree-SNE
plotClusteringTreeSNE <- function(clusterings, tSNE, count.filter = 0,
                                  prop.filter = 0.1, plot.sides = FALSE) {

    library("ggplot2")

    nodes <- getTreeNodes(clusterings, tSNE)
    edges <- getTreeEdges(clusterings, nodes, tSNE = TRUE) %>%
        dplyr::filter(TransCount > count.filter) %>%
        dplyr::filter(TransPropTo > prop.filter)

    top <- ggplot(as_tibble(tSNE), aes(x = tSNE_1, y = tSNE_2)) +
        geom_point(alpha = 0.1) +
        geom_point(data = nodes,aes(x = tSNE1, y = tSNE2,
                                     colour = factor(Res), size = Size)) +
        geom_segment(data = edges,
                     aes(x = tSNE1From, y = tSNE2From,
                         xend = tSNE1To, yend = tSNE2To,
                         alpha = TransPropTo, colour = factor(FromRes)),
                     arrow = arrow(length = unit(0.02, "npc")),
                     size = 1) +
        scale_size(range = c(3, 15)) +
        viridis::scale_colour_viridis(discrete = TRUE) +
        guides(size = guide_legend(title = "Size",
                                   title.position = "top"),
               colour = guide_legend(title = "Resolution",
                                     title.position = "top"),
               alpha = guide_legend(title = "Proportion",
                                    title.position = "top", nrow = 2)) +
        cowplot::theme_cowplot()

    if (plot.sides) {
        side1 <- plotTreeSNESide(tSNE, nodes, edges, side = 1)
        side2 <- plotTreeSNESide(tSNE, nodes, edges, side = 2)
        sides <- cowplot::plot_grid(side1, side2, ncol = 1)
        legend <- cowplot::get_legend(top + theme(legend.position = "bottom"))
        treeSNE <- cowplot::plot_grid(top + theme(legend.position = "none"),
                                      sides, ncol = 2, rel_widths = c(2, 1))
        treeSNE <- cowplot::plot_grid(treeSNE, legend, ncol = 1,
                                      rel_heights = c(4, 1))
    } else {
        treeSNE <- top
    }

    return(treeSNE)
}

#' Plot clustering tree-SNE
#'
#' Plot a clustering tree projected onto t-SNE coordinates
#'
#' @param seurat Seurat object that has been clustered at different resolutions
#' @param count.filter Threshold for the number of cells associated with an edge
#' @param prop.filter Threhold for the to cluster proportion associated with an
#'        edge
#' @param plot.sides Whether to add plots from the point of view of the t-SNE
#'        dimensions
#' @param levels Resolutions to included. If NULL all are used.
#'
#' @return ggplot of the clustering tree-SNE
plotClusteringTreeSNESeurat <- function(seurat, count.filter = 0,
                                        prop.filter = 0.1, plot.sides = FALSE,
                                        levels = NULL) {

    clusterings <- seurat@meta.data %>% dplyr::select(dplyr::contains("res."))
    tSNE <- seurat@dr$tsne@cell.embeddings

    if (!is.null(levels)) {
        clusterings <- clusterings[, paste0("res.", levels)]
    }

    plotClusteringTreeSNE(clusterings, tSNE, count.filter = count.filter,
                          prop.filter = prop.filter, plot.sides = plot.sides)
}

#' Plot clustering tree-SNE side
#'
#' Plot a clustering tree-SNE where the x-axis is a t-SNE dimension and the
#' y-axis is clustering resolution
#'
#' @param tSNE data.frame giving the t-SNE coordinates of each
#' @param nodes data.frame describing the nodes of the clustering tree
#' @param edges data.frame describing the edges of the clustering tree
#' @param side t-SNE dimension to use as the x-axis
#'
#' @return ggplot of the side view of a clustering tree-SNE
plotTreeSNESide <- function(tSNE, nodes, edges, side = 1) {

    nodes$ResFactor <- factor(nodes$Res)
    edges$FromResFactor <- factor(edges$FromRes)

    ggplot(as_tibble(tSNE), aes_string(x = paste0("tSNE_", side), y = 1)) +
        geom_point(alpha = 0.1) +
        geom_point(data = nodes,
                   aes_string(x = paste0("tSNE", side), y = "Res",
                              colour = "ResFactor", size = "Size")) +
        geom_segment(data = edges,
                     aes_string(x = paste0("tSNE", side, "From"),
                                y = "FromRes",
                                xend = paste0("tSNE", side, "To"),
                                yend = "ToRes",
                                alpha = "TransPropTo",
                                colour = "FromResFactor"),
                     arrow = arrow(length = unit(0.02, "npc")),
                     size = 1) +
        scale_size(range = c(1, 8)) +
        viridis::scale_colour_viridis(discrete = TRUE) +
        scale_y_reverse() +
        ylab("Resolution") +
        cowplot::theme_cowplot() +
        theme(legend.position = "none")
}
