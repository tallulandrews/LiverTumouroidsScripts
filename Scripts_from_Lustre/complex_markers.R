ctp_fast_AUC <- function (expression_vec, truth) {
    R = rank(expression_vec)
    N1 = sum(truth)
    N2 = sum(!truth)
    U2 = sum(R[!truth]) - N2 * (N2 + 1)/2
    if (N1 == 0) {
        return(c(0, 0, 0))
    }
    if (N2 == 0) {
        return(c(1, 1, 1))
    }
    AUC = 1 - U2/(N1 * N2)
    Q1 = AUC/(2 - AUC)
    Q2 = 2 * AUC^2/(1 + AUC)
    SE = sqrt((AUC * (1 - AUC) + (N1 - 1) * (Q1 - AUC^2) + (N2 - 
        1) * (Q2 - AUC^2))/(N1 * N2))
    return(c(max(0, AUC - 1.96 * SE), AUC, min(1, AUC + 1.96 * 
        SE)))
}

complex_markers <- function (expr_mat, labels, n_max = length(unique(labels)) -1, strict_only = FALSE) {
    if (length(labels) != length(expr_mat[1, ])) {
        stop("Length of labels does not match number of cells.")
    }
    if (n_max > (length(unique(labels)) - 1)) {
        stop("n_max must be less than the number of labels")
    }
    label_counts = factor_counts(labels)
    exclude = names(label_counts)[label_counts < 2]
    if (length(exclude) > 0) {
        warning(paste("Warning: Excluding", length(exclude), 
            "groups with less than 2 samples."))
        keep = !(labels %in% exclude)
        expr_mat <- expr_mat[, keep]
        labels <- labels[keep]
        labels <- factor(labels)
    }
    if (min(label_counts) < 10) {
        print("Warning: Small groups (n < 10) may bias marker gene results.")
    }
    gene_cluster_means <- function(mat, groups) {
        MAT <- as.matrix(mat)
        x <- split(seq(ncol(MAT)), groups)
        result <- sapply(x, function(a) rowMeans(MAT[, a]))
        return(result)
    }
    ranked_matrix <- t(apply(expr_mat, 1, rank))
    gene_cluster_ranks <- gene_cluster_means(ranked_matrix, labels)
    cluster_priority <- t(apply(-1 * gene_cluster_ranks, 1, rank))
    label_columns <- sort(unique(labels))
    gene_auc <- function(g) {
        if (sum(cluster_priority[g, ] == cluster_priority[g, 
            1]) == length(cluster_priority[g, ])) {
            group = rep(0, times = length(label_columns))
            pval = -1
            auc = -1
            return(c(auc, group, pval))
        }
        get_auc_ci <- function(n) {
            g_groups <- colnames(cluster_priority)[which(cluster_priority[g, 
                ] <= n)]
            if (length(g_groups) < length(unique(labels)) & length(g_groups) > 
                0) {
                return(ctp_fast_AUC(expr_mat[g, ], labels %in% 
                  g_groups))
            }
            else {
                return(c(0, 0.5, 1))
            }
        }
        auc_tab <- sapply(1:n_max, get_auc_ci)
        top = which(auc_tab[2, ] == max(auc_tab[2, ]))
        sec = which(auc_tab[2, ] == max(auc_tab[2, -top]))
        if (n_max > 1 & min(auc_tab[1, top]) <= max(auc_tab[2, 
            sec]) & strict_only) {
            group = rep(0, times = length(label_columns))
            pval = -1
            auc = -1
        }
        else {
            n = max(top)
            g_groups = colnames(cluster_priority)[which(cluster_priority[g, 
                ] <= n)]
            if (n_max == 1) {
                group = paste(g_groups, collapse = "+")
            }
            else {
                group = as.numeric(label_columns %in% g_groups)
            }
            auc = auc_tab[2, n]
            pval = wilcox.test(expr_mat[g, !(labels %in% g_groups)], 
                expr_mat[g, (labels %in% g_groups)])$p.value
        }
        return(c(auc, group, pval))
    }
    out <- sapply(1:length(expr_mat[, 1]), gene_auc)
    out_matrix = as.data.frame(t(out))
    rownames(out_matrix) <- rownames(expr_mat)
    if (n_max == 1) {
        colnames(out_matrix) <- c("AUC", "Group", "p.value")
    }
    else {
        colnames(out_matrix) <- c("AUC", as.character(label_columns), 
            "p.values")
    }
    out_matrix[, 1] = as.numeric(as.character(out_matrix[, 1]))
    out_matrix[, 3] = as.numeric(as.character(out_matrix[, 3]))
    out_matrix$q.value = out_matrix$p.value * length(unique(labels)) * 
        length(expr_mat[, 1])
    out_matrix$q.value[out_matrix$q.value < 0] = -1
    out_matrix$q.value[out_matrix$q.value > 1] = 1
    return(out_matrix)
}

