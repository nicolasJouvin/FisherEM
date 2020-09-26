plot.bfem <- function(res, type = "subspace", ...) {
  gg = switch (type,
    "subspace" = plot_subspace(res, ...),
    "elbo" = plot_bound(res),
    "crit" = plot_crit(res, ...)
  )
  gg
}


# Plot Y projected on the plot.dims dimensions of the latent space
# along with gaussian ellipses corresponding to p(\mu_k)
plot_subspace <- function(res, alpha_levels = c(0.95), plot.dims = c(1,2), 
                          show.ellipses = T, show.uncertainty = T, 
                          size = 2, uncertainty.factor = 1) {
  
  if (sum(plot.dims > res$d) != 0) stop('Plot dimensions must be < d')
  X_est = res$proj
  df = as.data.frame(X_est[,plot.dims])
  df$Cluster = as.factor(res$cl)
  ndim = length(plot.dims)
  
  if(ndim == 2) {
    gg = ggplot2::ggplot(df, aes(x = V1, y = V2, col = Cluster, shape = Cluster))
    
    if (show.uncertainty) {
      # inspired from Mclust uncertainty, adapted for ggplot2
      z = res$var_param$tau
      uncertainty <- 1 - apply(z, 1, max)
      u = (uncertainty - min(uncertainty))/(max(uncertainty) - min(uncertainty)
                                            + sqrt(.Machine$double.eps))
      b <- bubble(u, cex = uncertainty.factor * c(size, size+2), alpha = c(1, 0.7))    
      
      gg = gg + 
        geom_point(size = b$cex,
                   alpha = max(b$alpha) - b$alpha + min(b$alpha))
    } else {
      gg = gg +
        ggplot2::geom_point(size = size)
    }
    
    
    
    # -- code for ellipse
    if (show.ellipses) {
      alpha_levels = sort(alpha_levels)
      names(alpha_levels) <- alpha_levels ## to get id column in result
      contour_data = NULL
      for (k in 1:res$k) {
        m <- res$var_param$Varmeank[plot.dims,k]
        sigma <- res$param$Sigmak[plot.dims,plot.dims,k]
        contour_data = plyr::ldply(alpha_levels, 
                                   function(level) 
                                     ellipse::ellipse(level = level, x = sigma, centre = m))
        
        contour_data$Cluster = as.factor(k)
        gg = gg +
          ggplot2::geom_path(data=contour_data, aes(x, y, group = .id, col = Cluster), 
                    linetype='dotdash',
                    alpha = rep(length(alpha_levels):1,each = 100)/1.5 ,
                    show.legend = F)
      }
      
      means = as.data.frame(t(res$var_param$Varmeank))
      means$Cluster = as.factor(1:res$k)
      gg = gg + 
        ggplot2::geom_point(data=means, aes(x = V1, y = V2),shape = 3, size = 5, show.legend = F)
    }
    
    # Colorblind friendly palette
    gg = gg +
      ggplot2::scale_color_brewer(palette="Set2") +
      ggplot2::xlab(paste0('U', plot.dims[1])) + ggplot2::ylab(paste0('U', plot.dims[2]))
    return(gg)
  } else {
    base::message('ndim > 2 not implemented yet')
    return(NULL)
  }
}

plot_bound = function(res) {
  elbos = res$elbos
  df = data.frame(iteration = 1:length(elbos), elbo = elbos)
  gg = ggplot2::ggplot(df, aes(x=iteration, y=elbo)) +
    ggplot2::geom_point(size=1.5) + 
    ggplot2::geom_line(size=1, alpha = 0.7, linetype='dashed') +
    ggplot2::xlab('Iteration') + ggplot2::ylab('Elbo') +
    ggplot2::theme(text = element_text(size=20)) 
  gg
}


plot_crit = function(res, crit = NULL) {
  color_palette = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  shape_palette = factor(1:12)
  if(is.null(crit)) crit = res$crit
  df = res$allCriteria
  df$K = as.integer(df$K)
  gg = ggplot2::ggplot(df, 
                       ggplot2::aes_string(x="K", y=crit, shape = "model", 
                                  linetype = "model", col = "model")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ylab(base::casefold(crit, upper = T)) +
    scale_shape_manual(values = shape_palette) +
    scale_linetype_manual(values = shape_palette) +
    scale_color_manual(values = color_palette)
  
  gg
}
