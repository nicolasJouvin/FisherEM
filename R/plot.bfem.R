## For CRAN check ...
## https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when/12429344#12429344
utils::globalVariables(names = c('elbo', 'iteration', 'Cluster',
                                 'x', 'y', '.id', 'sigma', 'V1', 'V2'))


#' Plotting function
#'
#' Utility function to plot the results of the BFEM algorithm. The S3 plot
#' function is a wrapper function over the 3 other functions
#' @param res The results of \code{\link{bfem}}.
#' @param type The plot type: \itemize{ \item "subspace" (default) - Uses
#'   \code{plot_subspace()} to plot the projected data \item "criterion" - Uses
#'   \code{plot_crit()} to plot the criterion value. \item "elbo" - Uses
#'   \code{plot_bound()} to plot the variational lower bound evolution. }
#' @param ... Additional parameter to pass to corresponding functions:
#' \itemize{
#' \item crit - Used to specify which criterion should be plotted. Possible values are "aic", "bic" and 'icl. The default is the criterion used in the algorithm.
#' \item alpha_levels - A vector giving the desired Gaussian ellipses level set. Default to 0.95.
#' \item plot.dims - The dimension to be plotted. Default to the first two dimensions.
#' \item show.ellipses - Should Gaussian ellipses be plotted. Default to TRUE
#' \item show.uncertainty - Should uncertainty be plotted. A point is considered uncertain if its posterior probability of membership is peaked toward 2 or more clusters. Graphically, it can be displayed with a bigger point size depending on the uncertainty level, bigger points being more uncertain.
#' \item size - The point size.
#' \item cex.uncertainty - The multiplicative factor for the basic point size controlling the size of uncertain points.
#' }
#' @return a ggplot2 plot object
#' @export
#'
#' @examples
#' \donttest{
#' data(iris)
#' }
plot.bfem <- function(res, type = "subspace", ...) {
  gg = switch (type,
    "subspace" = plot_subspace(res, ...),
    "elbo" = plot_bound(res),
    "criterion" = plot_crit(res, ...)
  )
  gg
}




#--- Plot results from BFEM


bubble <- function(x, cex = c(0.2, 3), alpha = c(0.1, 1)) 
{
  # Uncertainty computing function from Mclust
  # Return size and transparency for points
  x <- as.vector(x)
  cex <- cex[!is.na(cex)]
  alpha <- alpha[!is.na(alpha)]
  x <- (x - min(x))/(max(x) - min(x) + sqrt(.Machine$double.eps))
  n <- length(x)
  r <- sqrt(x/pi)
  r <- (r - min(r, na.rm = TRUE))/
    (max(r, na.rm = TRUE) - min(r, na.rm = TRUE) + sqrt(.Machine$double.eps))
  cex <- r * diff(range(cex)) + min(cex)
  alpha <- x * diff(range(alpha)) + min(alpha)
  return(list(cex = cex, alpha = alpha))
}

#' @describeIn plot.bfem Plot Y projected on the `plot.dims` dimensions of the latent space
# along with gaussian ellipses corresponding to p(\mu_k)
#' @export
plot_subspace <- function(res, alpha_levels = c(0.95), plot.dims = c(1,2), 
                          show.ellipses = T, show.uncertainty = T, 
                          size = 2, cex.uncertainty = 1) {
  
  if (sum(plot.dims > res$d) != 0) stop('Plot dimensions must be < d')
  X_est = res$proj
  df = as.data.frame(X_est[,plot.dims])
  df$Cluster = as.factor(res$cl)
  ndim = length(plot.dims)
  if (res$d == 1) {
    message('Not implemented for d=1.')
    return(NULL)
  }
  if(ndim == 2) {
    gg = ggplot2::ggplot(df, 
                         ggplot2::aes(x = V1, y = V2, col = Cluster, shape = Cluster))
    
    if (show.uncertainty) {
      # inspired from Mclust uncertainty, adapted for ggplot2
      z = res$P
      uncertainty <- 1 - apply(z, 1, max)
      u = (uncertainty - min(uncertainty))/(max(uncertainty) - min(uncertainty)
                                            + sqrt(.Machine$double.eps))
      b <- bubble(u, cex = cex.uncertainty * c(size, size+2), alpha = c(1, 0.7))    
      
      gg = gg + 
        ggplot2::geom_point(size = b$cex,
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
      for (k in 1:res$K) {
        m <- res$var_param$Varmeank[plot.dims,k]
        sigma <- res$param$Sigmak[plot.dims,plot.dims,k]
        contour_data = plyr::ldply(alpha_levels, 
                                   function(level) 
                                     ellipse::ellipse(level = level, x = sigma, centre = m))
        
        contour_data$Cluster = as.factor(k)
        gg = gg +
          ggplot2::geom_path(data=contour_data, 
                             ggplot2::aes(x, y, group = .id, col = Cluster), 
                    linetype='dotdash',
                    alpha = rep(length(alpha_levels):1,each = 100)/1.5 ,
                    show.legend = F)
      }
      
      means = as.data.frame(t(res$var_param$Varmeank))
      means$Cluster = as.factor(1:res$K)
      gg = gg + 
        ggplot2::geom_point(data=means, 
                            ggplot2::aes(x = V1, y = V2),shape = 3, size = 5, show.legend = F)
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

#' @describeIn plot.bfem plot the variational bound evolution 
#' @export
plot_bound = function(res) {
  elbos = res$elbos
  df = data.frame(iteration = 1:length(elbos), elbo = elbos)
  gg = ggplot2::ggplot(df, ggplot2::aes(x=iteration, y=elbo)) +
    ggplot2::geom_point(size=1.5) + 
    ggplot2::geom_line(size=1, alpha = 0.7, linetype='dashed') +
    ggplot2::xlab('Iteration') + ggplot2::ylab('Elbo') +
    ggplot2::theme(text = element_text(size=20)) 
  gg
}

#' @describeIn plot.bfem Plot the criterion result 
#' @export
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
