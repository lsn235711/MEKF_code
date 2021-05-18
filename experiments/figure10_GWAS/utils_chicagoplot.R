font.size <- 10
title.font.size <- 10
legend.font.size <- 10

plot_chicago <- function(window.chr, window.left, window.right, Discoveries, upside=FALSE, plot.title=NULL, overlay=FALSE) {

    if(is.null(plot.title)) {
        plot.title <- "(b) Chicago plot (KnockoffGWAS)"
    }
    
    # Extract knockoff discoveries within this window
    if(!is.null(Discoveries)) {
        Knockoffs.window <- Discoveries %>% filter(Method=="Knockoffs") %>%
            filter(CHR==window.chr, BP.min<=window.right, BP.max>=window.left)
        cat(sprintf("There are %d knockoff discoveries within this window.\n", nrow(Knockoffs.window)))
    } else {
        Knockoffs.window <- tibble()
    }

    ## Plot knockoff discoveries
    if(upside) {
        resolution.list <- c(6,5,4,3,2,1,0)
        resolution.heights <- seq(length(resolution.list))
        names(resolution.heights) <- resolution.list
        resolution.labels <- c("425", "208", "81", "41", "20", "3", "single-SNP")
    } else {
        resolution.list <- c(6,5,4,3,2,1,0) %>% rev
        resolution.heights <- seq(length(resolution.list))
        names(resolution.heights) <- resolution.list
        resolution.labels <- c("425", "208", "81", "41", "20", "3", "SNP") %>% rev()
    }

    if(overlay) {
        Knockoffs.window <- Knockoffs.window %>%
            ungroup() %>%
            group_by(Method, Resolution, CHR, Group, BP.min, BP.max, FDP.local) %>%
            summarise(r = max(r))
        alpha <- 0.6
    } else {
        alpha <- 0.2
    }
            
    if(nrow(Knockoffs.window)>0) {
        print(Knockoffs.window)
        p.knockoffs <- Knockoffs.window %>%
            mutate(Resolution=as.character(Resolution), Environments=factor(r, 1:4, 1:4)) %>%
            mutate(Height=resolution.heights[Resolution]) %>%
            ggplot() +
            geom_rect(aes(xmin=BP.min, xmax=BP.max, ymin=Height-0.5, ymax=Height+0.5, fill=Environments, color=Environments))

    } else {
        p.knockoffs <- ggplot(tibble()) + geom_blank()
    }
    
    if(overlay) {
        my_theme <- theme(panel.grid.minor.y = element_blank(),
                          axis.line=element_blank(),
                          panel.border=element_blank(),
                          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
                          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
                          text = element_text(size=font.size),
                          axis.title.y = element_text(size=title.font.size),
                          plot.title = element_text(size=title.font.size),
                          legend.text = element_text(size=legend.font.size),
                          legend.title = element_text(size=legend.font.size),
                          legend.key.height = unit(0.75,"line")
                          )
    } else {
        p.knockoffs <- p.knockoffs +
            facet_grid(Environments~., labeller = labeller(Environments = label_both))            
        my_theme <- theme(panel.grid.minor.y = element_blank(),
                          axis.line=element_blank(),
                          panel.border=element_blank(),
                          panel.grid.major.x = element_line(size = 0.2, colour = "darkgray"),
                          panel.grid.minor.x = element_line(size = 0.1, colour = "darkgray"),
                          text = element_text(size=font.size),
                          axis.title.y = element_text(size=title.font.size),
                          plot.title = element_text(size=title.font.size),
                          legend.text = element_text(size=legend.font.size),
                          legend.title = element_text(size=legend.font.size),
                          legend.key.height = unit(0.75,"line"),
                          legend.position="none"                          
                          )
    }

    #color.palette <- c("grey20", "blue3", "red3", "orange1")
    color.palette <-  c("#999999", "#E69F00", "#56B4E9", "#009E73")
    p.knockoffs <- p.knockoffs +
        ylab("Resolution (kb)") + xlab("") +
        coord_cartesian(xlim = c(window.left,window.right)) +
        scale_x_continuous(expand=c(0.01,0.01), labels=bp.labeler) +
        scale_y_continuous(limits=c(0.5,max(resolution.heights)+0.5),
                           labels=resolution.labels, breaks=resolution.heights) +
        scale_fill_manual(labels = c(1:4), values=alpha(color.palette, alpha)) +
        scale_color_manual(labels = c(1:4), values=alpha(color.palette,1)) +
        ggtitle(plot.title) +
        labs(fill = "Environments", color = "Environments") +
        xlab(sprintf("Physical position on chromosome %d (Mb)", window.chr)) +
        theme_bw() +
        my_theme

    return(p.knockoffs)
}
