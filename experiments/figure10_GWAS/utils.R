suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))
suppressMessages(library(GenomicRanges))

filter_meta <- function(Meta, phenotype) {
    if(phenotype=="bmi") {
        Meta.small <- Meta %>% 
            filter(Phenotype==phenotype) %>%
            filter(grepl(phenotype,`Reported trait`, ignore.case=T) | grepl("body mass index",`Reported trait`, ignore.case=T))
    } else if(phenotype=="cvd") {
        Meta.small <- Meta %>% 
            filter(Phenotype==phenotype) %>%
            filter(grepl(phenotype,`Reported trait`, ignore.case=T) | grepl("cardiovascular disease",`Reported trait`, ignore.case=T))
    } else if(phenotype=="sbp") {
        Meta.small <- Meta %>% 
            filter(Phenotype==phenotype) %>%
            filter(grepl(phenotype,`Reported trait`, ignore.case=T) | grepl("systolic blood pressure",`Reported trait`, ignore.case=T))
    } else if(phenotype=="respiratory") {
        Meta.small <- Meta %>% 
            filter(Phenotype==phenotype) %>%
            filter(grepl(phenotype,`Reported trait`, ignore.case=T) | grepl("Asthma",`Reported trait`, ignore.case=T) | grepl("rhinitis",`Reported trait`, ignore.case=T) | grepl("pulmonary",`Reported trait`, ignore.case=T) | grepl("lung",`Reported trait`, ignore.case=T) )
    } else {
        Meta.small <- Meta %>% 
            filter(Phenotype==phenotype) %>%
            filter(grepl(phenotype,`Reported trait`, ignore.case=T))
    }
    return(Meta.small)
}

single_liftover <- function(chr, bp.min, bp.max, chain) {
    gr <- GRanges(seqnames = sprintf("chr%s",chr), ranges = IRanges(start = bp.min, bp.max))
    gr.new.list <- range(liftOver(gr, chain))
    #gr.new <- do.call("rbind", lapply(gr.new.list, function(gr.new) {
    #    gr.new <- gr.new %>% 
    #        as_tibble() %>%
    #        filter(seqnames==sprintf("chr%s", chr)) %>%
    #        transmute(BP.min=min(start), BP.max=max(end))
    #    return(gr.new)
    #}))
    start <- sapply(start(gr.new.list), function(x) min(x))
    end <- sapply(end(gr.new.list), function(x) max(x))
    gr.new <- tibble(BP.min=start, BP.max=end) %>%
                  mutate(BP.min = ifelse(is.infinite(BP.min), NA, BP.min),
                         BP.max = ifelse(is.infinite(BP.max), NA, BP.max))
    
    if(length(gr)!=nrow(gr.new)) {
        message(sprintf("Liftover error, chr %s!", chr))
    }
    
    if(nrow(gr.new)>0) {
        out <- sprintf("%d-%d", gr.new$BP.min, gr.new$BP.max)
    } else {
        out <- "NA-NA"
    }
    
    return(out)
}

tidy_liftover <- function(Input, source="hg19", target="hg38") {
    # source = "hg19", target = "hg38"
    chain.dir <- "/home/msesia/Workspace/population_structure/meta_associations/liftover"
    if( (source=="hg19") && (target=="hg38") ) {
        chain.path <- sprintf("%s/hg19ToHg38.over.chain", chain.dir)
    } else {
        chain.path <- sprintf("%s/NAToNA.over.chain")
    }
    chain <- import.chain(chain.path)
    
    Output <- Input
    
    chr.list <- Output$CHR
    Output$BP.lifted <-rep("NA-NA", length(chr.list))
    
    chr.list <- sort(unique(Input$CHR))
    for(chr in chr.list) {
        idx.chr <- which(Input$CHR==chr)
        bp.min <- Input$BP.min[idx.chr]
        bp.max <- Input$BP.max[idx.chr]
        Output$BP.lifted[idx.chr] <- single_liftover(chr, bp.min, bp.max, chain)
        message(sprintf("Converted genomic ranges on chromosome %s.", chr), appendLF=FALSE)
     }
    
    Output <- Output %>%
        separate(BP.lifted, "-", into=c("BP.min.lifted","BP.max.lifted"), convert=TRUE)

    return(Output)
}
    
confirm_with_meta <- function(Discoveries, Meta, gap=0) {
    phenotypes <- unique(Discoveries$Phenotype)
    
    output <- lapply(phenotypes, function(phenotype) {
        
        Meta.small <- Meta
        
        Discoveries.2 <- Discoveries %>%
            filter(Phenotype==phenotype)
        
        df.confirmed <- Discoveries.2 %>% 
            left_join(Meta.small, by = c("CHR", "Phenotype")) %>%
            filter(BP<=BP.max.lifted+gap, BP>=BP.min.lifted-gap) %>%
            group_by(Phenotype, Method, Resolution, CHR, Group, r, BP.min, BP.max, SNPs.group, BP.min.lifted, BP.max.lifted) %>%
            summarise(Associations=n(), P=min(P), SNP=paste(unique(SNP),collapse=",")) %>%
            ungroup()

        col.names <- intersect(colnames(Discoveries.2), colnames(df.confirmed))
        df.nonconfirmed <- Discoveries.2 %>%
            anti_join(df.confirmed, by = col.names) %>%
            select_at(vars(col.names)) %>%
            mutate(Associations=0, P=NA, SNP="") %>%
            ungroup()

        df.output <- rbind(df.confirmed, df.nonconfirmed) %>%
            arrange(CHR, BP.min)
        
        return(df.output)
    })
    output <- do.call("rbind", output)
    
    return(output)
}

power_with_meta <- function(Discoveries, Meta, resolution = "425 kb", gap=0) {
    phenotypes <- unique(Discoveries$Phenotype)
    
    output <- lapply(phenotypes, function(phenotype) {
        
        #Meta.small <- filter_meta(Meta, phenotype)
        Meta.small <- Meta
                                  
        Meta.distinct <- Meta.small %>% 
            filter(Phenotype==phenotype, CHR %in% 1:22, CHR!="Mapping not available'") %>%
            distinct(Phenotype, CHR, SNP, BP, .keep_all=TRUE)
        
        Discoveries.2 <- Discoveries %>%
            filter(Phenotype==phenotype, Resolution==resolution)
        
        df.found <- Meta.distinct %>% 
            left_join(Discoveries.2, by = c("CHR", "Phenotype")) %>%
            filter(BP<=BP.max.lifted+gap, BP>=BP.min.lifted-gap) %>%
            select_at(vars(colnames(Meta.distinct))) %>%
            mutate(Found=TRUE, Resolution=resolution) %>%
            ungroup()

        df.notfound <- Meta.distinct %>% 
            anti_join(df.found, by = c("Phenotype", "CHR", "SNP", "BP")) %>%
            select_at(vars(colnames(Meta.distinct))) %>%
            mutate(Found=FALSE, Resolution=resolution) %>%
            ungroup()
        
        df.out <- rbind(df.found, df.notfound) %>% arrange(P)
        return(df.out)
    })
    output <- do.call("rbind", output)
    
    return(output)
}
                  
hinge_exp <- function(x, c = 2) {
    out <- c * log(1 / (c * (1 - x)))
    out[out <= 0 ] <- 0
    return(out)
}


accum_test <- function(x, alpha = .2, c = 2, strict = FALSE) {
    if(strict) {
        accum_val <- (c+cumsum(hinge_exp(x, c))) / seq(2, 1+length(x))
    } else {
        accum_val <- cumsum(hinge_exp(x, c)) / seq(1, length(x))
    }
    below <- which(accum_val <= alpha)
    if(length(below)==0) {
        return(c())
    } else {
        return(1:max(below))
    }
}