library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(readr)
library(ggbeeswarm)
theme_set(theme_bw())
options(stringsAsFactors=F)

library(argparser)
p <- arg_parser("ddOWL mutatation allele phasing and plotting tools, v0.1 - Nils Koelling")
p <- add_argument(p, "FAMILIES", help="families file")
p <- add_argument(p, "SNPS", help="snps mask")
p <- add_argument(p, "--prefix", help="output filename prefix")

main = function(argv) {
    if(missing(argv)) {
        argv = commandArgs(trailingOnly = TRUE)
    }

    args <- parse_args(p, argv)

    if(is.na(args$prefix)) {
        args$prefix = args$FAMILIES
    }

    families = read_csv(args$FAMILIES) %>%
        rename(sample = BC) %>%
        mutate(rel_sam = sprintf('%s (%s)', relationship, sample))

    snps = families %>%
        distinct(FamilyID) %>%
        rowwise %>%
        do({
            fam = .$FamilyID

            read_tsv(sprintf(args$SNPS, fam), col_types = cols(
                POS = col_integer(),
                ID = col_character(),
                REF = col_character(),
                ALT = col_character()
            )) %>%
            rename(CHROM = `#CHROM`) %>%
            mutate(
                snp_name = sprintf('%s_%d', CHROM, POS),
                snp_name_or_id = ifelse(!is.na(ID), ID, name),
                FamilyID = fam
            )
        }) %>%
        distinct(FamilyID, snp_name, .keep_all=TRUE)

    #QC PLOTS
    rs = read_csv(sprintf('%s.read_stats.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample')) %>%
        mutate(
            fam_rel_sam = paste(FamilyID, rel_sam),
            mismatch_rate = mismatches / length
        )
    #rs %>% head %>% print

    fn = sprintf('%s.read_stats.pdf', args$prefix)
    pdf(fn)
    print(
        ggplot(rs, aes(fam_rel_sam, fill=FamilyID))
        + geom_bar(stat='count', position='dodge')
        + coord_flip()
        + xlab('Family relationship (sample)')
    )
    print(
        ggplot(rs, aes(fam_rel_sam, length, colour=FamilyID))
        + geom_boxplot()
        #+ geom_violin(aes(fill=FamilyID))
        + xlab('Family relationship (sample)')
        + coord_flip()
    )
    print(
        ggplot(rs, aes(fam_rel_sam, mismatch_rate, colour=FamilyID))
        + geom_boxplot()
        #+ geom_violin(aes(fill=FamilyID))
        + xlab('Family relationship (sample)')
        + coord_flip()
    )
    print(
        ggplot(rs, aes(fam_rel_sam, mapping_quality, colour=FamilyID))
        + geom_boxplot()
        #+ geom_violin(aes(fill=FamilyID))
        + xlab('Family relationship (sample)')
        + coord_flip()
    )
    print(
        ggplot(rs, aes(fam_rel_sam, mean_baseq, colour=FamilyID))
        + geom_boxplot()
        #+ geom_violin(aes(fill=FamilyID))
        + xlab('Family relationship (sample)')
        + coord_flip()
    )
    dev.off()
    message(fn)

    #background
    bc = read_csv(sprintf('%s.read_background_calls.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample'))
    var_ranges = bc %>%
        distinct(CHROM, POS)
    var_ranges = GenomicRanges::GRanges(var_ranges$CHROM, IRanges::IRanges(var_ranges$POS, var_ranges$POS))
    var_ranges

    GenomicRanges::elementMetadata(var_ranges)$ref = as.character(BSgenome::getSeq(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38, var_ranges))

    bcs = bc %>% left_join(
        var_ranges %>% data.frame %>% mutate(CHROM = as.character(seqnames)) %>% select(CHROM, POS = start, ref),
        by = c('CHROM', 'POS')
    ) %>% mutate(
        match = snp_call == ref
    )

    cat('REF match percentage:\n')
    print(table(bcs$match) / nrow(bcs))

    bcsp = bcs %>%
        replace_na(list(snp_call_quality = 0)) %>%
        mutate(
            qual_bin = cut(snp_call_quality, 20)
        ) %>%
        group_by(FamilyID, sample, qual_bin) %>%
        summarise(p_error_empirical = 1 - sum(match) / n()) %>%
        mutate(
            qual_bin_clean = gsub('\\(|\\]', '', qual_bin)
        ) %>%
        tidyr::separate(qual_bin_clean, sep = ',', into = c('start', 'end')) %>%
        mutate(
            start = as.numeric(start),
            end = as.numeric(end),
            mid = start + (end - start) / 2,
            p_error_theoretical = 10^(-mid/10),
            start = ifelse(start == min(start), -Inf, start),
            end = ifelse(end == max(end), Inf, end)
        )

    fn = sprintf('%s.quality_scores.pdf', args$prefix)
    pdf(fn)

    print(
        ggplot(bcs, aes(sample, fill = match))
        + geom_bar(stat = 'count', position = 'fill')
    )

    print(
        ggplot(bcs %>% group_by(snp_name) %>% filter(n() > 10) %>% ungroup, aes(snp_name, fill = match))
        + geom_bar(stat = 'count', position = 'fill')
        + coord_flip()
        + facet_grid(~sample)
    )

    print(
        ggplot(bcs %>% group_by(snp_name) %>% filter(n() > 10) %>% ungroup, aes(snp_name, fill = match))
        + geom_bar(stat = 'count', position = 'stack')
        + coord_flip()
        + facet_grid(~sample)
    )

    print(
        ggplot(bcs %>% replace_na(list(snp_call_quality = 0)), aes(snp_call_quality, fill = match))
        + geom_histogram()
        + facet_grid(ref~sample)
    )

    print(
        ggplot(bcs %>% replace_na(list(snp_call_quality = 0)), aes(match, snp_call_quality))
        + geom_boxplot()
        + facet_grid(ref~sample)
    )

    print(
        bcs %>%
            replace_na(list(snp_call_quality = 0)) %>%
            mutate(
                qual_bin = cut(snp_call_quality, 20)
            ) %>%
            ggplot
            + aes(qual_bin, fill = match)
            + geom_bar(stat = 'count', position = 'fill')
            + coord_flip()
            + facet_grid(ref~sample)
    )

    print(
        ggplot(bcsp, aes(p_error_theoretical, p_error_empirical, colour=FamilyID)) + geom_point() + geom_abline() + facet_wrap(~sample)
    )

    print(
        ggplot(bcsp, aes(p_error_theoretical, p_error_empirical, colour=FamilyID)) + geom_point() + scale_y_log10() + scale_x_log10() + geom_abline() + facet_wrap(~sample)
    )

    dev.off()
    message(fn)


    #allele counts
    ac = read_csv(sprintf('%s.allele_counts.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample')) %>%
        left_join(snps, by=c('FamilyID', 'snp_name'))
    #ac %>% head %>% print

    fn = sprintf('%s.allele_counts.pdf', args$prefix)
    pdf(fn)
    for(family in unique(ac$FamilyID)) {
        fam_ac = ac %>% filter(FamilyID == family)

        allele_totals = fam_ac %>%
    #        filter(!snp_allele %in% c('_FLAG', '_DEL')) %>%
            group_by(snp_name_or_id, snp_allele) %>%
            summarise(total = sum(count)) %>%
            mutate(
                rank = rank(-total, ties.method = 'first')
            )

        summarised_counts = fam_ac %>%
            group_by(rel_sam, relationship, snp_name_or_id) %>%
            inner_join(allele_totals, by=c('snp_name_or_id', 'snp_allele')) %>%
            mutate(
                allele_simple = ifelse(rank <= 4, snp_allele, '_OTHER'),
                fraction = count / sum(count)
            ) %>%
            group_by(sample, relationship, snp_name_or_id, allele_simple) %>%
            summarise(
                count = sum(count),
                fraction = sum(fraction)
            )

        plots = summarised_counts %>%
            group_by(snp_name_or_id) %>%
            do(
                plot_fractions = (
                    ggplot(., aes(sample, fraction, fill=allele_simple))
                    + geom_bar(stat='identity', position='stack')
                    + scale_fill_brewer(type = 'qual', palette = 'Dark2')
                    + ggtitle(family, sprintf('%s fraction', .$snp_name_or_id))
                    + scale_y_continuous(labels=scales::percent)
                    + xlab('Sample')
                    + ylab('Percent of reads')
                    + coord_flip()
                    + theme(legend.position = 'none')
                ),
                plot_counts = (
                    ggplot(., aes(sample, count, fill=allele_simple))
                    + geom_bar(stat='identity', position='stack')
                    + scale_fill_brewer(type = 'qual', palette = 'Dark2', name='Call')
                    + ggtitle(family, sprintf('%s count', .$snp_name_or_id))
                    + xlab('Sample')
                    + ylab('Reads')
                    + coord_flip()
                    + theme(legend.position = 'bottom')
                )
            ) %>%
            do(
                grid = plot_grid(.$plot_fractions, .$plot_counts, nrow=2)
            )
        lapply(plots$grid, print)
    }
    dev.off()
    message(fn)

    #Actual phasing
    pe = read_csv(sprintf('%s.phase_evidence.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample')) %>%
        left_join(snps, by=c('FamilyID', 'snp_name'))
    #pe %>% head %>% print

    fn = sprintf('%s.phase_evidence.pdf', args$prefix)
    pdf(fn, width=10, height=5)
    for(family in unique(pe$FamilyID)) {
        r_full_df = pe %>%
            filter(FamilyID == family) %>%
            filter(snp_allele != '_NONE')

        for(snp in unique(r_full_df$snp_name_or_id)) {
            snp_data = r_full_df %>%
                filter(snp_name_or_id == snp)

            #use top3 for snp allele type if we have no maternal/paternal info
            if(all(is.na(snp_data$snp_allele_type))) {
                #FIXME: we do not actually need this, because we simplify to REF/ALT/_OTHER anyway!
                top_alleles = snp_data %>%
                    group_by(sample, relationship) %>%
                    mutate(rank = rank(-count, ties.method = 'first')) %>%
                    group_by(snp_allele) %>%
                    summarise(meanrank = mean(rank)) %>%
                    top_n(4, wt = meanrank)
                #print(top_alleles)

                #now filter down to top alleles and set allele type to raw allele
                snp_data = snp_data %>%
                    inner_join(top_alleles, by = 'snp_allele') %>%
                    mutate(snp_allele_type = snp_allele)
            }

            snp_data = snp_data  %>%
                filter(!is.na(snp_allele_type) & !is.na(mutation_allele_type)) %>%
                mutate(
                    snp_allele_type = factor(snp_allele_type),
                    mutation_allele_type = factor(mutation_allele_type, levels = c('mutation', 'wild-type'))
                )

            #print(snp_data)

            snp_data = snp_data %>%
                complete(
                    nesting(sample, relationship),
                    snp_allele_type, mutation_allele_type,
                    fill = list(count = 0)
                )
            #print(snp_data)

            if(nrow(snp_data) > 0) {
                print(
                    do.call(plot_grid,
                        c(
                            snp_data %>%
                                group_by(sample, relationship) %>%
                                do(plot = {
                                    (
                                        ggplot(., aes(mutation_allele_type, snp_allele_type, fill = count, label = count))
                                        + geom_tile()
                                        + geom_text(aes(colour = count > mean(count)), size = 3)
                                        + scale_y_discrete(drop=FALSE)
                                        + scale_x_discrete(drop=FALSE)
                                        + coord_fixed()
                                        + scale_colour_manual(guide = FALSE, values = c(`TRUE` = "black", `FALSE` = "white"))
                                        #+ viridis::scale_fill_viridis(name = 'Number of reads')
                                        + scale_fill_continuous(low = 'darkblue', high = 'yellow')
                                        + ggtitle(sprintf('%s %s', family, unique(.$relationship)), sprintf('%s / %s', snp, unique(.$sample)))
                                        + theme(
                                            axis.title = element_blank(),
                                            legend.position = 'bottom',
                                            legend.text = element_text(size = 6)
                                        )
                                    )
                                }) %>%
                                .[['plot']],
                            nrow = 1
                        )
                    )
                )
            }
        }
    }

    dev.off()
    message(fn)
}


if(!interactive()) {
    main()
}