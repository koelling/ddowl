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
p <- add_argument(p, "--prefix", help="filename prefix")

phred2p = function(x) 10^(-x/10)
p2phred = function(x) -log10(x)*10

xx = c(0.0005, 0.003, 0.9, 1)
stopifnot(all.equal(phred2p(p2phred(xx)), xx))

alleles = c('paternal', 'maternal')

phred2pcorrect = function(x) 1-phred2p(x)

joint_snp_prop = function(snp_call_qualities, snp_call_types, expected_call) {
    other_call = alleles[alleles != expected_call]

    prod(c(
        phred2pcorrect(snp_call_qualities[snp_call_types == expected_call]),
        1-phred2pcorrect(snp_call_qualities[snp_call_types == other_call])
    ))
}

calc_read_p = function(expected_mutation_allele, mutation_call_type, mutation_call_quality, snp_call_qualities, snp_call_types) {
    expected_wt_allele = alleles[alleles != expected_mutation_allele]

    r = (
        #assuming the mutation call is correct
        (1-phred2p(mutation_call_quality)) *
        ifelse(
            mutation_call_type == 'wild-type',
            joint_snp_prop(snp_call_qualities, snp_call_types, expected_wt_allele),
            joint_snp_prop(snp_call_qualities, snp_call_types, expected_mutation_allele)
        )
        +
        #... or assuming the mutation call is wrong
        (phred2p(mutation_call_quality)) *
        ifelse(
            mutation_call_type == 'wild-type',
            joint_snp_prop(snp_call_qualities, snp_call_types, expected_mutation_allele),
            joint_snp_prop(snp_call_qualities, snp_call_types, expected_wt_allele)
        )
    )

#     if(length(r) > 1) {
#         str(list(expected_mutation_allele, mutation_call_type, mutation_call_quality, snp_call_qualities, snp_call_types, r))
#         stop()
#     }

    r
}

calc_read_likelihoods = function(reads) {
    reads %>%
            filter(!is.na(mutation_call_type) & mutation_call_type %in% c('wild-type', 'mutation')) %>%
            filter(!is.na(snp_call_type) & snp_call_type %in% c('paternal', 'maternal')) %>%
            group_by(FamilyID, sample, read, mutation_call_type, mutation_call_quality) %>%
            summarise(
                n_snps = n(),
                n_paternal = sum(snp_call_type == 'paternal'),
                n_maternal = sum(snp_call_type == 'maternal'),
                p_if_mut_mat = calc_read_p('maternal', mutation_call_type[1], mutation_call_quality[1], snp_call_quality, snp_call_type),
                p_if_mut_pat = calc_read_p('paternal', mutation_call_type[1], mutation_call_quality[1], snp_call_quality, snp_call_type),
                bf = p_if_mut_mat / p_if_mut_pat,
                logbf = log(p_if_mut_mat) - log(p_if_mut_pat),
                type = ifelse(p_if_mut_mat > 0.5, 'mat', ifelse(p_if_mut_pat > 0.5, 'pat', 'other'))
            )
}

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

    is = read_csv(sprintf('%s.informative_snps.csv', args$prefix)) %>%
        mutate(multi = grepl(',', paternal) | grepl(',', maternal))
    rmc = read_csv(sprintf('%s.read_mutation_calls.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample'))
    rsc = read_csv(sprintf('%s.read_snp_calls.csv', args$prefix)) %>%
        inner_join(families, by=c('FamilyID', 'sample'))
    rc = rsc
    #shuffle:
    #rmc = rmc %>% group_by(sample) %>% mutate(read = sample(read)) %>% ungroup

    fn = sprintf('%s.phasing.pdf', args$prefix)
    pdf(fn)
    for(family in unique(rmc$FamilyID)) {
        cat(family, ':\n')
        my_is = is %>% filter(FamilyID == family)
        my_rmc = rmc %>% filter(FamilyID == family)
        my_rsc = rsc %>% filter(FamilyID == family)
        my_rc = rc %>% filter(FamilyID == family)

        proband_rsc = rsc %>% filter(relationship == 'proband')
        proband_rc = proband_rsc %>%
            filter(snp_name %in% my_is$name[!my_is$multi]) %>% #FIXME: remove this
            select(FamilyID, sample, read, snp_name, snp_call_type, snp_call_quality) %>%
            left_join(rmc %>% select(FamilyID, sample, read, mutation_call_type, mutation_call_quality), by = c('FamilyID', 'sample', 'read'))
        print(table(proband_rc$mutation_call_type, proband_rc$snp_call_type))

        by_read = proband_rc %>%
            filter(!is.na(mutation_call_type)) %>%
            filter(!is.na(snp_call_type)) %>%
            group_by(FamilyID, sample, read, mutation_call_type, mutation_call_quality) %>% summarise(
                snps = list(snp_name),
                snp_calls = list(snp_call_type),
                snp_qualities = list(snp_call_quality),
                inconsistent = length(unique(snp_call_type)) != 1,
                n_snps = n(),
                n_paternal = sum(snp_call_type == 'paternal'),
                n_maternal = sum(snp_call_type == 'maternal'),
                q_paternal = sum(snp_call_quality[snp_call_type == 'paternal']),
                q_maternal = sum(snp_call_quality[snp_call_type == 'maternal']),
            )

        cat('Calculating foreground...\n')
        read_likelihood = calc_read_likelihoods(proband_rc)
        print(read_likelihood %>% head)

        cat('Reads with p(X) > 0.5:')
        print(table(read_likelihood$type))

        median_logbf = median(read_likelihood$logbf)
        cat('Median logBF for maternal > paternal: ', median_logbf, '\n')
        total_logbf = sum(read_likelihood$logbf)
        cat('Overall logBF for maternal > paternal: ', total_logbf, '\n')

        cat('Wilcoxon test of p(Mat) over p(Pat):\n')
        print(wilcox.test(read_likelihood$p_if_mut_mat, read_likelihood$p_if_mut_pat, paired = T))
        binom = binom.test(x = sum(read_likelihood$p_if_mut_mat > read_likelihood$p_if_mut_pat), n = nrow(read_likelihood))
        cat('Sign test of p(Mat) over p(Pat):\n')
        print(binom.test(x = sum(read_likelihood$p_if_mut_mat > read_likelihood$p_if_mut_pat), n = nrow(read_likelihood)))

        cat('Calculating background...\n')
        read_likelihood_bg = sapply(1:100, function(i) {
            cat(i, ' ')
            #note the  %>% mutate(read = sample(read)
            background_rc = proband_rsc %>%
                filter(snp_name %in% my_is$name[!my_is$multi]) %>% #FIXME: remove this
                select(FamilyID, sample, read, snp_name, snp_call_type, snp_call_quality) %>%
                left_join(rmc %>% select(FamilyID, sample, read, mutation_call_type, mutation_call_quality) %>% mutate(read = sample(read)), by = c('FamilyID', 'sample', 'read'))

            background_likelihood = calc_read_likelihoods(background_rc)

            c(
                bg_median_logbf = median(background_likelihood$logbf),
                bg_total_logbf = sum(background_likelihood$logbf),
                bg_p_wilcox = wilcox.test(background_likelihood$p_if_mut_mat, background_likelihood$p_if_mut_pat, paired = T)$p.value,
                bg_p_binom = binom.test(x = sum(background_likelihood$p_if_mut_mat > background_likelihood$p_if_mut_pat), n = nrow(background_likelihood))$p.value
            )
        }) %>%
        t %>%
        data.frame

        cat('.\n')

        cat('Emperical p: Median logBF for maternal > paternal: ', sum(abs(median_logbf) > abs(read_likelihood_bg$bg_median_logbf)) / (nrow(read_likelihood_bg)+1), '\n')
        cat('Emperical p: Overall logBF for maternal > paternal: ', sum(abs(total_logbf) > abs(read_likelihood_bg$bg_total_logbf)) / (nrow(read_likelihood_bg)+1), '\n')
        cat('Emperical p: Sign test of p(Mat) over p(Pat): ', sum(binom$p.value > read_likelihood_bg$bg_p_binom) / (nrow(read_likelihood_bg)+1), '\n')


        print(
            ggplot(my_rmc, aes(factor(1), fill=rel_sam)) + geom_bar(position = 'dodge') + coord_flip()
            + ggtitle(family, 'mutation site reads')
        )

        print(
            ggplot(my_rmc, aes(ifelse(is.na(mutation_call_type), 'NA', mutation_call_type), fill=rel_sam)) + geom_bar(position = 'dodge') + coord_flip()
            + xlab('mutation site call')
            + ggtitle(family, 'mutation site calls')
        )

        print(
            ggplot(my_rmc, aes(rel_sam, mutation_call_quality, fill=ifelse(is.na(mutation_call_type), 'NA', mutation_call_type))) + geom_violin(position = 'dodge') + scale_fill_discrete(name='')
            + ggtitle(family, 'mutation site call quality')
        )

        print(
            ggplot(my_rsc, aes(snp_name, fill=rel_sam)) + geom_bar(position = 'dodge') + coord_flip()
            + ggtitle(family, 'snp read count')
        )

        print(
            ggplot(my_rsc, aes(snp_name, fill=snp_call_type)) + geom_bar(position = 'dodge') + coord_flip() + facet_grid(~rel_sam)
            + ggtitle(family, 'snp calls')
        )

        print(
            ggplot(rsc, aes(snp_name, snp_call_quality, fill=snp_call_type)) + geom_violin(position = 'dodge') + scale_fill_discrete(name='')
             + coord_flip() + facet_grid(~rel_sam)
             + ggtitle(family, 'snp call quality')
        )

        print(
            ggplot(by_read, aes(n_snps)) + geom_bar() + facet_grid(~mutation_call_type) + ggtitle(family, 'n snps by mutation type')
        )

        print(
            ggplot(by_read, aes(inconsistent)) + geom_bar() + facet_grid(~mutation_call_type) + ggtitle(family, 'snp consistency by mutation type')
        )

        print(
            ggplot(by_read, aes(n_paternal, n_maternal, color = mutation_call_type))
            + geom_jitter()
            + geom_abline()
            + ggtitle(family, 'Number of SNP alleles by mutation allele')
        )

        print(
            ggplot(read_likelihood, aes(logbf)) +
            geom_histogram()
            + ggtitle(family, 'Per-read logBF for maternal > paternal')
        )

        print(
            ggplot(read_likelihood_bg, aes(bg_median_logbf))
            + geom_histogram()
            + geom_vline(xintercept = median_logbf, color = 'red')
            + ggtitle(family, 'foreground vs. background median logbf')
        )

        print(
            ggplot(read_likelihood_bg, aes(bg_total_logbf))
            + geom_histogram()
            + geom_vline(xintercept = total_logbf, color = 'red')
            + ggtitle(family, 'foreground vs. background total logbf')
        )


    }
    dev.off()
    message(fn)
}


if(!interactive()) {
    main()
}