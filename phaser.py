import os, sys
import collections
import pprint
import pandas as pd
import pysam

class Call:
    def __init__(self, call, quality = None, is_error = False):
        self.call = call
        self.quality = quality
        self.is_error = is_error
        self.is_indel = len(call) > 1


def get_call_for_pileup_read(pileup_read, ref = None):
    if pileup_read.alignment.mapping_quality < 20:
        return Call('_MAPQ', is_error = True)
    elif pileup_read.alignment.is_secondary or pileup_read.alignment.is_supplementary or pileup_read.alignment.is_qcfail:
        return Call('_FLAG', is_error = True)
    elif pileup_read.alignment.is_duplicate:
        return Call('_DUPE', is_error = True)
    elif pileup_read.indel > 0:
        quals = pileup_read.alignment.query_qualities[pileup_read.query_position:pileup_read.query_position+pileup_read.indel+1]
        return Call(
            pileup_read.alignment.query_sequence[pileup_read.query_position:pileup_read.query_position+pileup_read.indel+1],
            1.0 * sum(quals) / len(quals)
        )
    elif pileup_read.indel < 0:
        #print(ref, pileup_read.indel, len(ref), ref[0:-abs(pileup_read.indel)])
        #if abs(pileup_read.indel) < len(ref):
        #    return ref[0:-abs(pileup_read.indel)]
        #else:
        #    return '_DEL'
        #hacky way to handle deletions...
        return Call(
            '%s-%d' % (pileup_read.alignment.query_sequence[pileup_read.query_position], abs(pileup_read.indel)),
            sum(pileup_read.alignment.query_qualities[pileup_read.query_position:pileup_read.query_position+2]) / 2.0
        )
    elif pileup_read.is_del:
        return Call('_DEL', is_error = True)
    elif pileup_read.is_refskip:
        return Call('_SKIP', is_error = True)
    else:
        return Call(pileup_read.alignment.query_sequence[pileup_read.query_position], pileup_read.alignment.query_qualities[pileup_read.query_position])

def get_read_calls(samfile, chrom, pos1, ref = None, max_depth = 1e7, calculate_stats = False):
    read_calls = {}
    read_stats = {}

    for pileup_column in samfile.pileup(chrom, pos1-1, pos1, max_depth = max_depth, stepper = 'nofilter', truncate=True):
        print(chrom, pos1, '->', pileup_column.reference_name, pileup_column.reference_pos, ' - found ', pileup_column.nsegments, 'alignments')
        for pileup_read in pileup_column.pileups:
            #ignore secondary and supplementary alignments
            if pileup_read.alignment.is_secondary:
                continue
            if pileup_read.alignment.is_supplementary:
                continue

            assert not pileup_read.alignment.query_name in read_calls, 'encountered multiple alignments for single read?'
            read_calls[pileup_read.alignment.query_name] = get_call_for_pileup_read(pileup_read, ref)

            if calculate_stats:
                read_stats[pileup_read.alignment.query_name] = {
                    'length': pileup_read.alignment.infer_query_length(),
                    'mismatches': pileup_read.alignment.get_tag('NM'),
                    'mapping_quality': pileup_read.alignment.mapping_quality,
                    'mean_baseq': 1.0 * sum(pileup_read.alignment.query_qualities) / len(pileup_read.alignment.query_qualities)
                }

    if calculate_stats:
        return read_calls, read_stats
    else:
        return read_calls

"Get counts of how often each call was observed at each SNP"
def get_call_counts(samfile, snps):
    snp_call_counts = {}
    for snp in snps.itertuples():
        snp_call_counts[snp.name] = collections.Counter()
        for pileup_column in samfile.pileup(snp.CHROM, snp.POS-1, snp.POS, max_depth = 1e4, stepper = 'nofilter', truncate=True):
            for pileup_read in pileup_column.pileups:
                call = get_call_for_pileup_read(pileup_read) #, snp.REF
                snp_call_counts[snp.name][call.call] += 1
    return snp_call_counts

def get_allele_type(allele, snp):
    if allele in snp.paternal and allele in snp.maternal:
        return 'shared'
    elif allele in snp.maternal:
        return 'maternal'
    elif allele in snp.paternal:
        return 'paternal'
    else:
        return None

def get_mutation_allele_type(allele, mutation):
    if allele == mutation['REF_processed']:
        return 'wild-type'
    elif allele == mutation['ALT_processed']:
        return 'mutation'
    else:
        return None

def process_family(fam, fam_rows, bam_mask, snps_mask, subsample = None):
    print()
    print(fam, 'STARTING')

    assert 'proband' in fam_rows['relationship'].values, 'need at least one proband'

    snps_fn = snps_mask % fam
    if not os.path.isfile(snps_fn):
        raise Exception('%s: %s missing!' % (fam, snps_fn))

    snps = pd.read_csv(snps_fn, sep='\t', dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str})
    snps.rename(columns = { '#CHROM': 'CHROM' }, inplace=True)
    snps.loc[~snps['CHROM'].str.startswith('chr'), 'CHROM'] = ['chr' + c for c in snps.loc[~snps['CHROM'].str.startswith('chr'), 'CHROM']]
    snps['name'] = ['%s_%d' % (snp.CHROM, snp.POS) for snp in snps.itertuples()]
    #the "calls" we actually get back from the pileup are just a single base, so if we have a deletion
    #such as GA>G, what we actually see if a G with indel == 0 or indel == 1.
    #thus, we need to adjust the REF/ALT we actually expect to see
    #this is stupid, what we should really do is to process the pileups in a smarter way...
    snps['REF_processed'] = [snp.REF[0] if len(snp.REF) > 1 else snp.REF for snp in snps.itertuples()]
    snps['ALT_processed'] = ['%s-%d' % (snp.REF[0], len(snp.REF) - len(snp.ALT)) if len(snp.REF) > 1 else snp.ALT for snp in snps.itertuples()]
    print(snps)

    mutation = snps[snps.ID == 'mutation']
    assert len(mutation) == 1, 'only one mutation allowed'
    mutation = mutation.iloc[0]

    background_snps_list = []
    for offset in [10, 50, 100]:
        for sign in [1, -1]:
            background_snps_list.append(
                pd.DataFrame([
                    {
                        'CHROM': snp.CHROM,
                        'POS': snp.POS + sign * offset,
                        'name': '{}_{}'.format(snp.name, sign * offset)
                    } for snp in snps.itertuples()
                ])
            )

    background_snps = pd.concat(background_snps_list)
    background_snps = background_snps[background_snps.POS > 0]
    background_snps.drop_duplicates(inplace = True)
    background_snps.set_index(['name'], inplace = True)

    print(fam, 'Using', len(background_snps), 'background SNPs:')
    print(background_snps)

    #get allele counts for all SNPs
    snp_sample_counts = collections.defaultdict(dict)
    for sample in fam_rows.itertuples():
        fn = bam_mask % sample.BC
        print(fam, 'Getting allele counts for', sample.BC)
        with pysam.AlignmentFile(fn, "rb") as samfile:
            snp_sample_counts[sample.BC] = get_call_counts(samfile, snps)
    print(fam, 'Checked', len(snp_sample_counts), 'SNPs.')

    #make a dataframe, the stupid way
    rows = []
    for sample, x in snp_sample_counts.items():
        for snp_name, xx in x.items():
            for snp_allele, count in xx.items():
                rows.append({
                    'count': count,
                    'snp_name': snp_name,
                    'snp_allele': snp_allele,
                    'sample': sample
                })
    all_allele_counts = pd.DataFrame(rows).set_index(['sample', 'snp_name', 'snp_allele'])
    print(all_allele_counts.head())

    if 'mother' in fam_rows['relationship'].values and 'father' in fam_rows['relationship'].values:
        assert (fam_rows['relationship'] == 'proband').sum() == 1, 'only one proband sample allowed'

        print('Found parents, finding informative SNPs...')
        informative_snp_dicts = []
        for snp in snps.itertuples():
            #we don't want the mutation
            if snp.ID == 'mutation':
                continue

            gt_calls = {}
            for sample in fam_rows.itertuples():
                relationship = sample.relationship
                sample_counts = snp_sample_counts[sample.BC]
                snp_counts = sample_counts[snp.name]

                #only consider REF + ALT here, ignore the rest
                total_counts = snp_counts[snp.REF_processed] + snp_counts[snp.ALT_processed]

                if snp_counts[snp.REF_processed] > 0.25 * total_counts and snp_counts[snp.ALT_processed] > 0.25 * total_counts:
                    gt_calls[relationship] = 'het'
                elif snp_counts[snp.REF_processed] > 0.6 * total_counts:
                    gt_calls[relationship] = snp.REF_processed
                elif snp_counts[snp.ALT_processed] > 0.6 * total_counts:
                    gt_calls[relationship] = snp.ALT_processed
                else:
                    gt_calls[relationship] = None

                print(snp.name, relationship, snp_counts[snp.REF_processed], snp_counts[snp.ALT_processed], total_counts, gt_calls[relationship])

            print(snp.name, gt_calls)

            if 'mother' in gt_calls and 'father' in gt_calls and 'proband' in gt_calls \
            and gt_calls['proband'] == 'het' \
            and gt_calls['mother'] != gt_calls['father'] \
            and gt_calls['mother'] is not None \
            and gt_calls['father'] is not None:
                snp_dict = snp._asdict()
                snp_dict['maternal'] = [gt_calls['mother']] if gt_calls['mother'] != 'het' else [snp.REF_processed, snp.ALT_processed]
                snp_dict['paternal'] = [gt_calls['father']] if gt_calls['father'] != 'het' else [snp.REF_processed, snp.ALT_processed]
                informative_snp_dicts.append(snp_dict)

        informative_snps = pd.DataFrame(informative_snp_dicts)
        if len(informative_snps) == 0:
            print(fam, 'No informative SNPs')
            return None, None, None
        informative_snps = informative_snps[list(snps.columns) + ['maternal', 'paternal']] #remove Index, ensure right order
    else:
        print(fam, 'Assuming that all non-mutation SNPs are informative')
        informative_snps = snps.loc[snps.ID != 'mutation'].copy()
        informative_snps['maternal'] = [[] for _ in range(len(informative_snps))]
        informative_snps['paternal'] = [[] for _ in range(len(informative_snps))]

    print(fam, 'Found', len(informative_snps), 'informative SNPs:')
    print(informative_snps)

    sample_informative_allele_counts = {}
    sample_mutation_read_stats = {}
    sample_read_mutation_calls = {}
    sample_read_snp_calls = {}
    sample_read_background_calls = {}

    phase_evidence = []
    for sample in fam_rows.itertuples(): #include all -- [rows['relationship'] == 'proband']
        fn = bam_mask % sample.BC
        print(fn)

        #get calls for each read at the mutation and each informative snp
        informative_calls = {}
        sample_read_background_calls_rows = []
        with pysam.AlignmentFile(fn, "rb") as samfile:
            mutation_calls, mutation_read_stats = get_read_calls(samfile, mutation['CHROM'], mutation['POS'], calculate_stats = True)

            #get calls for background snps (should all be ref and high q)
            for snp in background_snps.itertuples():
                for read, background_call in get_read_calls(samfile, snp.CHROM, snp.POS).items():
                    sample_read_background_calls_rows.append({
                        'read': read,
                        'snp_name': snp.Index, #name = Index
                        'snp_call': background_call.call,
                        'snp_call_quality': background_call.quality,
                    })

            #get calls for actual informative snps
            for snp in informative_snps.itertuples():
                informative_calls[snp.name] = get_read_calls(samfile, snp.CHROM, snp.POS)

        sample_read_background_calls[sample.BC] = pd.DataFrame(sample_read_background_calls_rows) \
            .set_index(['read']) \
            .join(background_snps, on = 'snp_name')

        #should we subsample reads? just subsample on mutation calls, we will simply ignore the reads for which we don't have a mutation call
        if subsample is not None:
            n_subsample = ceil(subsample * len(mutation_calls))
            print('subsampling', len(mutation_calls), 'mutation_calls to', n_subsample, '(', subsample, ')')
            mutation_calls = dict(random.sample(mutation_calls.items(), k = n_subsample))
            mutation_read_stats = dict(random.sample(mutation_read_stats.items(), k = n_subsample))
            assert len(mutation_calls) == n_subsample

        #process mutation read stats
        sample_mutation_read_stats[sample.BC] = pd.DataFrame(list(mutation_read_stats.values()))

        #summarize basecounts for mutation and each SNP
        sample_read_mutation_calls_rows = []
        site_basecounts = collections.Counter()
        for read, mutation_base in mutation_calls.items():
            simplified_call = mutation_base.call
            if mutation_base.is_error:
                simplified_call = '_FILTERED'
            elif mutation_base.is_indel and not mutation_base.call in [mutation['REF_processed'], mutation['ALT_processed']]:
                simplified_call = '_OTHER_INDEL'
            site_basecounts[(mutation['name'],simplified_call)] += 1

            sample_read_mutation_calls_rows.append({
                'read': read,
                'mutation_call': mutation_base.call,
                'mutation_call_type': get_mutation_allele_type(mutation_base.call, mutation = mutation),
                'mutation_call_simplified': simplified_call,
                'mutation_call_quality': mutation_base.quality,
            })
        sample_read_mutation_calls[sample.BC] = pd.DataFrame(sample_read_mutation_calls_rows) \
            .set_index(['read'])

        for snp in informative_snps.itertuples():
            snp_calls = informative_calls[snp.name]
            for snp_base in snp_calls.values():
                simplified_call = snp_base.call
                if snp_base.is_error:
                    simplified_call = '_FILTERED'
                elif snp_base.is_indel and not snp_base.call in [snp.REF_processed, snp.ALT_processed]:
                    simplified_call = '_OTHER_INDEL'
                site_basecounts[(snp.name,simplified_call)] += 1
        sample_informative_allele_counts[sample.BC] = pd.DataFrame(
            [{'count': x, 'snp_name': name, 'snp_allele': allele} for ((name, allele), x) in site_basecounts.items()],
        ).set_index(['snp_name', 'snp_allele'])

        #count haplotypes for each informative SNP
        sample_read_snp_calls_rows = []
        mutation_allele_counts = {}

        for snp in informative_snps.itertuples():
            snp_calls = informative_calls[snp.name]
            mutation_allele_counts[snp.name] = collections.defaultdict(collections.Counter)

            for read, mutation_base in mutation_calls.items():
                simplified_mutatation_call = mutation_base.call

                #only allow REF/ALT for mutation
                if not mutation_base.call in [mutation['REF_processed'], mutation['ALT_processed']]:
                    simplified_mutatation_call = '_OTHER'

                #did we cover the SNP? if so, get the call
                simplified_snp_call = '_NONE'
                if read in snp_calls:
                    snp_base = snp_calls[read]

                    #only allow REF/ALT for SNP
                    simplified_snp_call = snp_base.call
                    if not snp_base.call in [snp.REF_processed, snp.ALT_processed]:
                        simplified_snp_call = '_OTHER'

                    #note each call separately for later testing
                    sample_read_snp_calls_rows.append({
                        'read': read,
                        'snp_name': snp.name,
                        'snp_call': snp_base.call,
                        'snp_call_type': get_allele_type(snp_base.call, snp = snp),
                        'snp_call_quality': snp_base.quality,
                    })

                #also build up a dict of dicts with counts for each individual snp for simple heatmaps
                mutation_allele_counts[snp.name][simplified_mutatation_call][simplified_snp_call] += 1

        sample_read_snp_calls[sample.BC] = pd.DataFrame(sample_read_snp_calls_rows) \
            .set_index(['read', 'snp_name'])

        #post-process counts into dataframe
        for snp in informative_snps.itertuples():
            pair_counts = mutation_allele_counts[snp.name]

            df = pd.DataFrame.from_dict(pair_counts)
            df.reset_index(inplace = True)
            df.fillna(0, inplace=True)
            df.rename(columns={'index': 'snp_allele'}, inplace=True)
            df = df.melt(id_vars=['snp_allele'], var_name='mutation_allele', value_name='count').astype({'count': int})

            #add annotation about paternal/maternal allele
            df['snp_allele_type'] = df['snp_allele'].apply(get_allele_type, snp = snp)
            #add annotation about mutation/wild-type
            df['mutation_allele_type'] = df['mutation_allele'].apply(get_mutation_allele_type, mutation = mutation)

            #add general info
            df['sample'] = sample.BC
            df['snp_name'] = snp.name

            df = df[['sample', 'snp_name', 'snp_allele', 'mutation_allele', 'count', 'snp_allele_type', 'mutation_allele_type']] #, 'fraction'
            df.set_index(['sample', 'snp_name', 'snp_allele', 'mutation_allele'], inplace=True)
            phase_evidence.append(df)

        print(fam, 'processed sample:', sample)

    print(fam, 'DONE')

    return (
        informative_snps,
        all_allele_counts,
        pd.concat(sample_mutation_read_stats, names=['sample']) if len(sample_mutation_read_stats) > 0 else None,
        pd.concat(sample_informative_allele_counts, names=['sample']) if len(sample_informative_allele_counts) > 0 else None,
        pd.concat(sample_read_background_calls, names=['sample']) if len(sample_read_background_calls_rows) > 0 else None,
        pd.concat(sample_read_mutation_calls, names=['sample']) if len(sample_read_mutation_calls_rows) > 0 else None,
        pd.concat(sample_read_snp_calls, names=['sample']) if len(sample_read_snp_calls_rows) > 0 else None,
        pd.concat(phase_evidence) if len(phase_evidence) > 0 else None,
    )

def write_csv(d, fn):
    d.to_csv(fn)
    print('Wrote', len(d), 'rows:')
    print(fn)

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description = "ddOWL mutatation allele phasing and plotting tools, v0.1 - Nils Koelling",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )

    #specify parameters
    parser.add_argument("-v", "--version", help="print version and exit", action="store_true")
    parser.add_argument("--first", help="only first family", action="store_true")
    parser.add_argument("--debug", help="enable debug output", action="store_true")
    parser.add_argument("--family", help="specifiv family to run")
    parser.add_argument("FAMILIES", help="CSV file with family info")
    parser.add_argument("SNPS", help="mask for SNP CSV filenames (%s replaced by family ID)")
    parser.add_argument("BAMS", help="mask for BAM filenames (%s replaced by sample ID)")
    args = parser.parse_args()

    families = pd.read_csv(args.FAMILIES)
    assert len(families) > 0, 'no families'
    print(families.head())

    if args.family:
        print('Only analysing family {}'.format(args.family))
        families = families[families['FamilyID'] == args.family]
        assert len(families) > 0, 'no families after filtering'
        print(families.head())

        fn_prefix = '{}.{}'.format(args.FAMILIES, args.family)
    else:
        fn_prefix = args.FAMILIES

    if args.debug:
        fn_prefix = '{}.DEBUG'.format(fn_prefix)
    if args.first:
        fn_prefix = '{}.FIRST'.format(fn_prefix)

    informative_snp_dict = {}
    read_background_calls_dict = {}
    read_mutation_calls_dict = {}
    read_snp_calls_dict = {}
    phase_evidence_dict = {}
    read_stat_dict = {}
    allele_count_dict = {}
    informative_allele_count_dict = {}
    for family, group in families.groupby('FamilyID'):
        try:
            (
                fam_informative_snps,
                fam_allele_counts,
                fam_read_stats,
                fam_informative_allele_counts,
                fam_read_background_calls,
                fam_read_mutation_calls,
                fam_read_snp_calls,
                fam_phase_evidence
            ) = process_family(
                family,
                group,
                bam_mask = args.BAMS,
                snps_mask = args.SNPS,
            )

            if fam_informative_snps is not None:
                informative_snp_dict[family] = fam_informative_snps
            if fam_allele_counts is not None:
                allele_count_dict[family] = fam_allele_counts
            if fam_read_stats is not None:
                read_stat_dict[family] = fam_read_stats
            if fam_informative_allele_counts is not None:
                informative_allele_count_dict[family] = fam_informative_allele_counts
            if fam_read_background_calls is not None:
                read_background_calls_dict[family] = fam_read_background_calls
            if fam_read_mutation_calls is not None:
                read_mutation_calls_dict[family] = fam_read_mutation_calls
            if fam_read_snp_calls is not None:
                read_snp_calls_dict[family] = fam_read_snp_calls
            if fam_phase_evidence is not None:
                phase_evidence_dict[family] = fam_phase_evidence
        except Exception as e:
            #print(e)
            raise

        if args.first:
            print('--first is set: Stopping after first family!')
            break

    informative_snps = pd.concat(informative_snp_dict, names=['FamilyID'])
    write_csv(informative_snps, '%s.informative_snps.csv' % (fn_prefix))

    allele_counts = pd.concat(allele_count_dict, names=['FamilyID'])
    write_csv(allele_counts, '%s.allele_counts.csv' % (fn_prefix))

    read_stats = pd.concat(read_stat_dict, names=['FamilyID'])
    write_csv(read_stats, '%s.read_stats.csv' % (fn_prefix))

    informative_allele_counts = pd.concat(informative_allele_count_dict, names=['FamilyID'])
    write_csv(informative_allele_counts, '%s.informative_allele_counts.csv' % (fn_prefix))

    write_csv(pd.concat(read_background_calls_dict, names=['FamilyID']), '%s.read_background_calls.csv' % (fn_prefix))
    write_csv(pd.concat(read_mutation_calls_dict, names=['FamilyID']), '%s.read_mutation_calls.csv' % (fn_prefix))
    write_csv(pd.concat(read_snp_calls_dict, names=['FamilyID']), '%s.read_snp_calls.csv' % (fn_prefix))

    phase_evidence = pd.concat(phase_evidence_dict, names=['FamilyID'])
    write_csv(phase_evidence, '%s.phase_evidence.csv' % (fn_prefix))

if __name__ == '__main__':
    sys.exit(main())
