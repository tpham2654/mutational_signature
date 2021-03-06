import copy
import re
import loadmutsigs
import os
import sys
from os.path import isfile
from pprint import pprint

import pyfaidx
from pandas import DataFrame, isnull, read_table
if sys.version_info[0] < 3: 
    from StringIO import StringIO
else:
    from io import StringIO

def compute_mutational_signature_enrichment(mutation_file_path,
                                                   reference_file_path,
                                                   signature_number=2,
                                                   upper_fasta=True,
                                                   chromosome_format='ID',
                                                   regions=None,
                                                   ids=None,
                                                   verbose=False):
    """
    Compute APOBEC mutational signature enrichment.
    Arguments:
        mutation_file_path (str | iterable): file_path(s) to mutation files
            (.vcf file | .vcf.gz file | .maf file)
        reference_file_path (str): file path to referece genome (.fasta file |
            .fa file)
        signature_number (int or list [signumber int,doreversecomplement bool]): Mutation signature from http://cancer.sanger.ac.uk/cosmic/signatures
        upper_fasta (bool): whether to read all .fasta file seqeunces as
            upper case
        chromosome_format (str): 'ID' | 'chrID'
        regions (dict):
        ids (None | iterable): count all variants if None or only variants
            with IDs present in the iterable (thus ignores all variants with IDs
            if an empty iterable)
        verbose (bool)
    Returns:
        DataFrame: (n_mutation + n_motifs counted, n_mutation_file, )
    """

    if regions is None:
        regions = {}

    # If only 1 file_path is passed, put it in a list
    if isinstance(mutation_file_path, str):
        mutation_file_path = [mutation_file_path]

    # Load reference genome (mutation files must use the same reference)
    if not isfile(reference_file_path):
        raise ValueError(
            '.fasta file {} doesn\'t exist.'.format(reference_file_path))

    reference_index_file_path = '{}.fai'.format(reference_file_path)
    if not isfile(reference_index_file_path):
        raise ValueError(
            '.fai file {} doesn\'t exist.'.format(reference_index_file_path))

    fasta = pyfaidx.Fasta(
        reference_file_path, sequence_always_upper=upper_fasta)

    if not fasta:
        raise ValueError('Loaded nothing from .fasta file.')

    print('Loaded .fasta file: {}.'.format(list(fasta.keys())))

    span = 20

    # Set up mutational signature and their weights shown in COSMIC figure
    
    if signature_number==-1:
        ss = {
        'TCA ==> TGA': 1,
        'TCA ==> TTA': 1,
        'TCT ==> TGT': 1,
        'TCT ==> TTT': 1,
        # Reverse complement
        'TGA ==> TCA': 1,
        'TGA ==> TAA': 1,
        'AGA ==> ACA': 1,
        'AGA ==> AAA': 1,
        }
    else:
	    if type(signature_number) is list:
	    	ss = loadmutsigs.makesigdict(signature_number[0],signature_number[1])
	    else:
	        ss = loadmutsigs.makesigdict(signature_number)
    
    # Identify what to count to compute enrichment
    signature_mutations, control_mutations, signature_b_motifs, control_b_motifs = _identify_what_to_count(
        ss)

    # Count
    samples = {}
    for i, fp in enumerate(mutation_file_path):

        # Get sample ID
        id_ = fp.split('/')[-1].split('.')[0]
        if id_ in samples:
            raise ValueError('{} duplicated'.format(id_))
        print('({}/{}) {} ...'.format(i + 1, len(mutation_file_path), id_))

        # Count
        samples[id_] = _count(
            fp,
            fasta,
            span,
            signature_mutations,
            control_mutations,
            signature_b_motifs,
            control_b_motifs,
            chromosome_format=chromosome_format,
            regions=regions,
            ids=ids,
            verbose=verbose)

    # Tabulate results
    df = DataFrame(samples)
    df.columns.name = 'Sample'

    df.ix['TCW'] = df.ix[['TCA', 'TCT']].sum()
    df.ix['TCW ==> TGW'] = df.ix[['TCA ==> TGA', 'TCT ==> TGT']].sum()
    df.ix['TCW ==> TTW'] = df.ix[['TCA ==> TTA', 'TCT ==> TTT']].sum()

    df.ix['WGA'] = df.ix[['AGA', 'TGA']].sum()
    df.ix['WGA ==> WCA'] = df.ix[['AGA ==> ACA', 'TGA ==> TCA']].sum()
    df.ix['WGA ==> WAA'] = df.ix[['AGA ==> AAA', 'TGA ==> TAA']].sum()

    n_signature_mutations = df.ix[signature_mutations.keys()].apply(
        lambda s: s * signature_mutations[s.name]['weight'], axis=1).sum()

    n_control_mutations = df.ix[control_mutations.keys()].apply(
        lambda s: s * control_mutations[s.name]['weight'], axis=1).sum()

    n_signature_b_mutations = df.ix[signature_b_motifs.keys()].apply(
        lambda s: s * signature_b_motifs[s.name]['weight'], axis=1).sum()

    n_control_b_mutations = df.ix[control_b_motifs.keys()].apply(
        lambda s: s * control_b_motifs[s.name]['weight'], axis=1).sum()

    mse = (n_signature_mutations / n_control_mutations) / (
        n_signature_b_mutations / n_control_b_mutations)
    if type(signature_number) is list:
        df.ix['Mutation Signature'] = signature_number[0]
    else:
        df.ix['Mutation Signature'] = signature_number
    df.ix['Mutational Signature Enrichment'] = mse.fillna(0)

    return df.sort_index()

def _identify_what_to_count(signature_mutations):
    """
    Identigy what mutations and/or motifs to count to compute mutational
        signature enrichment:
            signature_mutation_enrichment = [(n_signature_motif_mutations) /
            (n_changing_signature_motif_mutations)] /
            [(n_signature_motifs_in_context) /
            (n_changing_signature_motifs_in_context)]
    Arguments:
        signature_mutations (dict): {signature_mutation: weight, ...}
    Returns:
        dict:
        dict:
        dict:
        dict:
    """

    # Signature mutations
    s_mutations = {}

    # Control mutations
    c_mutations = {}

    # Signature before-motifs
    s_b_motifs = {}

    # Control before-motifs
    c_b_motifs = {}

    for m, w in signature_mutations.items():

        # Signature mutations

        # Get before & after motifs, which must be the same length
        m_split = m.split('==>')
        b_m, a_m = [m_.strip() for m_ in m_split]
        if len(b_m) != len(a_m):
            raise ValueError(
                'Before ({}) & after ({}) motifs differ in length.'.format(
                    b_m, a_m))

        s_mutations[m] = {
            'before': b_m,
            'after': a_m,
            'n': 0,
            'weight': w,
        }

        # Control mutations

        # Get changing-before & -after motifs
        differences = [i for i in range(len(b_m)) if b_m[i] != a_m[i]]
        c_s_i = min(differences)
        c_e_i = max(differences)
        c_b_m = b_m[c_s_i:c_e_i + 1]
        c_a_m = a_m[c_s_i:c_e_i + 1]

        k = '{} ==> {}'.format(c_b_m, c_a_m)
        if k in c_mutations:
            c_mutations[k]['weight'] = max(c_mutations[k]['weight'], w)
        else:
            c_mutations[k] = {
                'before': c_b_m,
                'after': c_a_m,
                'n': 0,
                'weight': w,
            }

        # Signature before-motifs

        k = b_m
        if k in s_b_motifs:
            s_b_motifs[k]['weight'] = max(s_b_motifs[k]['weight'], w)
        else:
            s_b_motifs[k] = {
                'n': 0,
                'weight': w,
            }

        # Control before-motifs

        k = c_b_m
        if k in c_b_motifs:
            c_b_motifs[k]['weight'] = max(c_b_motifs[k]['weight'], w)
        else:
            c_b_motifs[k] = {
                'n': 0,
                'weight': w,
            }

    print('s_mutations:')
    pprint(s_mutations)

    print('\nc_mutations:')
    pprint(c_mutations)

    print('\ns_b_motifs:')
    pprint(s_b_motifs)

    print('\nc_b_motifs:')
    pprint(c_b_motifs)

    return s_mutations, c_mutations, s_b_motifs, c_b_motifs

	
def _detect_maf_cols(mutation_file_path):
    expected_columns = ["Chromosome",
    "Start_Position",
    "dbSNP_RS",
    "Reference_Allele",
    "Tumor_Seq_Allele2"]
    df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1')
    df_cols = list(df.columns.values)
    df_unique_cols = set(df_cols)
    df_cols_in_expected = [c for c in df_unique_cols if c.lower() in [e_c.lower() for e_c in expected_columns]]
    df_cols_in_expected_lower = [c.lower() for c in df_cols_in_expected]
    if len(df_cols_in_expected)==len(expected_columns):
        expected_columns_lower = [c.lower() for c in expected_columns]
        df_cols_ordered = []
        for i in range(0, len(expected_columns_lower)):
            df_cols_ordered.append(df_cols_in_expected[df_cols_in_expected_lower.index(expected_columns_lower[i])])
        return df[df_cols_ordered]
    else:
        raise ValueError('Unknown formatted MAF: {}.'.format(mutation_file_path))

def _count(mutation_file_path,
           fasta,
           span,
           signature_mutations,
           control_mutations,
           signature_b_motifs,
           control_b_motifs,
           chromosome_format='ID',
           regions=None,
           ids=None,
           verbose=False):
    """
    Count motifs.
    Arguments:
        mutation_file_path (str | iterable): file_path(s) to mutation files
            (.vcf file | .vcf.gz file | .maf file)
        fasta (pyfaidx handle):
        span (int):
        signature_mutations (dict):
        control_mutations (dict):
        signature_b_motifs (dict):
        control_b_motifs (dict):
        chromosome_format (str): 'ID' | 'chrID'
        regions (dict):
        ids (None | iterable): count all variants if None or only variants
            with IDs present in the iterable (thus ignores all variants with IDs if
            an empty iterable)
        verbose (bool):
    Returns:
        dict:
    """

    if regions is None:
        regions = {}

    # Load mutation file
    if mutation_file_path.endswith(('.vcf', '.vcf.gz')):
        vcf_lines = []
        with open(mutation_file_path, 'r') as f:
            vcf_lines = [l for l in f if not l.startswith('##')]
        df = read_table(StringIO(os.linesep.join(vcf_lines)), 
            encoding='ISO-8859-1').iloc[:, [0, 1, 2, 3, 4]]
    elif mutation_file_path.endswith(('.maf','.maf.txt')):
        '''df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1').iloc[:, [4, 5, 13, 10, 12]]'''
        df = _detect_maf_cols(mutation_file_path)
    else:
        raise ValueError(
            'Unknown mutation_file_path: {}.'.format(mutation_file_path))

    # Get ready to count mutations &/| motifs
    s_mutations = copy.deepcopy(signature_mutations)
    c_mutations = copy.deepcopy(control_mutations)
    s_b_motifs = copy.deepcopy(signature_b_motifs)
    c_b_motifs = copy.deepcopy(control_b_motifs)

    # Initialize counts
    n_mutation_in_region = 0
    n_kept_mutations_in_region = 0
    n_mutation_analyzed = 0
    n_spanning_bases = 0

    # Use dict for faster ID look up
    if ids is not None:
        ids = {id_: None for id_ in ids}
	print df.empty
    for i, (chr_, pos, id_, ref, alt) in df.iterrows():
        #print "inside for loop"
        # Name chromosome
        chr_ = str(chr_)
        if chromosome_format == 'chrID':
            if not chr_.startswith('chr'):
                chr_ = 'chr{}'.format(chr_).replace('MT', 'M')
        elif chromosome_format == 'ID':
            if chr_.startswith('chr'):
                chr_ = chr_.replace('chr', '').replace('M', 'MT')
        else:
            raise ValueError(
                'Unknown chromosome_format: {}.'.format(chromosome_format))

        # Shift position
        pos = int(pos) - 1

        # Skip if variant is not in the specified regions
        if regions:

            r = regions.get(chr_)

            if not r or not any([s < pos < e for s, e in r]):
                if verbose:
                    print('\t{}:{} not in regions.'.format(chr_, pos))
                continue

        n_mutation_in_region += 1

        # Filter variant with IDs
        if ids is not None and not isnull(id_) and id_.startswith('rs'):
            if id_ in ids:
                if verbose:
                    print('\tKeep variant with ID {}.'.format(id_))
            else:
                if verbose:
                    print('\tSkip variant with ID {}.'.format(id_))
                continue

        n_kept_mutations_in_region += 1

        # Skip if there is no reference information
        if chr_ not in fasta.keys():
            print('\tChromosome {} not in .fasta file.'.format(chr_))
            #pprint(fasta.keys())
            continue

        # Skip if variant is not a SNP
        #print ref
        #print alt
        if isinstance(ref, basestring) and isinstance(alt, basestring):
            #has problems with  MAF made with R maftools.icgcSimpleMutationToMAF for samples like PBCA-DE_DO35598._nodup.tsv.maf because ref or alt is boolean (True), not string.
            if not (1 == len(ref) == len(alt)) or ref == '-' or alt == '-':
                if verbose:
                    print('\tSkip non-SNP variant {} ==> {}.'.format(ref, alt))
                continue

        if ref != fasta[chr_][pos].seq:
            #TCGA-04-1337-01 has problems with this line when you use 1000genomes FASTA
            print('\tReferences mismatch: {}:{} {} != ({}){}({}).'.format(
                chr_, pos, ref, *fasta[chr_][pos - 1:pos + 2].seq))
            continue

        # Analyze
        n_mutation_analyzed += 1

        # Check if this mutation matches any signature mutation
        for m, d in s_mutations.items():

            # Get before & after signature motifs
            b_m = d.get('before')
            a_m = d.get('after')

            # Get changing-before & -after motifs
            differences = [i for i in range(len(b_m)) if b_m[i] != a_m[i]]
            c_s_i = min(differences)
            c_e_i = max(differences)
            c_b_m = b_m[c_s_i:c_e_i + 1]
            c_a_m = a_m[c_s_i:c_e_i + 1]

            # Check if the chaning-before & -after motifs match REF & ALT
            # Check if the surrounding sequences are the same
            if c_b_m == ref:
                if c_a_m == alt:
                    if b_m == fasta[chr_][pos - c_s_i:
                                          pos + len(b_m) - c_e_i].seq:
                        d['n'] += 1

        # Check if this mutation matches any control mutation
        for m, d in c_mutations.items():

            # Check if the before & after motifs match REF & ALT
            if d.get('before') == ref and d.get('after') == alt:
                d['n'] += 1

        # Get mutation-spanning sequences and strip 'N's
        
        
        start_pos = pos - span
        end_pos = pos + span + 1
        if pos < span:
            start_pos = 0
            end_pos = 21
        #print str(chr_) + " " + str(pos) + " [" + str(start_pos) + " - " + str(end_pos) + "]"
        
        span_seq = fasta[chr_][start_pos:end_pos].seq
        #print span_seq
        if "N" in span_seq:
            span_seq = span_seq.strip('N')

        if re.findall('[^AaCcGgTt]', span_seq):
            print('\t{} (centered on {}:{}) contains unknown nucleotide.'.
                  format(span_seq, chr_, pos))

        n_spanning_bases += len(span_seq)

        # Count signature's changing-before motifs in the spanning sequences
        for m in s_b_motifs:
            s_b_motifs[m]['n'] += span_seq.count(m)

        # Count control's changing-before motifs in the spanning sequences
        for m in c_b_motifs:
            c_b_motifs[m]['n'] += span_seq.count(m)
	counts = {
        'N Entry in Mutation File': 1,
        'N Mutation in Region': n_mutation_in_region,
        'N Kept Mutation in Region': n_kept_mutations_in_region,
        'N Mutation Analyzed': n_mutation_analyzed,
        'N Spanning Bases': n_spanning_bases,
        }
    try:
        i
    except NameError:     
        i = 0
    counts = {
        'N Entry in Mutation File': i + 1,
        'N Mutation in Region': n_mutation_in_region,
        'N Kept Mutation in Region': n_kept_mutations_in_region,
        'N Mutation Analyzed': n_mutation_analyzed,
        'N Spanning Bases': n_spanning_bases,
    }
    pprint(counts)
    counts.update({m: d['n'] for m, d in s_mutations.items()})
    counts.update({m: d['n'] for m, d in c_mutations.items()})
    counts.update({m: d['n'] for m, d in s_b_motifs.items()})
    counts.update({m: d['n'] for m, d in c_b_motifs.items()})

    return counts