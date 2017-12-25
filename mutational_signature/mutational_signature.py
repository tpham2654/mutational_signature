import copy
import re
from os.path import isfile
from pprint import pprint

import pyfaidx
from pandas import DataFrame, isnull, read_table

# TODO: refactor


def compute_apobec_mutational_signature_enrichment(mutation_file_path,
                                                   reference_file_path,
                                                   upper_fasta=True,
                                                   chromosome_format='ID',
                                                   regions={},
                                                   ids=None,
                                                   verbose=False):
    """
    Compute APOBEC mutational signature enrichment.
    Arguments:
        mutation_file_path (str | iterable): file_path(s) to mutation files
            (.vcf file | .vcf.gz file | .maf file)
        referece_file_path (str): file path to referece genome (.fasta file |
            .fa file)
        upper_fasta (bool): whether to read all .fasta file seqeunces as
            upper case
        chromosome_format (str): 'ID' | 'chrID'
        regions (dict):
        ids (None | iterable): count all variants if None or only variants
            with IDs present in the iterable (thus ignores all variants with IDs
            if an empty iterable)
        verbose (bool)
    Returns:
        DataFrame: (n_mutations + n_motifs counted, n_mutation_files)
    """

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

    # Identigy what to count to compute enrichment
    signature_mutations, control_mutations, signature_b_motifs, control_b_motifs = _identify_what_to_count(
        ss)

    # Count
    samples = {}
    for i, fp in enumerate(mutation_file_path):

        # TODO: generalize
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

    amse = (n_signature_mutations / n_control_mutations) / (
        n_signature_b_mutations / n_control_b_mutations)

    df.ix['APOBEC Mutational Signature Enrichment'] = amse.fillna(0)

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


def _count(mutation_file_path,
           fasta,
           span,
           signature_mutations,
           control_mutations,
           signature_b_motifs,
           control_b_motifs,
           chromosome_format='ID',
           regions={},
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

    # Load mutation file
    if mutation_file_path.endswith(('.vcf', '.vcf.gz')):
        df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1').iloc[:, [0, 1, 2, 3, 4]]

    elif mutation_file_path.endswith('.maf'):
        df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1').iloc[:, [4, 5, 13, 10, 12]]

    # Get ready to count mutations &/| motifs
    s_mutations = copy.deepcopy(signature_mutations)
    c_mutations = copy.deepcopy(control_mutations)
    s_b_motifs = copy.deepcopy(signature_b_motifs)
    c_b_motifs = copy.deepcopy(control_b_motifs)

    # Initialize counts
    n_mutations_in_region = 0
    n_kept_mutations_in_region = 0
    n_mutations_analyzed = 0
    n_spanning_bases = 0

    # Use dict for faster ID look up
    if ids is not None:
        ids = {id_: None for id_ in ids}

    for i, (chr_, pos, id_, ref, alt) in df.iterrows():

        # TODO: remove
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

        n_mutations_in_region += 1

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
            continue

        # Skip if variant is not a SNP
        if not (1 == len(ref) == len(alt)) or ref == '-' or alt == '-':
            if verbose:
                print('\tSkip non-SNP variant {} ==> {}.'.format(ref, alt))
            continue

        if ref != fasta[chr_][pos].seq:
            print('\tRefereces mismatch: {}:{} {} != ({}){}({}).'.format(
                chr_, pos, ref, *fasta[chr_][pos - 1:pos + 2].seq))
            continue

        # Analyze
        n_mutations_analyzed += 1

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
        span_seq = fasta[chr_][pos - span:pos + span + 1].seq.strip('N')

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
        'N Entries in Mutation File': i + 1,
        'N Mutations in Region': n_mutations_in_region,
        'N Kept Mutations in Region': n_kept_mutations_in_region,
        'N Mutations Analyzed': n_mutations_analyzed,
        'N Spanning Bases': n_spanning_bases,
    }
    pprint(counts)
    counts.update({m: d['n'] for m, d in s_mutations.items()})
    counts.update({m: d['n'] for m, d in c_mutations.items()})
    counts.update({m: d['n'] for m, d in s_b_motifs.items()})
    counts.update({m: d['n'] for m, d in c_b_motifs.items()})

    return counts
