import copy
import re
from os.path import isfile
from pprint import pprint

import pyfaidx
from pandas import DataFrame, isnull, read_table


def get_apobec_mutational_signature_enrichment(mutation_file_path,
                                               reference_file_path,
                                               upper_fasta=True,
                                               chromosome_format='ID',
                                               regions={},
                                               ids=None,
                                               verbose=False):
    """
    Compute APOBEC mutational signature enrichment.
    :param mutation_file_path: str or iterable; of file_path(s) to mutations
    (.VCF | .VCF.GZ | .MAF)
    :param referece_file_path: str; file_path to referece genome (.FASTA | .FA)
    :param upper_fasta: bool; read all .FASTA seqeunces as uppe case or not
    :param chromosome_format: str; 'ID' | 'chrID'
    :param regions: dict
    :param ids: None | iterable; count all variants if None or only variants
        with IDs present in the iterable (thus ignores all variants with IDs if
        an empty iterable)
    :param verbose: bool
    :return: DataFrame; (n_mutations + n_motifs counted, n_mutation_files)
    """

    # If only 1 file_path is passed, put it in a list
    if isinstance(mutation_file_path, str):
        mutation_file_path = [mutation_file_path]

    # Load reference genome (mutation files must use the same reference)
    reference_index_file_path = '{}.fai'.format(reference_file_path)

    if not isfile(reference_index_file_path):
        raise ValueError('.FAI index {} doens\'t exist.'.format(
            reference_index_file_path))

    fasta = pyfaidx.Fasta(
        reference_file_path, sequence_always_upper=upper_fasta)

    if not fasta:
        raise ValueError('Loaded nothing from the reference genome.')

    if verbose:
        print('Loaded reference genome: {}.'.format(list(fasta.keys())))

    span = 20

    # Set up mutational signature and their weights shown in COSMIC
    ss = {
        'tCa ==> tGa': 0.3,
        'tCa ==> tTa': 0.5,
        'tCt ==> tGt': 0.4,
        'tCt ==> tTt': 0.3,
        # 'tCc ==> tGc': 0.1,
        # 'tCc ==> tTc': 0.1,
        # Reverse complement
        'tGa ==> tCa': 0.3,
        'tGa ==> tAa': 0.5,
        'aGa ==> aCa': 0.4,
        'aGa ==> aAa': 0.3,
        'gGa ==> gCa': 0.1,
        'gGa ==> gAa': 0.1,
    }

    # Identigy what to count to compute enrichment
    signature_mutations, control_mutations, signature_b_motifs, control_b_motifs = _identify_what_to_count(
        ss, verbose=verbose)

    # Count
    samples = {}
    for i, fp in enumerate(mutation_file_path):

        # Get sample ID
        id_ = fp.split('/')[-1].split('.')[0]
        if id_ in samples:
            raise ValueError('{} duplicated'.format(id_))
        print('({}/{}) {} ...'.format(i + 1, len(mutation_file_path), id_))

        # Count
        samples[id_] = count(
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


def _identify_what_to_count(signature_mutations, verbose=False):
    """
    Identigy what mutations and/or motifs to count to compute mutational
    signature enrichment:
        signature_mutation_enrichment = [(n_signature_motif_mutations) /
        (n_changing_signature_motif_mutations)] /
        [(n_signature_motifs_in_context) /
        (n_changing_signature_motifs_in_context)]
    :param signature_mutations: dict; {signature_mutation: weight, ...}
    :param verbose: bool
    :return: dict & dict & dict & dict
    """

    # Signature mutations
    s_mutations = {}
    for m, w in signature_mutations.items():

        # Get before & after motifs, which must be the same length
        m_split = m.split('==>')
        b_m, a_m = [m_.strip() for m_ in m_split]
        if len(b_m) != len(a_m):
            raise ValueError(
                'Before ({}) & after ({}) motifs differ in length.'.format(
                    b_m, a_m))

        s_mutations[m] = {
            'before': b_m.upper(),  # Before motif
            'after': a_m.upper(),  # After motif
            'n': 0,  # Mutation count
            'change_start_i': min([
                m.start() for m in re.finditer('[A-Z]+', b_m)
            ]),  # Chaning-motif start index
            'change_end_i': max([m.end() for m in re.finditer('[A-Z]+', b_m)
                                 ]),  # Changing-motif end index
            'weight': w
        }
    if verbose:
        print('s_mutations:')
        pprint(s_mutations)

    # Control mutations
    c_mutations = {}
    for d in s_mutations.values():

        # Get before & after motifs
        b_m = d.get('before')
        a_m = d.get('after')

        # Get changing-before & -after motifs
        c_s_i, c_e_i = d.get('change_start_i'), d.get('change_end_i')
        c_b_m = b_m[c_s_i:c_e_i]
        c_a_m = a_m[c_s_i:c_e_i]

        k = '{} ==> {}'.format(c_b_m, c_a_m)
        if k in c_mutations:
            c_mutations[k]['weight'] = max(c_mutations[k]['weight'],
                                           d.get('weight'))
        else:
            c_mutations[k] = {
                'before': c_b_m,
                'after': c_a_m,
                'n': 0,
                'weight': d.get('weight'),
            }
    if verbose:
        print('\nc_mutations:')
        pprint(c_mutations)

    # Signature before-motifs
    s_b_motifs = {}
    for d in s_mutations.values():
        k = d.get('before').lower()
        if k in s_b_motifs:
            s_b_motifs[k]['weight'] = max(s_b_motifs[k]['weight'],
                                          d.get('weight'))
        else:
            s_b_motifs[k] = {'n': 0, 'weight': d.get('weight')}
    if verbose:
        print('\ns_b_motifs:')
        pprint(s_b_motifs)

    # Control before-motifs
    c_b_motifs = {}
    for d in c_mutations.values():
        k = d.get('before').lower()
        if k in c_b_motifs:
            c_b_motifs[k]['weight'] = max(c_b_motifs[k]['weight'],
                                          d.get('weight'))
        else:
            c_b_motifs[k] = {'n': 0, 'weight': d.get('weight')}
    if verbose:
        print('\nc_b_motifs:')
        pprint(c_b_motifs)

    return s_mutations, c_mutations, s_b_motifs, c_b_motifs


def count(mutation_file_path,
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
    Count.
    :param mutation_file_path: str; file_path to mutations (.VCF | .VCF.GZ |
    .MAF)
    :param fasta: pyfaidx handle
    :param span: int
    :signature_mutations: dict
    :control_mutations: dict
    :signature_b_motifs: dict
    :control_b_motifs: dict
    :param chromosome_format: str; 'ID' | 'chrID'
    :param regions: dict
    :param ids: None | iterable; count all variants if None or only variants
        with IDs present in the iterable (thus ignores all variants with IDs if
        an empty iterable)
    :param verbose: bool
    :return: dict
    """

    # Load mutation file
    if mutation_file_path.endswith('.vcf') or mutation_file_path.endswith(
            '.vcf.gz'):
        df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1').iloc[:, [0, 1, 2, 3, 4]]

    elif mutation_file_path.endswith('.maf'):
        df = read_table(
            mutation_file_path, comment='#',
            encoding='ISO-8859-1').iloc[:, [4, 5, 13, 10, 12]]

    # Get ready to count mutations and/or motifs
    s_mutations = copy.deepcopy(signature_mutations)
    c_mutations = copy.deepcopy(control_mutations)
    s_b_motifs = copy.deepcopy(signature_b_motifs)
    c_b_motifs = copy.deepcopy(control_b_motifs)

    # Evaluate each row
    n_mutations_in_region = 0
    n_kept_mutations_in_region = 0
    n_mutations_analyzed = 0
    n_spanning_bases = 0

    # Use dict for faster ID look up
    if ids is not None:
        ids = {id_: None for id_ in ids}

    for i, (chr_, pos, id_, ref, alt) in df.iterrows():

        # TODO: Remove
        # Name chromosome
        chr_ = str(chr_)
        if chromosome_format == 'chrID':
            if not chr_.startswith('chr'):
                chr_ = 'chr{}'.format(chr_).replace('MT', 'M')
        elif chromosome_format == 'ID':
            if chr_.startswith('chr'):
                chr_ = chr_.replace('chr', '').replace('M', 'MT')
        else:
            raise ValueError('Unknown chromosome_format {}.'.format(
                chromosome_format))

        # TODO: Rationalize
        # Shift position
        pos = int(pos) - 1

        if regions:
            # Skip if variant is not in the specified regions
            r = regions.get(chr_)
            if not r or not any([s < pos < e for s, e in r]):
                if verbose:
                    print('\t{}:{} not in regions.'.format(chr_, pos))
                continue

        n_mutations_in_region += 1

        # Filter variant with IDs
        if ids is not None and not isnull(id_) and id_.startswith('rs'):
            if len(ids) and id_ in ids:
                if verbose:
                    print('\tKeep variant with ID {}.'.format(id_))
            else:
                if verbose:
                    print('\tSkip variant with ID {}.'.format(id_))
                continue

        n_kept_mutations_in_region += 1

        if chr_ not in fasta.keys():
            # Skip if there is no reference information
            if verbose:
                print('\tChromosome {} not in .FASTA.'.format(chr_))
            continue

        if not (1 == len(ref) == len(alt)) or ref == '-' or alt == '-':
            # Skip if variant is not a SNP
            if verbose:
                print('\tSkip non-SNP variant {} ==> {}.'.format(ref, alt))
            continue

        if ref != fasta[chr_][pos].seq:
            # if verbose:
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
            c_s_i, c_e_i = d.get('change_start_i'), d.get('change_end_i')
            c_b_m = b_m[c_s_i:c_e_i]
            c_a_m = a_m[c_s_i:c_e_i]

            # Check if the chaning-before motif matches the ref
            if c_b_m == ref:

                # Check if the surrounding sequences are the same
                if b_m == fasta[chr_][pos - c_s_i:pos + len(b_m) - c_e_i +
                                      1].seq:

                    # Check if the changing-after motif matches the alt
                    if c_a_m == alt:

                        # Matched a signature mutation, so increment count
                        d['n'] += 1

        # Check if this mutation matches any control mutation
        for m, d in c_mutations.items():

            # Get before & after control motifs
            b_m = d.get('before')
            a_m = d.get('after')

            # Check if the before motif matches the ref
            if b_m == ref:

                # Check if the after motif matches the alt
                if a_m == alt:

                    # Matched a control mutation, so increment count
                    d['n'] += 1

        # Get mutation-spanning sequences and strip 'N's
        span_seq = fasta[chr_][pos - span:pos + span + 1].seq.strip('N')
        if re.findall('[^ACGT]', span_seq):
            if verbose:
                print('\t{} (centered on {}:{}) contains unknown nucleotide.'.
                      format(span_seq, chr_, pos))
            continue

        n_spanning_bases += len(span_seq)

        # Count signature's changing-before motifs in the spanning sequences
        for m in s_b_motifs:
            s_b_motifs[m]['n'] += span_seq.count(m.upper())
        # Count control's changing-before motifs in the spanning sequences
        for m in c_b_motifs:
            c_b_motifs[m]['n'] += span_seq.count(m.upper())

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
