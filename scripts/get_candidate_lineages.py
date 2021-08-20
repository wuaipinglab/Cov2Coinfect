import os
from Bio.Seq import Seq
import pandas as pd

DIRPATH = '/SSD/yexiao/co_infection/'
genome_path = DIRPATH + 'data/SARS_CoV_2.csv'
lineage_10_path = DIRPATH + 'data/lineage_10.txt'
lineage_75_path = DIRPATH + 'data/lineage_75.txt'
ngs_dir = DIRPATH + 'project/ngs/'
candidate_dir = DIRPATH + 'project/candidate/'

if not os.path.exists(candidate_dir):
    os.makedirs(candidate_dir)


def read_dataframe(df_fc, tp):
    info = {}
    for row in df_fc.index:
        p = df_fc.loc[row, 'Reference Position']
        l = df_fc.loc[row, 'Length']
        a = df_fc.loc[row, 'Allele']
        f = df_fc.loc[row, 'Frequency']
        ps = []
        for n in range(0, l):
            ps.append(str(p + n))
        if tp == 'mutation':
            info[(','.join(ps))] = [a, f]
        else:
            info[(','.join(ps))] = [a * l, f]
    return info


def get_peptide_positions(tp):
    peptide_positions_dict_fc = {}
    for n in range(0, len(sites)):
        site = sites[n]
        if tp == 'mutation':
            base = mutations[m][0][n]
        else:
            base = '-'
        if site in genome['genomePos'].values:
            p_ps = genome[genome['genomePos'] == site]['peptidePos'].values
            for p_p in p_ps:
                peptide_positions_dict_fc.setdefault(p_p, []).append(str(site) + '/' + base)
    return peptide_positions_dict_fc


def generate_mutation_info(m_or_d, fc_aa_mut, tp):
    aa_ref = genome[genome['peptidePos'] == p_p]['aa'].values[0]
    if fc_aa_mut != aa_ref:
        if fc_aa_mut == '*':
            fc_aa_mut = 'stop'
        product = genome[genome['peptidePos'] == p_p]['product'].values[0]
        aa_position = str(genome[genome['peptidePos'] == p_p]['aaPos'].values[0])
        info = product + '_' + aa_ref + aa_position + fc_aa_mut
        if tp == 'mutation':
            m_or_d_dict = mutations
        else:
            m_or_d_dict = deletions
        if info not in sample_mutations:
            sample_mutations[info] = {
                'site': [int(m_or_d.split(',')[0])],
                'frequency': m_or_d_dict[m_or_d][1]
            }
        else:
            sample_mutations[info]['site'].append(int(m_or_d.split(',')[0]))
            sample_mutations[info]['frequency'] += m_or_d_dict[m_or_d][1]


if __name__ == '__main__':
    genome = pd.read_csv(genome_path, index_col=0)

    lineage_10 = {}
    with open(lineage_10_path) as f:
        for line in f.readlines():
            lineage_10[line.split(',')[0]] = {'feature': line.strip().split(',')[2:]}

    lineage_75 = {}
    with open(lineage_75_path) as f:
        for line in f.readlines():
            lineage_75[line.split(',')[0]] = line.strip().split(',')[2:]

    for file in os.listdir(ngs_dir):
        if file.endswith('.csv'):
            file_path = ngs_dir + file
            df = pd.read_csv(file_path)
            filter_count = df['Count'] > 20
            filter_frequency = df['Frequency'] > 10
            filter_quality = df['Average quality'] > 20
            is_snv = df['Type'] == 'SNV'
            is_mnv = df['Type'] == 'MNV'
            is_deletion = df['Type'] == 'Deletion'

            sample_mutations = {}

            # --- mutation calling ---
            df_mutation = df[filter_count
                             & filter_quality
                             & filter_frequency
                             & (is_snv | is_mnv)][
                ['Reference Position', 'Type', 'Length', 'Reference', 'Allele', 'Frequency']]

            mutations = read_dataframe(df_mutation, tp='mutation')

            for m in mutations:
                sites = [int(s) for s in m.split(',')]
                peptide_positions_dict = get_peptide_positions(tp='mutation')
                for p_p in peptide_positions_dict:
                    codon_sites = genome[genome['peptidePos'] == p_p]['genomePos'].values
                    mutation_dict = {}
                    for s_b in peptide_positions_dict[p_p]:
                        mutation_dict[int(s_b.split('/')[0])] = s_b.split('/')[1]
                    for c_s in codon_sites:
                        if c_s not in mutation_dict:
                            mutation_dict[c_s] = genome[genome['genomePos'] == c_s]['nucleotide'].values[0]

                    mutation_dict = dict(sorted(mutation_dict.items(), key=lambda x: x[0]))
                    aa_mut = str(Seq(''.join(mutation_dict.values())).translate())
                    generate_mutation_info(m, aa_mut, tp='mutation')

            # --- deletion calling ---
            df_deletion = df[filter_count
                             & filter_quality
                             & filter_frequency
                             & is_deletion][
                ['Reference Position', 'Type', 'Length', 'Reference', 'Allele', 'Frequency']]

            deletions = read_dataframe(df_deletion, tp='deletion')

            ambiguous_sites = genome[genome['genomePos'].duplicated()]['genomePos'].values
            for d in deletions:
                sites = [int(s) for s in d.split(',')]
                no_ambiguous_site = True
                for a_s in ambiguous_sites:
                    if a_s in sites:
                        no_ambiguous_site = False
                if not no_ambiguous_site:
                    continue
                peptide_positions_dict = get_peptide_positions(tp='deletion')
                if len(peptide_positions_dict) > 0:
                    if len(deletions[d][0]) % 3 != 0:
                        continue
                    p_ps = list(peptide_positions_dict.keys())
                    p_p_f = p_ps[0]
                    p_p_l = p_ps[-1]
                    if len(peptide_positions_dict[p_p_f]) == 3:
                        for p_p in peptide_positions_dict:
                            generate_mutation_info(d, 'del', tp='deletion')
                    else:
                        for p_p in peptide_positions_dict:
                            if p_p != p_p_f:
                                generate_mutation_info(d, 'del', tp='deletion')
                            else:
                                sites_new = []
                                if len(peptide_positions_dict[p_p_f]) == 2:
                                    sites_new.append(int(peptide_positions_dict[p_p_f][0].split('/')[0]) - 1)
                                    sites_new.append(int(peptide_positions_dict[p_p_l][0].split('/')[0]) + 1)
                                    sites_new.append(int(peptide_positions_dict[p_p_l][0].split('/')[0]) + 2)
                                else:
                                    sites_new.append(int(peptide_positions_dict[p_p_f][0].split('/')[0]) - 2)
                                    sites_new.append(int(peptide_positions_dict[p_p_f][0].split('/')[0]) - 1)
                                    sites_new.append(int(peptide_positions_dict[p_p_l][0].split('/')[0]) + 2)
                                bases_new = []
                                for n in range(0, 3):
                                    bases_new.append(genome[genome['genomePos'] == sites_new[n]]
                                                     ['nucleotide'].values[0])

                                aa_mut = str(Seq(''.join(bases_new)).translate())
                                generate_mutation_info(d, aa_mut, tp='deletion')

            # --- get candidate lineages ---
            lineage_info = lineage_10

            for l in lineage_info:
                m_f = {}
                for m in list(sample_mutations.keys()):
                    if m in lineage_info[l]['feature']:
                        m_f[m] = sample_mutations[m]['frequency']
                lineage_info[l]['sample_mutation'] = m_f

            candidate_lineages = []
            for l in lineage_info:
                if len(lineage_info[l]['sample_mutation']) / len(lineage_info[l]['feature']) > 0.3 \
                        and len(lineage_info[l]['sample_mutation']) > 6:
                    candidate_lineages.append(l)

            result_dict = {}
            for l in candidate_lineages:
                for m in lineage_info[l]['sample_mutation']:
                    result_dict.setdefault(m, []).append(l)

            with open(candidate_dir + file, 'w') as f:
                f.write('mutation'
                        + ',' + 'position'
                        + ',' + 'frequency'
                        + ',' + 'lineage'
                        + ',' + 'proportion'
                        + ',' + 'feature_threshold' + '\n')
                for m in sample_mutations:
                    pos = str(genome[(genome['product'] == m.split('_')[0])
                                     & (genome['aaPos'] == int("".join(list(filter(
                                      str.isdigit, m.split('_')[1])))))]['peptidePos'].values[0])
                    fre = str(sample_mutations[m]['frequency'])
                    if m in result_dict:
                        for l in result_dict[m]:
                            pro = str(len(lineage_info[l]['sample_mutation'])) + '/' + str(len(lineage_info[l]['feature']))
                            if m in lineage_75[l]:
                                f_t = 'FV-75'
                            else:
                                f_t = 'FV-10'
                            f.write(m + ',' + pos + ',' + fre + ',' + l + ',' + pro + ',' + f_t + '\n')
                    else:
                        f.write(m + ',' + pos + ',' + fre + ',' + 'N.D.' + ',' + 'N.D.' + ',' + 'N.D.' + '\n')
