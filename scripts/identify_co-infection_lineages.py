import os
import numpy as np
import pandas as pd
from scipy import stats

DIRPATH = '/SSD/yexiao/co_infection/'
lineage_10_path = DIRPATH + 'data/lineage_10.txt'
candidate_dir = DIRPATH + 'examples/candidate/'

defined_1_detail_dir = DIRPATH + 'examples/defined/defined_1/detail/'
defined_1_summary_dir = DIRPATH + 'examples/defined/defined_1/summary/'
defined_2_detail_dir = DIRPATH + 'examples/defined/defined_2/detail/'
defined_2_summary_dir = DIRPATH + 'examples/defined/defined_2/summary/'
defined_3_detail_dir = DIRPATH + 'examples/defined/defined_3/detail/'
defined_3_summary_dir = DIRPATH + 'examples/defined/defined_3/summary/'
defined_u_detail_dir = DIRPATH + 'examples/defined/defined_u/detail/'
defined_u_summary_dir = DIRPATH + 'examples/defined/defined_u/summary/'

dir_list = [defined_1_detail_dir,
            defined_1_summary_dir,
            defined_2_detail_dir,
            defined_2_summary_dir,
            defined_3_detail_dir,
            defined_3_summary_dir,
            defined_u_detail_dir,
            defined_u_summary_dir]

for d in dir_list:
    if not os.path.exists(d):
        os.makedirs(d)

ALL_MUTATIONS_NUM = 92571


def read_file(fc_file_path):
    with open(fc_file_path) as fc_f:
        fc_mutations = {}
        fc_lineages = {}
        for line in fc_f.readlines():
            if not line.startswith('mutation'):
                l = line.split(',')[3]
                if l != 'N.D.':
                    m = line.split(',')[0]
                    f = float(line.split(',')[2])
                    if m not in fc_mutations:
                        fc_mutations[m] = {'lineage': [l], 'frequency': f}
                    else:
                        fc_mutations[m]['lineage'].append(l)
                    if l not in fc_lineages:
                        fc_lineages[l] = {'mutation': [m]}
                    else:
                        fc_lineages[l]['mutation'].append(m)

    return fc_mutations, fc_lineages


def identify_infection(f_l, unique):
    freq = []
    if unique:
        for m in ls_more_than_1_unique_mutation[f_l]:
            freq.append(mutations[m]['frequency'])
    else:
        for m in lineages[f_l]['mutation']:
            freq.append(mutations[m]['frequency'])

    for m in lineages[f_l]['mutation']:
        mutations[m]['lineage'].remove(f_l)
        if len(mutations[m]['lineage']) == 0:
            del mutations[m]

    freq_std = np.std(freq, ddof=0)
    if freq_std < 20 \
            and len(lineages[f_l]['mutation']) > 6 \
            and len(lineages[f_l]['mutation']) / len(lineage_10[f_l]['feature']) > 0.3:
        mean = np.mean(freq)
        infections[f_l] = mean
        for m in lineages[f_l]['mutation']:
            if m in mutations:
                mutations[m]['frequency'] -= mean
                if mutations[m]['frequency'] <= 10:
                    for l in mutations[m]['lineage']:
                        lineages[l]['mutation'].remove(m)
                        if len(lineages[l]['mutation']) == 0:
                            del lineages[l]
                    del mutations[m]
    del lineages[f_l]


if __name__ == '__main__':
    lineage_10 = {}
    with open(lineage_10_path) as f:
        for line in f.readlines():
            lineage_10[line.split(',')[0]] = {'feature': line.strip().split(',')[2:], 'count': int(line.split(',')[1])}

    for file in os.listdir(candidate_dir):
        if file.endswith('.csv'):

            mutations, lineages = read_file(candidate_dir + file)

            infections = {}

            # --- get enriched lineages ---
            while len(lineages) != 0:
                for l in lineages:
                    p = stats.hypergeom.logsf(len(lineages[l]['mutation']) - 1, ALL_MUTATIONS_NUM,
                                              len(lineage_10[l]['feature']), len(mutations))
                    lineages[l]['p-value'] = p

                ls_with_unique_mutation = {}
                for m in mutations:
                    if not m.endswith('del'):
                        ls = mutations[m]['lineage']
                        if len(ls) == 1:
                            ls_with_unique_mutation.setdefault(ls[0], []).append(m)

                ls_more_than_1_unique_mutation = {}
                for l in ls_with_unique_mutation:
                    if len(ls_with_unique_mutation[l]) >= 3:
                        ls_more_than_1_unique_mutation[l] = ls_with_unique_mutation[l]

                if len(ls_more_than_1_unique_mutation) > 0:
                    ls_more_than_1_unique_mutation = dict(
                        sorted(ls_more_than_1_unique_mutation.items(), key=lambda x: x[0])
                    )
                    ls_more_than_1_unique_mutation = dict(
                        sorted(ls_more_than_1_unique_mutation.items(), key=lambda x: lineages[x[0]]['p-value'])
                    )
                    first_l = list(ls_more_than_1_unique_mutation.keys())[0]
                    identify_infection(first_l, unique=True)

                else:
                    ls_shared_mutation = dict(sorted(lineages.items(), key=lambda x: x[0]))
                    ls_shared_mutation = dict(
                        sorted(ls_shared_mutation.items(), key=lambda x: x[1]['p-value'])
                    )
                    first_l = list(ls_shared_mutation.keys())[0]
                    same_mutation_lineages = {}
                    for l in ls_shared_mutation:
                        if ls_shared_mutation[l]['mutation'] == ls_shared_mutation[first_l]['mutation']:
                            same_mutation_lineages[l] = lineage_10[l]['count']
                    same_mutation_lineages = dict(
                        sorted(same_mutation_lineages.items(), key=lambda x: x[1], reverse=True)
                    )
                    first_l = list(same_mutation_lineages.keys())[0]
                    identify_infection(first_l, unique=False)

            # --- filter unqualified lineages ---
            df_mutation = pd.read_csv(candidate_dir + file)
            df_mutation.loc[~df_mutation['lineage'].isin(infections),
                            ['lineage', 'proportion', 'feature_threshold']] = 'N.D.'
            df_mutation.drop_duplicates(ignore_index=True, inplace=True)
            df_infection = df_mutation[df_mutation['lineage'].isin(infections)]
            df_nd = df_mutation[~df_mutation['mutation'].isin(df_infection['mutation'])]
            df_final = pd.concat([df_infection, df_nd], axis=0, ignore_index=True)

            shared_lineages = []
            for i in list(infections.keys())[::-1]:
                df_unique = df_final.drop_duplicates(subset=['mutation'], keep=False)
                if len(df_unique[df_unique['lineage'] == i]) < 3:
                    df_final.drop(df_final[(df_final['lineage'] == i)
                                           & (~df_final.index.isin(df_unique.index))].index, inplace=True)
                    df_final.loc[df_final['lineage'] == i, ['lineage', 'proportion', 'feature_threshold']] = 'N.D.'
                    shared_lineages.append(i)
            for l in shared_lineages:
                infections.pop(l)

            for i in df_final.index:
                if df_final.loc[i, 'lineage'] != 'N.D.':
                    group = (list(infections.keys()).index(df_final.loc[i, 'lineage']) + 1) * 2
                    if df_final.loc[i, 'feature_threshold'] == 'FV-75':
                        group -= 1
                else:
                    group = 0
                df_final.loc[i, 'group'] = str(group)

            df_final.sort_values('position', ignore_index=True, inplace=True)

            # --- get defined lineages ---
            total_frequency = 0
            for l in infections:
                total_frequency += infections[l]

            nd_proportion = len(set(df_final[df_final['lineage'] == 'N.D.']['mutation']))/len(set(df_final['mutation']))

            if 80 <= total_frequency <= 120 and nd_proportion <= 0.3 and 0 < len(infections) <= 3:

                if len(infections) == 1:
                    df_final.to_csv(defined_1_detail_dir + file)
                    with open(defined_1_summary_dir + file, 'w') as f:
                        f.write('lineage' + ',' + 'frequency' + '\n')
                        for i in infections:
                            f.write(i + ',' + str(infections[i]) + '\n')

                if len(infections) == 2:
                    df_final.to_csv(defined_2_detail_dir + file)
                    with open(defined_2_summary_dir + file, 'w') as f:
                        f.write('lineage' + ',' + 'frequency' + '\n')
                        for i in infections:
                            f.write(i + ',' + str(infections[i]) + '\n')

                if len(infections) == 3:
                    df_final.to_csv(defined_3_detail_dir + file)
                    with open(defined_3_summary_dir + file, 'w') as f:
                        f.write('lineage' + ',' + 'frequency' + '\n')
                        for i in infections:
                            f.write(i + ',' + str(infections[i]) + '\n')

            else:
                df_final.to_csv(defined_u_detail_dir + file)
                with open(defined_u_summary_dir + file, 'w') as f:
                    f.write('lineage' + ',' + 'frequency' + '\n')
                    for i in infections:
                        f.write(i + ',' + str(infections[i]) + '\n')
