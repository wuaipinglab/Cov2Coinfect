import os
import numpy as np
import pandas as pd
from scipy import stats

DIRPATH = '/SSD/yexiao/co_infection/'
lineage_10_path = DIRPATH + 'data/lineage_10.txt'
lineage_75_path = DIRPATH + 'data/lineage_75.txt'
candidate_dir = DIRPATH + 'project/candidate/'

defined_dir = DIRPATH + 'project/defined/'
defined_0_detail_dir = defined_dir + 'defined_0/detail/'
defined_0_summary_dir = defined_dir + 'defined_0/summary/'
defined_1_detail_dir = defined_dir + 'defined_1/detail/'
defined_1_summary_dir = defined_dir + 'defined_1/summary/'
defined_cd_detail_dir = defined_dir + 'defined_cd/detail/'
defined_cd_summary_dir = defined_dir + 'defined_cd/summary/'
defined_cp_detail_dir = defined_dir + 'defined_cp/detail/'
defined_cp_summary_dir = defined_dir + 'defined_cp/summary/'
defined_u_detail_dir = defined_dir + 'defined_u/detail/'
defined_u_summary_dir = defined_dir + 'defined_u/summary/'

dir_list = [
    defined_0_detail_dir,
    defined_0_summary_dir,
    defined_1_detail_dir,
    defined_1_summary_dir,
    defined_cd_detail_dir,
    defined_cd_summary_dir,
    defined_cp_detail_dir,
    defined_cp_summary_dir,
    defined_u_detail_dir,
    defined_u_summary_dir
]

for d in dir_list:
    if not os.path.exists(d):
        os.makedirs(d)

ALL_MUTATIONS_NUM = 92571


def read_file(file_path_fc):
    with open(file_path_fc) as f_fc:
        mutations_fc = {}
        lineages_fc = {}
        for line in f_fc.readlines():
            if not line.startswith('mutation'):
                l = line.split(',')[3]
                if l != 'N.D.':
                    m = line.split(',')[0]
                    f = float(line.split(',')[2])
                    if m not in mutations_fc:
                        mutations_fc[m] = {'lineage': [l], 'frequency': f}
                    else:
                        mutations_fc[m]['lineage'].append(l)
                    if l not in lineages_fc:
                        lineages_fc[l] = {'mutation': [m]}
                    else:
                        lineages_fc[l]['mutation'].append(m)

    return mutations_fc, lineages_fc


def identify_infection(f_l, unique):
    freq = []
    if unique:
        for m in ls_mt2_unique_mutation[f_l]:
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


def output_result(detail_path, summary_path):
    df_final.to_csv(detail_path + file)
    with open(summary_path + file, 'w') as f:
        f.write('lineage' + ',' + 'frequency' + '\n')
        for i in infections:
            f.write(i + ',' + str(infections[i]) + '\n')


if __name__ == '__main__':
    lineage_10 = {}
    with open(lineage_10_path) as f:
        for line in f.readlines():
            lineage_10[line.split(',')[0]] = {'feature': line.strip().split(',')[2:], 'count': int(line.split(',')[1])}

    lineage_75 = {}
    with open(lineage_75_path) as f:
        for line in f.readlines():
            lineage_75[line.split(',')[0]] = line.strip().split(',')[2:]

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

                ls_mt2_unique_mutation = {}
                for l in ls_with_unique_mutation:
                    if len(ls_with_unique_mutation[l]) >= 3:
                        ls_mt2_unique_mutation[l] = ls_with_unique_mutation[l]

                if len(ls_mt2_unique_mutation) > 0:
                    ls_mt2_unique_mutation = dict(sorted(ls_mt2_unique_mutation.items(), key=lambda x: x[0]))
                    ls_mt2_unique_mutation = dict(sorted(ls_mt2_unique_mutation.items(), key=lambda x: lineages[x[0]]['p-value']))
                    first_l = list(ls_mt2_unique_mutation.keys())[0]
                    identify_infection(first_l, unique=True)

                else:
                    ls_shared_mutation = dict(sorted(lineages.items(), key=lambda x: x[0]))
                    ls_shared_mutation = dict(sorted(ls_shared_mutation.items(), key=lambda x: x[1]['p-value']))
                    first_l = list(ls_shared_mutation.keys())[0]
                    same_mutation_lineages = {}
                    for l in ls_shared_mutation:
                        if ls_shared_mutation[l]['mutation'] == ls_shared_mutation[first_l]['mutation']:
                            same_mutation_lineages[l] = lineage_10[l]['count']
                    same_mutation_lineages = dict(sorted(same_mutation_lineages.items(), key=lambda x: x[1], reverse=True))
                    first_l = list(same_mutation_lineages.keys())[0]
                    identify_infection(first_l, unique=False)

            # --- filter unqualified lineages ---
            df_mutation = pd.read_csv(candidate_dir + file)
            df_mutation.loc[~df_mutation['lineage'].isin(infections), ['lineage', 'proportion', 'feature_threshold']] = 'N.D.'
            df_mutation.drop_duplicates(ignore_index=True, inplace=True)
            df_infection = df_mutation[df_mutation['lineage'].isin(infections)]
            df_nd = df_mutation[~df_mutation['mutation'].isin(df_infection['mutation'])]
            df_final = pd.concat([df_infection, df_nd], axis=0, ignore_index=True)

            shared_lineages = []
            for i in list(infections.keys())[::-1]:
                df_unique = df_final.drop_duplicates(subset=['mutation'], keep=False)
                if len(df_unique[df_unique['lineage'] == i]) < 3:
                    df_final.drop(df_final[(df_final['lineage'] == i) & (~df_final.index.isin(df_unique.index))].index, inplace=True)
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

            if len(set(df_final['mutation'])) != 0:
                if len(infections) == 1:
                    output_result(defined_1_detail_dir, defined_1_summary_dir)
                elif len(infections) > 0:
                    nd_proportion = len(set(df_final[df_final['lineage'] == 'N.D.']['mutation'])) / len(set(df_final['mutation']))
                    if 80 <= total_frequency <= 120 and nd_proportion <= 0.3:
                        output_result(defined_cd_detail_dir, defined_cd_summary_dir)
                    else:
                        output_result(defined_cp_detail_dir, defined_cp_summary_dir)
                else:
                    lineages_u = {}
                    for l in lineage_10:
                        lineages_u[l] = {'observed': []}
                        for m in lineage_10[l]['feature']:
                            if m in df_final['mutation'].values:
                                lineages_u[l]['observed'].append(m)

                        p = stats.hypergeom.logsf(len(lineages_u[l]['observed']) - 1, ALL_MUTATIONS_NUM,
                                                  len(lineage_10[l]['feature']), len(set(df_final['mutation'])))
                        lineages_u[l]['p-value'] = p

                    lineages_u = dict(sorted(lineages_u.items(), key=lambda x: x[1]['p-value']))
                    lineage_u = list(lineages_u.keys())[0]

                    observed_proportion = str(len(lineages_u[lineage_u]['observed'])) + '/' + str(len(lineage_10[lineage_u]['feature']))

                    df_final.loc[df_final['mutation'].isin(lineage_10[lineage_u]['feature']),
                                 ['lineage', 'proportion', 'feature_threshold', 'group']] = lineage_u, observed_proportion, 'FV-10', '2'
                    df_final.loc[df_final['mutation'].isin(lineage_75[lineage_u]),
                                 ['lineage', 'proportion', 'feature_threshold', 'group']] = lineage_u, observed_proportion, 'FV-75', '1'

                    mean = np.mean(df_final[df_final['mutation'].isin(lineage_10[lineage_u]['feature'])]['frequency'])
                    infections[lineage_u] = mean
                    output_result(defined_u_detail_dir, defined_u_summary_dir)

            else:
                output_result(defined_0_detail_dir, defined_0_summary_dir)
