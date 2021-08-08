import pandas as pd

DIRPATH = '/SSD/yexiao/co_infection/'
variant_surveillance_path = DIRPATH + 'data/variant_surveillance.tsv'
lineage_10_path = DIRPATH + 'data/lineage_10.txt'
lineage_75_path = DIRPATH + 'data/lineage_75.txt'
mutation_num_path = DIRPATH + 'data/mutation_num.txt'

df = pd.read_csv(variant_surveillance_path, delimiter='\t').dropna(axis=0, subset=['AA Substitutions', 'Pango lineage'])

mutations = []
lineages = {}
for i in df.index:
    l = df.loc[i, 'Pango lineage']
    if l not in lineages:
        lineages[l] = {'count': 1, 'mutations': {}}
    else:
        lineages[l]['count'] += 1

    substitutions = df.loc[i, 'AA Substitutions'][1:-1].split(',')
    for s in substitutions:
        if 'ins' not in s and s != '':
            mutations.append(s)
            lineages[l]['mutations'][s] = lineages[l]['mutations'].get(s, 0) + 1

mutations = list(set(mutations))

with open(mutation_num_path, 'w') as f:
    f.write('Num'+'\n')
    f.write(str(len(mutations))+'\n')

lineages = dict(sorted(lineages.items(), key=lambda x: x[0]))
defining_SNPs_75 = {}
defining_SNPs_10 = {}
for l in lineages:
    if l != 'None' and lineages[l]['count'] / len(df) > 0.0001:
        for m in lineages[l]['mutations']:
            if lineages[l]['mutations'][m] / lineages[l]['count'] > 0.1:
                defining_SNPs_10.setdefault(l, []).append(m)
                if l not in defining_SNPs_75:
                    defining_SNPs_75[l] = []
            if lineages[l]['mutations'][m] / lineages[l]['count'] > 0.75:
                defining_SNPs_75[l].append(m)

with open(lineage_10_path, 'w') as f:
    for l in defining_SNPs_10:
        f.write(l + ',' + str(lineages[l]['count']) + ',' + ','.join(defining_SNPs_10[l]) + '\n')

with open(lineage_75_path, 'w') as f:
    for l in defining_SNPs_75:
        f.write(l + ',' + str(lineages[l]['count']) + ',' + ','.join(defining_SNPs_75[l]) + '\n')
