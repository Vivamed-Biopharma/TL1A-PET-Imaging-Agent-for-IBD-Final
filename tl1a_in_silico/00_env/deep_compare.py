#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

BASE = os.path.join(os.path.dirname(os.path.dirname(__file__)))

# Inputs
dev = pd.read_csv(os.path.join(BASE, 'developability.csv'))
dar = pd.read_csv(os.path.join(BASE, 'dar.csv'))
para = pd.read_csv(os.path.join(BASE, 'paratope.csv'))
immu = pd.read_csv(os.path.join(BASE, 'immunogenicity.csv'))
manu = pd.read_csv(os.path.join(BASE, 'manufacturability.csv'))
rank = pd.read_csv(os.path.join(BASE, 'composite_ranking.csv'))

# Merge master
master = dev.merge(dar[['Clone','P_DAR_1_2','P_DAR_ge4','E_DAR']], on='Clone')\
             .merge(para, on='Clone')\
             .merge(immu, on='Clone')\
             .merge(manu[['Clone','AggProxyMax_VH','AggProxyMax_VL']], on='Clone')\
             .merge(rank, on='Clone')

# Lead sets
originals = [f'Fab{i:02d}' for i in range(1,13)]
new_leads = ['Fab79','Fab96','Fab122']
orig_leads = ['Fab04','Fab03','Fab09']
compare = new_leads + orig_leads

subset = master[master['Clone'].isin(compare)].copy()

# Normalize pillars to z-scores (better= higher except for risks where lower is better)
cols_higher = ['P_DAR_1_2','Paratope','Score']
cols_lower  = ['P_DAR_ge4','E_DAR','AggProxyMax_VH','AggProxyMax_VL','ImmBurden_VH','ImmBurden_VL']

for col in cols_higher:
    z = (master[col] - master[col].mean())/master[col].std(ddof=0)
    subset[col+'_z'] = z.reindex(subset.index)
for col in cols_lower:
    z = - (master[col] - master[col].mean())/master[col].std(ddof=0)
    subset[col+'_z'] = z.reindex(subset.index)

subset['Composite_Z'] = subset[[c+'_z' for c in cols_higher+cols_lower]].mean(axis=1)

# Export
out_csv = os.path.join(BASE, 'lead_deep_compare.csv')
subset_out = subset[['Clone','Score','P_DAR_1_2','P_DAR_ge4','E_DAR','Paratope','DR3_adj','AggProxyMax_VH','AggProxyMax_VL','ImmBurden_VH','ImmBurden_VL','Composite_Z']]
subset_out.to_csv(out_csv, index=False)

# Figure: radar-like bar plot per clone
melt_cols = ['P_DAR_1_2_z','P_DAR_ge4_z','E_DAR_z','Paratope_z','AggProxyMax_VH_z','AggProxyMax_VL_z','ImmBurden_VH_z','ImmBurden_VL_z','Score_z']
plot_df = subset[['Clone'] + melt_cols].melt('Clone', var_name='Metric', value_name='Z')
plt.figure(figsize=(10,5))
sns.barplot(data=plot_df, x='Metric', y='Z', hue='Clone')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.savefig(os.path.join(BASE, 'fig_lead_compare.png'), dpi=200)
print('Wrote', out_csv)
