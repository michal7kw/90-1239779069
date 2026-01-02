#!/usr/bin/env python3
import os
import pandas as pd
from pathlib import Path

output_dir = os.environ.get('OUTPUT_DIR', '.')
comparisons = ['Neg_vs_Parental', 'Pos_vs_Parental', 'KO_vs_Parental',
               'Pos_vs_Neg', 'KO_vs_Neg', 'KO_vs_Pos']
event_types = ['SE', 'A5SS', 'A3SS', 'MXE', 'RI']

summary = []
for comp in comparisons:
    for event in event_types:
        jc_file = Path(output_dir) / comp / f'{event}.MATS.JC.txt'
        if jc_file.exists():
            df = pd.read_csv(jc_file, sep='\t')
            sig = df[(df['FDR'] < 0.05) & (abs(df['IncLevelDifference']) > 0.1)]
            summary.append({
                'Comparison': comp,
                'Event_Type': event,
                'Total_Events': len(df),
                'Significant_FDR005_dPSI01': len(sig)
            })

summary_df = pd.DataFrame(summary)
summary_df.to_csv(Path(output_dir) / 'splicing_summary.csv', index=False)
print(summary_df.to_string(index=False))
