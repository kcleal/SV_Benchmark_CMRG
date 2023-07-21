import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns

data_ont = []
for caller in ('sniffles', 'cuteSV', 'delly', 'dysgu'):
	with open(f'truvari_{caller}/summary.json', 'r') as f:
		dta = json.load(f)
	d = {'caller': caller, 'platform': 'ONT', 'TP': None}
	d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
	d['TP'] = d['TP-base']
	del d['TP-base']
	data_ont.append(d)
print(pd.DataFrame().from_records(data_ont).round(4).to_markdown())
print()

data_pacbio = []
for caller in ('sniffles', 'cuteSV', 'delly', 'dysgu'):
	with open(f'truvari_{caller}_pacbio/summary.json', 'r') as f:
		dta = json.load(f)
	d = {'caller': caller, 'platform': 'PacBio', 'TP': None}
	d.update({k: v for k, v in dta.items() if k in ('precision', 'recall', 'f1', 'gt_concordance', 'TP-base', 'FP', 'FN')})
	d['TP'] = d['TP-base']
	del d['TP-base']
	data_pacbio.append(d)
print(pd.DataFrame().from_records(data_pacbio).round(4).to_markdown())

data = pd.DataFrame.from_records(data_ont + data_pacbio)
sns.scatterplot(data, x="recall", y="precision", hue="caller", style="platform", s=50)
plt.legend()
plt.savefig("benchmark_result.png")
# plt.show()
