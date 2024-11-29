import pandas as pd
import json
import matplotlib.pyplot as plt
import seaborn as sns
import glob


def process_benchmark_data(platform, callers, benchmarks):
    data = []
    p = 'results/' if len(glob.glob('results')) > 0 else '' 
    for caller in callers:
        for bench in benchmarks:
            if bench == 'CMRG':
                path = glob.glob(f'{p}truvari_{bench}_{platform}_{caller}*/summary.json')
            else:
                path = glob.glob(f'{p}truvari_{bench}_{platform}_{caller}*/refine.variant_summary.json')
            if not path:
                continue
            with open(path[0], 'r') as f:
                dta = json.load(f)
            d = {'caller': caller, 'platform': platform, 'Benchmark': bench, 'TP': None}
            metrics = ['precision', 'recall', 'f1', 'TP-base', 'FP', 'FN']
            if bench == 'CMRG':
                metrics.append('gt_concordance')
            d.update({k: v for k, v in dta.items() if k in metrics})
            d['TP'] = d['TP-base']
            del d['TP-base']
            data.append(d)
    if not data:
        return []
    return pd.DataFrame.from_records(data, index=None).sort_values('Benchmark')


def plot_benchmark_results(df, output_file="benchmark_result.png"):
    for bench in ('CMRG', 'GIAB'):
        plot = sns.relplot(
            data=df[df['Benchmark'] == bench], x="recall", y="precision",
            col="platform", hue="caller",
            kind="scatter"
        )
        plot.fig.suptitle(bench)
        plt.savefig(f'benchmark_result.{bench}.png')
        plt.close()


if __name__ == "__main__":

    # Process ONT data
    ont_callers = ['sniffles', 'cuteSV', 'severus', 'delly', 'ngsep', 'dysgu']
    ont_data = process_benchmark_data('ont', ont_callers, ['CMRG', 'GIAB'])
    if len(ont_data):
        with open('ont_results.md', 'w') as ont_file:
            ont_file.write("ONT Results:\n")
            ont_file.write(ont_data.round(4).to_markdown(index=False))

    # Process PacBio data
    pacbio_callers = ['sniffles', 'cuteSV', 'severus', 'delly', 'ngsep', 'sawfish', 'dysgu']
    pacbio_data = process_benchmark_data('pacbio', pacbio_callers, ['CMRG', 'GIAB'])

    if len(pacbio_data):
        with open('pacbio_results.md', 'w') as pacbio_file:
            pacbio_file.write("PacBio Results:\n")
            pacbio_file.write(pacbio_data.round(4).to_markdown(index=False))
     
    if len(ont_data) and len(pacbio_data):
        final_data = pd.concat([ont_data, pacbio_data])
        plot_benchmark_results(final_data)

