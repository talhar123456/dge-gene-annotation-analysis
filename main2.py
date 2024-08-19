import pandas as pd
from collections import Counter
from scipy.stats import hypergeom

# Function to parse the GO annotation file
def parse_go_annotations(file_path):
    go_annotations = []
    with open('human_GO.gaf', 'r') as f:
        for line in f:
            if line.startswith('!'):
                continue
            parts = line.strip().split('\t')
            if parts[0] == 'UniProtKB' and parts[8] == 'P':
                go_annotations.append({
                    'Accession': parts[1],
                    'Gene': parts[2],
                    'GO_ID': parts[4]
                })
    go_df = pd.DataFrame(go_annotations)
    print("First few rows of the DataFrame:")
    print(go_df.head())  # Debugging: Print the first few rows
    print("Columns of the DataFrame:")
    print(go_df.columns)  # Debugging: Print the column names
    return go_df

# Function to get common and rare GO terms
def common_rare_go_terms(go_annotations, de_genes, n=20):
    filtered_go = go_annotations[go_annotations['Gene'].isin(de_genes)]
    print(f"Filtered GO annotations for DE genes:\n{filtered_go.head()}")  # Debugging: Print filtered annotations
    go_counts = Counter(filtered_go['GO_ID'])
    most_common = go_counts.most_common(n)
    least_common = sorted(go_counts.items(), key=lambda x: (x[1], x[0]))[:n]
    return most_common, least_common

# Function to compute hypergeometric p-value
def hypergeometric_test(m, N, n, k):
    p_value = hypergeom.sf(k-1, N, m, n)
    return p_value

# Adjusted p-value function (from Exercise 1)
def adjust_pval(pvals):
    n = len(pvals)
    adj_pvals = [min(1, p * n) for p in pvals]
    return adj_pvals

# Function to perform enrichment analysis
def enrichment_analysis(go_annotations, de_genes, background_genes):
    go_counts = Counter(go_annotations['GO_ID'])
    results = []
    for go_id, count in go_counts.items():
        m = count
        N = len(background_genes)
        n = len(de_genes)
        k = len(go_annotations[(go_annotations['GO_ID'] == go_id) & (go_annotations['Gene'].isin(de_genes))])
        if k > 0:
            p_value = hypergeometric_test(m, N, n, k)
            results.append({'GO_ID': go_id, 'p_value': p_value})
    results_df = pd.DataFrame(results)
    print(f"Results DataFrame before adjustment:\n{results_df.head()}")  # Debugging: Print results before p-value adjustment
    if not results_df.empty:
        results_df['adj_p_value'] = adjust_pval(results_df['p_value'])
    return results_df.sort_values(by='adj_p_value').head(20)

# Example usage
go_annotations = parse_go_annotations('human_GO.gaf')
background_genes = go_annotations['Gene'].unique()
de_genes = ['DNAJC25-GNG10', 'IGKV3-7']  # Use actual gene names from the parsed GO annotations

# Get most and least common GO terms
most_common, least_common = common_rare_go_terms(go_annotations, de_genes, n=20)
print("Most common GO terms:", most_common)
print("Least common GO terms:", least_common)

# Perform enrichment analysis
results = enrichment_analysis(go_annotations, de_genes, background_genes)
print("Enrichment analysis results:")
print(results)
