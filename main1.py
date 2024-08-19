import pandas as pd
import numpy as np
from scipy.stats import ttest_ind

# Step 1: Read the data
file_path = 'GSE10718.tsv'
data = pd.read_csv(file_path, sep='\t')

# Step 2: Calculate Means and Variances
def calculate_means_variances(data, time_point):
    result = []
    time_columns = {
        '2h': ['Air_2h_A', 'Air_2h_B', 'Air_2h_C', 'Smoke_2h_A', 'Smoke_2h_B', 'Smoke_2h_C'],
        '4h': ['Air_4h_A', 'Air_4h_B', 'Air_4h_C', 'Smoke_4h_A', 'Smoke_4h_B', 'Smoke_4h_C'],
        '24h': ['Air_24h_A', 'Air_24h_B', 'Air_24h_C', 'Smoke_24h_A', 'Smoke_24h_B', 'Smoke_24h_C']
    }
    
    for index, row in data.iterrows():
        air_values = row[time_columns[time_point][:3]].values
        smoke_values = row[time_columns[time_point][3:]].values
        
        mean_air = np.mean(air_values)
        mean_smoke = np.mean(smoke_values)
        var_air = np.var(air_values, ddof=1)  # Sample variance
        var_smoke = np.var(smoke_values, ddof=1)  # Sample variance
        
        result.append((row['ProbeID'], row['Symbol'], mean_air, var_air, mean_smoke, var_smoke))
    
    return pd.DataFrame(result, columns=['ProbeID', 'Symbol', 'Mean_Air', 'Var_Air', 'Mean_Smoke', 'Var_Smoke'])

# Step 3: Log Fold Change Calculation
def get_LFC(mean_air, mean_smoke):
    return np.log2(mean_smoke / mean_air)

def add_logFC(df):
    df['LogFC'] = df.apply(lambda row: get_LFC(row['Mean_Air'], row['Mean_Smoke']), axis=1)
    return df

# Step 4: p-value Calculation
def get_pval(air_values, smoke_values):
    _, p_value = ttest_ind(air_values, smoke_values, equal_var=True)
    return p_value

def add_pvals(data, df, time_point):
    time_columns = {
        '2h': ['Air_2h_A', 'Air_2h_B', 'Air_2h_C', 'Smoke_2h_A', 'Smoke_2h_B', 'Smoke_2h_C'],
        '4h': ['Air_4h_A', 'Air_4h_B', 'Air_4h_C', 'Smoke_4h_A', 'Smoke_4h_B', 'Smoke_4h_C'],
        '24h': ['Air_24h_A', 'Air_24h_B', 'Air_24h_C', 'Smoke_24h_A', 'Smoke_24h_B', 'Smoke_24h_C']
    }
    
    p_values = []
    for index, row in data.iterrows():
        # Convert values to numeric type
        air_values = pd.to_numeric(row[time_columns[time_point][:3]], errors='coerce').dropna().values
        smoke_values = pd.to_numeric(row[time_columns[time_point][3:]], errors='coerce').dropna().values
        
        # Check if there are NaNs after conversion
        if np.any(np.isnan(air_values)) or np.any(np.isnan(smoke_values)):
            continue
        
        p_value = get_pval(air_values, smoke_values)
        p_values.append(p_value)
    
    df['p_value'] = p_values
    return df

# Step 5: Multiple Testing Correction
def adjust_pval(df):
    n = len(df)
    df['p_adjusted'] = df['p_value'].apply(lambda p: min(1.0, n * p))
    return df

# Step 6: Filter and Report
def get_top_probes(df, top_n=10):
    return df.nsmallest(top_n, 'p_adjusted')

# Main script to execute all steps
if __name__ == '__main__':
    # Calculate means and variances for each time point
    means_variances_2h = calculate_means_variances(data, '2h')
    means_variances_4h = calculate_means_variances(data, '4h')
    means_variances_24h = calculate_means_variances(data, '24h')
    
    # Add logFC to the means and variances data
    means_variances_2h = add_logFC(means_variances_2h)
    means_variances_4h = add_logFC(means_variances_4h)
    means_variances_24h = add_logFC(means_variances_24h)
    
    # Add p-values to the means and variances data
    means_variances_2h = add_pvals(data, means_variances_2h, '2h')
    means_variances_4h = add_pvals(data, means_variances_4h, '4h')
    means_variances_24h = add_pvals(data, means_variances_24h, '24h')
    
    # Adjust p-values using Bonferroni correction
    means_variances_2h = adjust_pval(means_variances_2h)
    means_variances_4h = adjust_pval(means_variances_4h)
    means_variances_24h = adjust_pval(means_variances_24h)
    
    # Get top 10 probes at each time point
    top_probes_2h = get_top_probes(means_variances_2h)
    top_probes_4h = get_top_probes(means_variances_4h)
    top_probes_24h = get_top_probes(means_variances_24h)
    
    # Print results
    print("Top 10 probes at 2h:")
    print(top_probes_2h[['ProbeID', 'Symbol', 'LogFC', 'p_value', 'p_adjusted']])
    
    print("\nTop 10 probes at 4h:")
    print(top_probes_4h[['ProbeID', 'Symbol', 'LogFC', 'p_value', 'p_adjusted']])
    
    print("\nTop 10 probes at 24h:")
    print(top_probes_24h[['ProbeID', 'Symbol', 'LogFC', 'p_value', 'p_adjusted']])
