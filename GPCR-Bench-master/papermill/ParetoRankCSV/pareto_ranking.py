# %%
import os
import pandas as pd

# %%
os.getcwd()
# %%
test = pd.read_csv("ParetoRanking_metrics_SMO.csv")
test
# %%
files = os.scandir()

pareto_files = []

for file in files:
    if file.name.endswith(".csv"):
        pareto_files.append(file.name)

print(pareto_files)
# %%
pareto_df = pd.concat([pd.read_csv(file) for file in pareto_files], ignore_index=True)

pareto_df
# %%
sort_pareto = pareto_df.copy()
# %%
sort_pareto.sort_values(
    by="Delta Linear Log10 AUC (x10)", inplace=True, ascending=False
)

sort_pareto
# %%
combined_data = pd.read_csv("../combined_data.csv")
combined_data
# %%
combined_data_pareto = pd.concat([combined_data, pareto_df], axis=0)
combined_data_pareto
# %%
combined_data_pareto.to_csv("combined_data_pareto.csv", index=False)

# %%
