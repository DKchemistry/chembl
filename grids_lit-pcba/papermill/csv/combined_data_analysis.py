# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotx

# %%
sns.set_style("darkgrid")
# plt.style.context(matplotx.styles.duftify(matplotx.styles.dracula))

# %%
#data = pd.read_csv("./ParetoRankCSV/combined_data_pareto.csv")
#data
data = pd.read_csv("combined_data.csv")
data.head()
# %%
# Grouping the data by 'Strain Energy Cutoff' and calculating the mean of the delta values
grouped_data = (
    data.groupby("Strain Energy Cutoff")[
        ["deltaEF1%", "deltaEF5%", "Delta Linear Log10 AUC (x10)", "deltaAUC"]
    ]
    .mean()
    .reset_index()
)

# Display the grouped data
display(grouped_data)

# %%

# %%
# for index, item_indexed in enumerate tracks the index and the item at that index
# here we are iterating over the columns of the grouped data
plt.figure(figsize=(20, 10))
for i, column in enumerate(
    ["deltaEF1%", "deltaEF5%", "Delta Linear Log10 AUC (x10)", "deltaAUC"], 1
):
    plt.subplot(2, 2, i)
    sns.barplot(
        x="Strain Energy Cutoff",
        y=column,
        data=grouped_data,
    )
    plt.xlabel("Strain Energy Cutoff")
    plt.ylabel(column)
    plt.xticks(rotation=45)
    plt.title(f"Mean {column} vs Strain Energy Cutoff")

plt.tight_layout()
plt.show()

# %%
sns.barplot(
    x="Strain Energy Cutoff", y="Delta Linear Log10 AUC (x10)", data=grouped_data
)
plt.title("Mean Delta Linear Log10 AUC (x10) vs Strain Energy Cutoff")
plt.xlabel("Strain Energy Cutoff")
plt.ylabel("Delta Linear Log10 AUC (x10)")
plt.xticks(rotation=45)
plt.show()

# %%
plt.figure(figsize=(20, 15))

# Define the metrics to plot
metrics = ["deltaEF1%", "deltaEF5%", "Delta Linear Log10 AUC (x10)", "deltaAUC"]

# Creating box plots for each metric
for i, metric in enumerate(metrics, 1):
    plt.subplot(2, 2, i)
    sns.boxplot(x="Strain Energy Cutoff", y=metric, data=data)
    plt.xticks(rotation=45)
    plt.title(f"Box Plot of {metric} by Strain Energy Cutoff")

# Adjust layout
plt.tight_layout()

# Show the plots without legends
plt.show()

# %%
# Filter the dataset for strain energy cutoff of 4
cutoff_4_data = data[data["Strain Energy Cutoff"] == "4"]

# Create a bar plot of Delta Linear Log10 AUC (x10) for each protein at cutoff 4
plt.figure(figsize=(10, 8))
sns.barplot(x="Delta Linear Log10 AUC (x10)", y="Protein", data=cutoff_4_data, ci=None)
plt.title("Delta Linear Log10 AUC (x10) for Each Protein at Strain Energy Cutoff 4")
plt.xlabel("Delta Linear Log10 AUC (x10)")
plt.ylabel("Protein")

plt.show()

# %%
data.columns

# %%
# Create a barplot plot of deltaEF1% for cutoff 4
plt.figure(figsize=(10, 8))
sns.barplot(x="deltaEF1%", y="Protein", data=cutoff_4_data, ci=None)
plt.title("deltaEF1% for Each Protein at Strain Energy Cutoff 4")
plt.xlabel("deltaEF1%")
plt.ylabel("Protein")

plt.show()


# %%
# Create a barplot plot of deltaEF5% for cutoff 4
plt.figure(figsize=(10, 8))
sns.barplot(x="deltaEF5%", y="Protein", data=cutoff_4_data, ci=None)
plt.title("deltaEF5% for Each Protein at Strain Energy Cutoff 4")
plt.xlabel("deltaEF5%")
plt.ylabel("Protein")

plt.show()

# %%
# filter the data for strain energy cut off of 5.0
cutoff_5_data = data[data["Strain Energy Cutoff"] == "5.0"]

# Create a bar plot of Delta Linear Log10 AUC (x10) for each protein at cutoff 5
plt.figure(figsize=(10, 8))
sns.barplot(x="Delta Linear Log10 AUC (x10)", y="Protein", data=cutoff_5_data, ci=None)
plt.title("Delta Linear Log10 AUC (x10) for Each Protein at Strain Energy Cutoff 5")
plt.xlabel("Delta Linear Log10 AUC (x10)")
plt.ylabel("Protein")

plt.show()


# %%
# Create a barplot plot of deltaEF1% for cutoff 5
plt.figure(figsize=(10, 8))
sns.barplot(x="deltaEF1%", y="Protein", data=cutoff_5_data, ci=None)
plt.title("deltaEF1% for Each Protein at Strain Energy Cutoff 5")
plt.xlabel("deltaEF1%")
plt.ylabel("Protein")

plt.show()

# %%
# Create a barplot plot of deltaEF5% for cutoff 5
plt.figure(figsize=(10, 8))
sns.barplot(x="deltaEF5%", y="Protein", data=cutoff_5_data, ci=None)
plt.title("deltaEF5% for Each Protein at Strain Energy Cutoff 5")
plt.xlabel("deltaEF5%")
plt.ylabel("Protein")

plt.show()

# %%
# filter the data for strain energy cut off of 7.0
cutoff_7_data = data[data["Strain Energy Cutoff"] == "7.0"]

# Create a bar plot of Delta Linear Log10 AUC (x10) for each protein at cutoff 5
plt.figure(figsize=(10, 8))
sns.barplot(x="Delta Linear Log10 AUC (x10)", y="Protein", data=cutoff_7_data, ci=None)
plt.title("Delta Linear Log10 AUC (x10) for Each Protein at Strain Energy Cutoff 7")
plt.xlabel("Delta Linear Log10 AUC (x10)")
plt.ylabel("Protein")

plt.show()

# %%
