import pandas as pd
import matplotlib.pyplot as plt

# Read the data into a DataFrame
df = pd.read_csv('../data/summary_statistics_new2.csv')

# Extract labels and values for the bar plot
labels = df.columns.tolist()
values = df.iloc[0].tolist()

# Create a simple bar plot
fig, ax = plt.subplots(figsize=(10, 6))
ax.bar(labels, values)

# Add labels and title
ax.set_ylabel('Count')
ax.set_title('Sample Data Overview')
plt.xticks(rotation=45, ha='right')

# Adjust layout for readability
plt.tight_layout()
plt.savefig('../data/general_stats_count.png', dpi=800)
