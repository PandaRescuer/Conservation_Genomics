import random
import pandas as pd
import numpy as np
from collections import defaultdict

# ========== Parameter settings ==========
N_REPEATS = 25          # Number of bootstrap resampling replicates
SAMPLE_SIZE = 5          # Number of individuals sampled per population
TOTAL_ALLELES = 2 * SAMPLE_SIZE  # =10
RANDOM_SEED = 42

# ========== File paths ==========
LOF_file_dir = 'LOF_snp_RMA'
DEL_file_dir = 'mis_snp_DELETERIOUS'
TOL_file_dir = 'mis_snp_TOLERATED'
Mis_file_dir = 'missense_snp_RMA'
SYN_file_dir = 'synonymous_snp_RMA'
CNE_file_dir = 'CNE_snp_RMA'

# output_path = 'merge_SFS_bootstrap_CI_20x.txt'

# mean_output_path = output_path.replace(
#     "_CI_20x.txt", "_mean_20x.txt"
# )

output_path = 'merge_SFS_bootstrap_CI_all.txt'

mean_output_path = output_path.replace(
    "_CI_all.txt", "_mean_all.txt"
)

if mean_output_path == output_path:
    mean_output_path = output_path + ".mean_only.txt"

file_dict = {
    'LOF': LOF_file_dir + '_filter_MAF.vcf',
    'Deleterous': DEL_file_dir + '_filter_MAF.vcf',
    'Tolerated': TOL_file_dir + '_filter_MAF.vcf',
    'Missense': Mis_file_dir + '_filter_MAF.vcf',
    'Synonymous': SYN_file_dir + '_filter_MAF.vcf',
    'CNE': CNE_file_dir + '_filter_MAF.vcf'
}
categories = ['CNE', 'Deleterous', 'LOF', 'Missense', 'Synonymous', 'Tolerated']

# ========== Hard-coded population sample IDs (individuals with depth > 20) ==========
SC_WLP_ids = {'GP01', 'GP02', 'GP03', 'GP04', 'GP05', 'GP06', 'GP07', 'GP09', 'GP10'}
QLI_ids = {'QLI05c', 'QLI06c', 'QLI07c', 'QLI04c', 'GP55'}
SC_CAP_ids = {'GP12', 'GP13', 'GP14', 'GP15', 'GP16', 
              'GP17', 'GP18', 'GP19', 'GP23', 'GP24',
              'GP25', 'GP27', 'GP28', 'GP29', 'GP30', 
              'GP31', 'GP32', 'GP33', 'GP34', 'GP35', 
              'GP38', 'GP40', 'GP43', 'GP44', 'GP46'}
MIX_CAP_ids = {'GP39','GP41','GP42', 'GP45','GP47', 'GP48', 'GP49', 'GP50', 'GP51', 'GP52', 'GP53'}

pop_to_ids = {
    'SC_WLP': SC_WLP_ids,
    'QLI': QLI_ids,
    'SC_CAP': SC_CAP_ids,
    'MIX_CAP': MIX_CAP_ids
}
pop_names = list(pop_to_ids.keys())

# ========== Read pop.txt to build the ID -> NAME mapping ==========
id_to_name = {}
with open('D:/重测序项目/测序部分/pop.txt', 'r') as f:
    header = f.readline().strip().split('\t')
    id_idx = 0          # ID column
    name_idx = 2        # NAME column
    for line in f:
        cols = line.strip().split('\t')
        if len(cols) <= max(id_idx, name_idx):
            continue
        sample_id = cols[id_idx]
        name = cols[name_idx]
        if sample_id and name:
            id_to_name[sample_id] = name

print(f"Loaded {len(id_to_name)} ID -> NAME mappings.")

# ========== Convert ID lists to NAME lists ==========
pop_to_names = {}
for pop, ids in pop_to_ids.items():
    names = set()
    for sid in ids:
        if sid in id_to_name:
            names.add(id_to_name[sid])
        else:
            print(f"Warning: ID {sid} has no matching NAME in pop.txt")
    pop_to_names[pop] = names
    print(f"Number of NAME entries for population {pop}: {len(names)}")

# ========== Read the sample list from the VCF header (using NAME) ==========
def get_vcf_samples(vcf_path):
    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                return line.strip().split()[9:]
    return None

first_vcf = list(file_dict.values())[0]
samples_vcf = get_vcf_samples(first_vcf)
if samples_vcf is None:
    raise ValueError("Unable to read samples from the VCF header")

print(f"\nTotal VCF samples: {len(samples_vcf)}, first 10: {samples_vcf[:10]}")

# ========== Build population indices (matching by NAME) ==========
pop_indices = {pop: [] for pop in pop_names}
for idx, sample in enumerate(samples_vcf):
    for pop, names in pop_to_names.items():
        if sample in names:
            pop_indices[pop].append(idx)
            break

print("\nNumber of samples matched per population:")
for pop in pop_names:
    print(f"  {pop}: {len(pop_indices[pop])}")
    if len(pop_indices[pop]) == 0:
        print(f"    Example expected NAMEs: {list(pop_to_names[pop])[:5]}")
        found = [n for n in pop_to_names[pop] if n in set(samples_vcf)]
        print(f"    Actually present in VCF: {found[:5]}")

# Check sample counts
for pop in pop_names:
    if len(pop_indices[pop]) == 0:
        raise ValueError(f"Population '{pop}' has no matching samples; please check the mapping.")
    if len(pop_indices[pop]) < SAMPLE_SIZE:
        print(f"Warning: population '{pop}' has only {len(pop_indices[pop])} samples; sampling with replacement will be used.")

# ========== Read GT data from all VCFs ==========
vcf_data = {}
for cat, path in file_dict.items():
    with open(path, 'r') as f:
        lines = f.readlines()
    start = 0
    for i, line in enumerate(lines):
        if line.startswith('#'):
            continue
        else:
            start = i
            break
    gt_rows = []
    for line in lines[start:]:
        cols = line.strip().split()
        gt_rows.append(cols[9:])
    vcf_data[cat] = gt_rows
    print(f"Category {cat}: read {len(gt_rows)} sites")

# ========== Initialize the dictionary that stores all resampling results ==========
# Structure: all_counts[pop][freq][cat] = list of values (length N_REPEATS)
all_counts = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

# ========== Perform multiple random resamplings ==========
random.seed(RANDOM_SEED)

freq_range = range(1, TOTAL_ALLELES + 1)  # 1~10

for rep in range(N_REPEATS):
    # Sample SAMPLE_SIZE individual indices for each population
    sub_indices = {}
    for pop in pop_names:
        indices = pop_indices[pop]
        if len(indices) < SAMPLE_SIZE:
            sub = random.choices(indices, k=SAMPLE_SIZE)
        else:
            sub = random.sample(indices, SAMPLE_SIZE)
        sub_indices[pop] = sub

    # Initialize the temporary per-replicate counts (all zeros)
    temp_counts = {pop: {freq: {cat: 0 for cat in categories} for freq in freq_range} for pop in pop_names}

    # Iterate over each category (variant type)
    for cat, gt_rows in vcf_data.items():
        for gts in gt_rows:
            # Count the number of alt alleles for each population at this site
            alt_num = {pop: 0 for pop in pop_names}
            for pop in pop_names:
                alt = 0
                for idx in sub_indices[pop]:
                    gt = gts[idx]
                    if gt.startswith('./.') or gt.startswith('.|.'):
                        continue
                    if gt.startswith('0/0') or gt.startswith('0|0'):
                        continue
                    elif gt.startswith('1/1') or gt.startswith('1|1'):
                        alt += 2
                    elif gt.startswith('0/1') or gt.startswith('0|1') or gt.startswith('1|0'):
                        alt += 1
                    # Other genotypes (e.g. 1/2) are ignored
                alt_num[pop] = alt
            # Add this site to the current replicate's counts
            for pop in pop_names:
                freq = alt_num[pop]
                if freq > 0:  # freq ranges from 1 to 10
                    temp_counts[pop][freq][cat] += 1

    # Store this replicate's results into all_counts
    for pop in pop_names:
        for freq in freq_range:
            for cat in categories:
                all_counts[pop][freq][cat].append(temp_counts[pop][freq][cat])

# ========== Compute means and 95% confidence intervals ==========
rows = []
for pop in pop_names:
    for freq in freq_range:
        daf = freq / TOTAL_ALLELES
        row = {'POP': pop, 'freq': freq, 'DAF': daf}
        for cat in categories:
            vals = all_counts[pop][freq][cat]  # list of length N_REPEATS
            mean_val = np.mean(vals)
            lower_val = np.percentile(vals, 2.5)
            upper_val = np.percentile(vals, 97.5)
            row[f'{cat}_mean'] = mean_val
            row[f'{cat}_lower'] = lower_val
            row[f'{cat}_upper'] = upper_val
        rows.append(row)

df = pd.DataFrame(rows)
df = df.sort_values(['POP', 'freq']).reset_index(drop=True)

# ========== Output results ==========
# Save the complete results.
df.to_csv(output_path, sep="\t", index=False)
print(f"\nFull results saved to: {output_path}")

# Save mean values without lower and upper bounds.
mean_columns = ["POP", "freq", "DAF"] + [
    f"{cat}_mean" for cat in categories
]
df_mean = df[mean_columns].copy()

df_mean.to_csv(mean_output_path, sep="\t", index=False)
print(f"Mean-only results saved to: {mean_output_path}")

print(df_mean.head())
