import collections

# Initiate dictionaries for each ancestry loci file
aa_ancestry_counts = collections.defaultdict(lambda: [0, 0, 0])
nhw_ancestry_counts = collections.defaultdict(lambda: [0, 0, 0])
hispanic_ancestry_counts = collections.defaultdict(lambda: [0, 0, 0])

# Read the file with the AA MAGENTA samples and their ancestries.
with open("AA_individual_ancestry_loci.txt") as f:
    for i, line in enumerate(f):
        sample_id, allele_id, anc1, anc2, gt = line.strip().split()
        aa_ancestry_counts[sample_id][int(anc1)] += 1
        aa_ancestry_counts[sample_id][int(anc2)] += 1

# Read the file with the NHW MAGENTA samples and their ancestries.
with open("NHW_individual_ancestry_loci.txt") as f:
    for i, line in enumerate(f):
        sample_id, allele_id, anc1, anc2, gt = line.strip().split()
        nhw_ancestry_counts[sample_id][int(anc1)] += 1
        nhw_ancestry_counts[sample_id][int(anc2)] += 1

# Read the file with the HISPANIC MAGENTA samples and their ancestries.
with open("HISPANIC_individual_ancestry_loci.txt") as f:
    for i, line in enumerate(f):
        sample_id, allele_id, anc1, anc2, gt = line.strip().split()
        hispanic_ancestry_counts[sample_id][int(anc1)] += 1
        hispanic_ancestry_counts[sample_id][int(anc2)] += 1

# Calculate the proportions for each ancestry loci file
aa_ancestry_proportions = collections.defaultdict(lambda: [0, 0, 0])
nhw_ancestry_proportions = collections.defaultdict(lambda: [0, 0, 0])
hispanic_ancestry_proportions = collections.defaultdict(lambda: [0, 0, 0])

# Calculate the proportions for AA
for sample_id, counts in aa_ancestry_counts.items():
    total = sum(counts)
    aa_ancestry_proportions[sample_id] = [count / total for count in counts]

# Calculate the proportions for NHW
for sample_id, counts in nhw_ancestry_counts.items():
    total = sum(counts)
    nhw_ancestry_proportions[sample_id] = [count / total for count in counts]

# Calculate the proportions for HISPANIC
for sample_id, counts in hispanic_ancestry_counts.items():
    total = sum(counts)
    hispanic_ancestry_proportions[sample_id] = [count / total for count in counts]

# Write the counts to a tab delimited txt file for AA
with open("AA_ancestry_counts.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, counts in aa_ancestry_counts.items():
        f.write(f"{sample_id}\t{counts[0]}\t{counts[1]}\t{counts[2]}\n")

# Write the proportions to a tab delimited txt file for AA
with open("AA_ancestry_proportions.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, proportions in aa_ancestry_proportions.items():
        f.write(f"{sample_id}\t{proportions[0]}\t{proportions[1]}\t{proportions[2]}\n")

# Write the counts to a tab delimited txt file for NHW
with open("NHW_ancestry_counts.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, counts in nhw_ancestry_counts.items():
        f.write(f"{sample_id}\t{counts[0]}\t{counts[1]}\t{counts[2]}\n")

# Write the proportions to a tab delimited txt file for NHW
with open("NHW_ancestry_proportions.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, proportions in nhw_ancestry_proportions.items():
        f.write(f"{sample_id}\t{proportions[0]}\t{proportions[1]}\t{proportions[2]}\n")

# Write the counts to a tab delimited txt file for HISPANIC
with open("HISPANIC_ancestry_counts.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, counts in hispanic_ancestry_counts.items():
        f.write(f"{sample_id}\t{counts[0]}\t{counts[1]}\t{counts[2]}\n")

# Write the proportions to a tab delimited txt file for HISPANIC
with open("HISPANIC_ancestry_proportions.txt", "w") as f:
    f.write("sample_ID\tPEL\tCEU\tYRI\n")
    for sample_id, proportions in hispanic_ancestry_proportions.items():
        f.write(f"{sample_id}\t{proportions[0]}\t{proportions[1]}\t{proportions[2]}\n")
