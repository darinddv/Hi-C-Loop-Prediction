from intervaltree import Interval, IntervalTree
import pandas as pd

def find_overlaps(motifs_df, chr_col_motifs, start_col_motifs, end_col_motifs,
                  search_df, chr_col_search, start_col_search, end_col_search):
    
    # Ensure chromosome columns are strings
    motifs_df[chr_col_motifs] = motifs_df[chr_col_motifs].astype(str)
    search_df[chr_col_search] = search_df[chr_col_search].astype(str)
    
    chromosomes = [str(i) for i in range(1, 25)]
    overlaps = []

    for chrom in chromosomes:
        motifs_chrom = motifs_df[motifs_df[chr_col_motifs] == chrom]
        search_chrom = search_df[search_df[chr_col_search] == chrom]

        tree = IntervalTree(Interval(row[start_col_motifs], row[end_col_motifs], index) for index, row in motifs_chrom.iterrows())

        for index, row in search_chrom.iterrows():
            overlapping_intervals = tree[row[start_col_search]:row[end_col_search]]
            overlaps.extend((index, interval.data) for interval in overlapping_intervals)

    if not overlaps:
        # Return an empty DataFrame with the appropriate columns if no overlaps are found
        print("No overlaps found")
        return pd.DataFrame(columns=list(search_df.columns) + list(motifs_df.columns))

    overlapping_regions = pd.DataFrame([(search_df.loc[i1], motifs_df.loc[i2]) for i1, i2 in overlaps], 
                                       columns=['search_df', 'motifs_df'])

    result = pd.concat([overlapping_regions['search_df'].apply(pd.Series).set_axis(search_df.columns, axis=1),
                        overlapping_regions['motifs_df'].apply(pd.Series).set_axis(motifs_df.columns, axis=1)], axis=1)
    
    return result



