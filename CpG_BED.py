import os
import pandas as pd
from pybedtools import BedTool
from datatable import dt, f, by, g, join, sort, update, ifelse

def find_overlapping_regions(bed1_df, bed2_df, bed1_cols, bed2_cols):
    bed_1 = BedTool.from_dataframe(bed1_df[bed1_cols].sort_values(bed1_cols))
    bed_2 = BedTool.from_dataframe(bed2_df[bed2_cols].sort_values(bed2_cols))
    return BedTool.to_dataframe(bed_1.intersect(bed_2, wa=True, wb=True))
    
def CpGs_to_features(met, ref, metcols, refcols, loffset=0, roffset=1, loffset2=0, roffset2=1):
    # met: datatable
    # ref: pandas
    # metcols: [chrom, pos, other columns that are not samples]
    # refcols: [chrom, start, end, name]
    cpgpos = met[:, [metcols[0], metcols[1]]].to_pandas()
    cpgpos['start'] = cpgpos[metcols[1]]+loffset
    cpgpos['end'] = cpgpos[metcols[1]]+roffset

    ref2 = ref.copy()
    ref2[refcols[1]] = ref2[refcols[1]] + loffset2
    ref2[refcols[2]] = ref2[refcols[2]] + roffset2
    cpgsoi = find_overlapping_regions(ref2, cpgpos, refcols, [metcols[0],'start','end'])
    cpgsoi.index = cpgsoi['score'].astype(str)+'@'+(cpgsoi['strand']-loffset).astype(str)

    cpgpos['name'] = cpgpos[metcols[0]].astype(str)+'@'+cpgpos['start'].astype(str)
    cpgpos['name'] = cpgpos['name'].map(cpgsoi['name'].to_dict())

    met_features = met.copy()
    del met_features[: , metcols]
    met_features['name'] = dt.Frame(cpgpos['name'])
    met_features = met_features[:, dt.mean(f[:]), by('name')]

    return met_features
