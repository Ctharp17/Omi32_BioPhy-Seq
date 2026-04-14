import pandas as pd

count_tables = []
for sample, fname in zip(snakemake.params.samples, snakemake.input.tsv):
    df = pd.read_table(fname, dtype={'geno': 'string'})
    df['geno'] = df['geno'].astype('string').str.zfill(13)
    count_tables.append(
        df.set_index('geno').rename(columns=lambda _: sample)
    )

(pd.concat(count_tables, axis=1)
  .fillna(0)
  .astype(int)
  .sort_index()
  .to_csv(snakemake.output.tsv, sep='\t'))
