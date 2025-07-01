import polars as pl


# data = pl.read_csv('clsa/23ME002_UdeM_SGTaliun_Baseline_CoPv7_Qx_PA_BS.csv')

cat_traits = pl.read_csv('cat_count.txt', has_header = False, separator = " ")


print(cat_traits)

counts = cat_traits.select(
    pl.col("column_2").mean().alias("Mean N Variables"),
    pl.col("column_2").median().alias("Median N Variables"),
    pl.col("column_3").mean().alias("Mean"),
    pl.col("column_3").median().alias("Median"),
    pl.col("column_3").mode().alias("Mode"),
    pl.col("column_3").count().alias("N")
)

print(counts)

counts.write_csv("cat_counts.csv")

dict = pl.read_excel('clsa/23ME002_UdeM_SGTaliun_Baseline_CoPv7_Qx_PA_BS-dictionary.xlsx', sheet_id = 0)

print(dict)
# counts_dict = counts.row(0).to_dict()

# results_df = pl.DataFrame({
#     "column_names": cat_traits["column_names"].to_list(),
#     "non_na_counts": [counts_dict[col] for col in cat_traits["column_names"]]
# })

# results_df.write_csv("cat_counts.csv")