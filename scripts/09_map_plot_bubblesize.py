import os
import struct
import binascii
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap
import geohash

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------
csv_path = "../data/full_card_metadata_aro_allfilters_metagenomes_WHOcategories_latlon.csv"

base_col = {
    "livestock":   "#fac723",
    "human":       "#936fac",
    "wastewater":  "#e95e50",
    "marine":      "#0cb2af",
    "freshwater":  "#a1c65d",
    "soil":        "#f29222",
}
keep = set(base_col)

usecols = ["WTK", "metagenome_category", "acc"]    # read just 3 columns to reduce memory usage
dtype   = {
	"metagenome_category": "category",   # 1 byte per row instead of a long str
	"acc": "object",                     # keep accessions as Python strings
}

out_dir = "../data/map_plots"
os.makedirs(out_dir, exist_ok=True)

# bubble-size parameters (your original linear formula)
min_size = 25              # pt²
scale    = 0.3            # additional pt² per extra sample

gh_prec = 2               # geohash precision (3 = ~150km, 4 = ~40km, 5 = ~5km)
# -------------------------------------------------------------------

## Functions

# Translate WTK to lon/lat
def wkb_to_lonlat(hexstr: str):
    """
    Decode a WKB POINT (hex) with or without SRID (EPSG:4326) → (lon, lat).

    Handles both:
        • 25-byte variant 0101… (little endian) or 0001… (big endian) with SRID
        • 21-byte variant without SRID
    Returns (np.nan, np.nan) if anything looks wrong.
    """
    try:
        b = binascii.unhexlify(hexstr.strip())
    except (binascii.Error, AttributeError):
        return np.nan, np.nan

    if len(b) < 21:         # too short to hold 2 doubles
        return np.nan, np.nan

    endian = "<" if b[0] == 1 else ">"
    try:
        if len(b) == 25:    # SRID header present
            x = struct.unpack(endian + "d", b[9:17])[0]
            y = struct.unpack(endian + "d", b[17:25])[0]
        else:               # plain WKB point
            x = struct.unpack(endian + "d", b[5:13])[0]
            y = struct.unpack(endian + "d", b[13:21])[0]
    except struct.error:
        return np.nan, np.nan
    return x, y

def main():
    df = pd.read_csv(csv_path, usecols=usecols, dtype=dtype, low_memory=False)

    df = df.drop_duplicates(subset="acc", keep="first")

    df = df[df["metagenome_category"].isin(keep)].copy()

    lonlat = df["WTK"].apply(wkb_to_lonlat)

    # unpack the tuples into two numeric columns
    df["lon"] = lonlat.apply(lambda t: t[0]).astype("float32")
    df["lat"] = lonlat.apply(lambda t: t[1]).astype("float32") 

    # drop rows where decoding failed
    df = df.dropna(subset=["lon", "lat"])

    # Geohash cells
    df["geo"] = [geohash.encode(lat, lon, precision=gh_prec)
                 for lat, lon in zip(df["lat"], df["lon"])]

    grouped = (
        df.groupby(["geo", "metagenome_category"], observed=True)
          .agg(lon=("lon", "mean"),
               lat=("lat", "mean"),
               n_samples=("acc", "size"))
          .reset_index()
    )

    # use the cell centre as the plotting coordinate
    grouped["size"] = min_size + grouped["n_samples"] * scale
    grouped["colour"] = grouped["metagenome_category"].map(base_col)

    fig, ax = plt.subplots(figsize=(13, 7))
    m = Basemap(projection="moll", lon_0=0, resolution="c")
    m.drawcoastlines(linewidth=0.4)
    m.fillcontinents(color="#ededed", lake_color="white")
    #m.drawparallels(np.arange(-90, 91, 30))
    #m.drawmeridians(np.arange(0, 361, 60))
    m.drawmapboundary(fill_color="white")

    for cat in base_col:
        sub = grouped[grouped["metagenome_category"] == cat]
        if sub.empty:
            continue
        x, y = m(sub["lon"].values, sub["lat"].values)
        m.scatter(
            x, y,
            s=sub["size"].values,
            c=sub["colour"].values,
            edgecolors="face",
            linewidths=0.2,
            alpha=0.6,
            zorder=5,
        )    

	# legend 1 : colour
    colour_handles = [
		mpatches.Patch(facecolor=base_col[c], edgecolor='none', label=c)
		for c in base_col
	]
    leg_colour = ax.legend(
		handles=colour_handles,
		title="Metagenome category",
		loc="lower left",
		fontsize="small",
		frameon=False,
	)
    ax.add_artist(leg_colour)

	# legend 2 : bubble size
    unique_counts = np.unique(grouped["n_samples"])
    if len(unique_counts) >= 3:
        counts_for_legend = [unique_counts.min(),
							 np.mean(unique_counts),
							 unique_counts.max()]
    else:
        counts_for_legend = unique_counts

    sizes_for_legend = min_size + np.array(counts_for_legend) * scale
    size_handles = [
		plt.scatter([], [], s=s, marker='o', color='grey',
					edgecolors='face', alpha=0.5,
					label=f"{int(c)}")
		for s, c in zip(sizes_for_legend, counts_for_legend)
	]
    leg_size = ax.legend(
		handles=size_handles,
		title="Unique accessions",
		loc="lower right",
		scatterpoints=1,
		fontsize="small",
		frameon=False,
	)
    ax.add_artist(leg_size)

    plt.title("AMR positive unique accessions by metagenome category")
    out = os.path.join(out_dir, f"map_plot_sizebubbles_geohash.png")
    fig.savefig(out, dpi=600)
    fig.savefig(out.replace(".png", ".svg"))


if __name__ == "__main__":
    main()
