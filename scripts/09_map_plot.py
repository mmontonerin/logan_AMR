import os
import struct
import binascii
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.basemap import Basemap

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
point_size = 12

out_dir = "../data/map_plots"
os.makedirs(out_dir, exist_ok=True)
# -------------------------------------------------------------------


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
    if not os.path.isfile(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path!r}")

    df = pd.read_csv(csv_path, low_memory=False, dtype=str)

    df = df[df["metagenome_category"].isin(keep)].copy()

    lonlat = df["WTK"].apply(wkb_to_lonlat)
    df["lon"] = lonlat.str[0]
    df["lat"] = lonlat.str[1]
    df = df.dropna(subset=["lon", "lat"])

    # one dot per unique positive accession
    df = df.sort_values("acc").drop_duplicates("acc")

    df["colour"] = df["metagenome_category"].map(base_col)

    plt.figure(figsize=(13, 7))
    m = Basemap(projection="moll", lon_0=0, resolution="c")
    m.drawcoastlines(linewidth=0.4)
    m.fillcontinents(color="lightgrey", lake_color="white")
    m.drawparallels(np.arange(-90, 91, 30))
    m.drawmeridians(np.arange(0, 361, 60))
    m.drawmapboundary(fill_color="white")

    x, y = m(df["lon"].values, df["lat"].values)
    m.scatter(
        x,
        y,
        s=point_size,
        c=df["colour"].values,
        edgecolors="k",
        linewidths=0.15,
        alpha=0.85,
        zorder=5,
    )

    legend_handles = [
        mpatches.Patch(color=base_col[c], label=c) for c in keep
    ]
    plt.legend(
        handles=legend_handles,
        loc="lower left",
        fontsize="small",
        title="Metagenome category",
    )

    plt.title("Positive unique accessions by metagenome category (Mollweide projection)")
    out = os.path.join(out_dir, f"map_plot.png")
    plt.savefig(out, dpi=600)
    plt.savefig(out.replace(".png", ".svg"))


if __name__ == "__main__":
    main()
