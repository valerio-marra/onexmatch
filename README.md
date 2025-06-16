
# onexmatch

`onexmatch` is a lightweight Python module for crossmatching two astronomical catalogs based on sky coordinates. It is designed for use in cosmological and survey data analysis workflows, such as matching sources between surveys like [J-PAS](https://www.j-pas.org) and [DESI](https://www.desi.lbl.gov).

---

## Features

- Angular crossmatching using `astropy.coordinates.SkyCoord`
- Match filtering by maximum angular separation
- Duplicate resolution by selecting the nearest match
- Customizable input columns (RA, DEC, ID)
- Outputs matched catalog in CSV format
- Generates diagnostic plots (histogram of separations, sky map)

---

## Installation

Clone the repository and install with:

```bash
git clone https://github.com/valerio-marra/onexmatch.git
cd onexmatch
pip install .
```

**Dependencies**:
- numpy
- pandas
- matplotlib
- astropy

Install them with:

```bash
pip install numpy pandas matplotlib astropy
```

---

## Usage


`my_survey` with `my_labels` is the **catalog** you are working on.
`your_survey` with `your_labels` is the **source** you want to crossmatch against, typically to retrieve additional properties for `my_survey` entries, such as spectroscopic redshifts or classifications.
For each object in the **source**, the closest object in the **catalog** is identified based on sky coordinates.
However, the same **catalog** object may be matched to multiple **source** objects.
In such cases, `onexmatch` retains only the closest **source** match for each **catalog** object.


Here's a basic usage example using synthetic data. A full example is available in `example/example.ipynb`.

```python
from onexmatch import onexmatch

matched_df = onexmatch(
    my_labels={
        # Provide either a file path:
        'file': 'cats/J-PAS_synthetic.csv',
        # or directly a pandas DataFrame:
        # 'df': my_df,
        'label': 'J-PAS',
        'id': 'TILE-NUMBER',
        'ra': 'RA',
        'dec': 'DEC'
    },
    your_labels={
        'file': 'cats/DESI_synthetic.csv',
        # 'df': your_df,
        'label': 'DESI',
        'id': 'TILE-ID',
        'ra': 'RA',
        'dec': 'DEC'
    },
    max_sep_arcsec=1.0,
    verbose=True,
    make_plot=True
)
```

This will:

- Match all sources in `your_labels` against `my_labels` within 1.0 arcseconds
- Save a CSV file: `onexmatch_J-PAS_DESI.csv`
- Generate a plot file: `onexmatch_J-PAS_DESI_sep_and_skyplot.pdf` or `.png` (depending on the size of the xmatch)

If file paths are used, outputs are saved in the same directory as `my_labels['file']`.  
If DataFrames are used, outputs are saved in the current working directory.

---

## Output Columns

The resulting DataFrame and CSV include 7 columns:

- RA and DEC columns from both catalogs, renamed as `RA_{label}` and `DEC_{label}`
- ID columns (renamed to `ID_{label}` only if both use the same original name)
- `separation_arcsec`: the angular separation between matched pairs

---


## License

MIT License

---

## Author

[Valerio Marra](http://marra.cosmo-ufes.org)  
[valerio.marra@me.com](mailto:valerio.marra@me.com)  
Federal University of Espírito Santo (UFES), Brazil
