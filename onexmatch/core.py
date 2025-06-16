import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import os

def onexmatch(my_labels, your_labels, max_sep_arcsec=1, verbose=True, make_plot=True):
    """
    Perform a sky coordinate crossmatch between two catalogs.

    Parameters:
        my_labels (dict): Dictionary with keys:
            'file' (str, optional): Path to the first catalog CSV file.
            'df' (pd.DataFrame, optional): DataFrame containing the first catalog.
            'label' (str): Short name for the survey (used in column renaming).
            'id' (str): Name of the unique ID column.
            'ra' (str): Name of the Right Ascension column (in degrees).
            'dec' (str): Name of the Declination column (in degrees).
            One of 'file' or 'df' must be provided.

        your_labels (dict): Same structure as my_labels, for the second catalog.

        max_sep_arcsec (float): Maximum angular separation for a match (arcseconds).
        verbose (bool): Whether to print match statistics and duplicate entries.
        make_plot (bool): Whether to generate and save diagnostic plots.

    Returns:
        matched_df (pd.DataFrame): Table of matched sources, including renamed RA, DEC, ID columns,
                                   and a separation column in arcseconds. Output is saved as CSV in
                                   the same directory as 'file', or in the current working directory
                                   if DataFrames are used as input.
    """

    # `my_survey` with `my_labels` is the **catalog** you are working on.
    # `your_survey` with `your_labels` is the **source** you want to crossmatch against, typically to retrieve additional properties for `my_survey` entries, such as spectroscopic redshifts or classifications.

    my_label, my_id, my_ra, my_dec = (
        my_labels['label'], my_labels['id'], my_labels['ra'], my_labels['dec']
    )
    your_label, your_id, your_ra, your_dec = (
        your_labels['label'], your_labels['id'], your_labels['ra'], your_labels['dec']
    )

    if 'df' in my_labels:
        my_df = my_labels['df']
        my_survey_file = None
    elif 'file' in my_labels:
        my_survey_file = my_labels['file']
        my_df = pd.read_csv(my_survey_file)
    else:
        raise ValueError("my_labels must contain either 'file' or 'df'.")

    if 'df' in your_labels:
        your_df = your_labels['df']
        your_survey_file = None
    elif 'file' in your_labels:
        your_survey_file = your_labels['file']
        your_df = pd.read_csv(your_survey_file)
    else:
        raise ValueError("your_labels must contain either 'file' or 'df'.")

    for col in [my_ra, my_dec, my_id]:
        if col not in my_df.columns:
            raise ValueError(f"Missing column '{col}' in {my_label}")

    for col in [your_ra, your_dec, your_id]:
        if col not in your_df.columns:
            raise ValueError(f"Missing column '{col}' in {your_label}")

    # your_survey is the source
    ra_source, dec_source = your_df[your_ra].values, your_df[your_dec].values
    # my_survey is the catalog
    ra_cat, dec_cat = my_df[my_ra].values, my_df[my_dec].values

    source = SkyCoord(ra=ra_source * u.deg, dec=dec_source * u.deg)
    cat = SkyCoord(ra=ra_cat * u.deg, dec=dec_cat * u.deg)
    idx, d2d, _ = source.match_to_catalog_sky(cat)

    # for each source object, one obtains the nearest cat object
    # idx and d2d have length equal to len(source)
    # idx[i] gives the index of the closest cat object to source[i].

    max_sep = max_sep_arcsec * u.arcsec
    sep_constraint = d2d < max_sep

    source_matches_indices = np.where(sep_constraint)[0]
    idx_cat = idx[sep_constraint]
    d2d_cat = d2d[sep_constraint]

    # Create a DataFrame for matched sources
    # the same my_survey object could be matched to two (or more) your_survey objects
    df_dedu = pd.DataFrame({
        f'{your_label}_idx': source_matches_indices,
        f'{my_label}_idx': idx_cat,
        'sep': d2d_cat.arcsec
    })

    # subset with duplicates according to my_survey
    duplicates_my_survey = df_dedu[df_dedu.duplicated(subset=f'{my_label}_idx', keep=False)]
    # keep only the closest match per my_survey object
    df_unique = df_dedu.sort_values('sep').drop_duplicates(subset=f'{my_label}_idx', keep='first')

    if my_id == your_id:
        my_id_new = 'ID_' + my_label
        your_id_new = 'ID_' + your_label
    else:
        my_id_new = my_id
        your_id_new = your_id

    # final matched rows
    # Rename RA and DEC columns to indicate origin
    your_matched = your_df.iloc[df_unique[f'{your_label}_idx'].values].rename(columns={
        your_ra: f'RA_{your_label}',
        your_dec: f'DEC_{your_label}',
        your_id: your_id_new
    }).reset_index(drop=True)

    my_matched = my_df.iloc[df_unique[f'{my_label}_idx'].values].rename(columns={
        my_ra: f'RA_{my_label}',
        my_dec: f'DEC_{my_label}',
        my_id: my_id_new
    }).reset_index(drop=True)

    # Reset index for alignment
    separation = df_unique['sep'].reset_index(drop=True)
    # Concatenate your_survey, my_survey, and separation
    matched_df = pd.concat([your_matched, my_matched, separation.rename('separation_arcsec')], axis=1)

    if my_survey_file is not None:
        output_dir = os.path.dirname(os.path.abspath(my_survey_file))
    else:
        output_dir = os.getcwd()
    output_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}.csv')
    matched_df.to_csv(output_path, index=False)

    if verbose:
        print(f"Number of objects in {my_label}: {len(my_df)}")
        print(f"Number of objects in {your_label}: {len(idx)}")
        print(f"Number of matches after imposing the constraint: {len(idx_cat)}")
        print(f"Number of matches after removing duplicates: {len(df_unique)}")
        print(f"Output file saved to: {output_path}")
        print(f"Median separation: {np.median(matched_df['separation_arcsec'].values):.2g} arcsec")
        print("")
        print(f"Duplicates according to {my_label}:")
        print(duplicates_my_survey.to_string(index=False))

    if make_plot:
        plt.rcParams.update({'font.size': 14})
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        tabplot_galaxies = matched_df['separation_arcsec'].values
        axes[0].set_title(f"Median separation: {np.median(tabplot_galaxies):.2g} arcsec — Total matches: {len(tabplot_galaxies)}")
        filtered_plot = d2d[d2d < max_sep * 3].arcsec
        bins = np.logspace(np.log10(filtered_plot.min()), np.log10(filtered_plot.max()), 100)
        axes[0].hist(filtered_plot, bins=bins, log=True)
        axes[0].set_ylabel('Counts')
        axes[0].set_xlabel('Distance (arcsec)')
        axes[0].axvline(max_sep.value, color='red', linestyle='dashed', label=f'Max sep: {max_sep.value:.2g} arcsec')
        axes[0].legend()
        axes[0].set_xscale('log')

        axes[1].scatter(matched_df[f'RA_{my_label}'], matched_df[f'DEC_{my_label}'], s=50, label=f'{my_label} (matched)', alpha=0.5)
        axes[1].scatter(matched_df[f'RA_{your_label}'], matched_df[f'DEC_{your_label}'], s=40, label=f'{your_label} (matched)', color='red', marker='x')
        axes[1].set_xlabel('RA [deg]')
        axes[1].set_ylabel('Dec [deg]')
        axes[1].set_title(f'Matched {your_label} and {my_label} Sources')
        axes[1].legend()
        axes[1].grid(True)
        axes[1].invert_xaxis()

        plt.tight_layout()
        if len(matched_df) >1e3:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.png')
            plt.savefig(plot_path, dpi=400)
        else:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.pdf')
            plt.savefig(plot_path)
        if verbose:
            print("")
            print(f"Plot saved to: {plot_path}")
        plt.show()

    return matched_df