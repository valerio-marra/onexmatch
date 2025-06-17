import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import os

def onexmatch(my_labels, your_labels, max_sep_arcsec=1, ambiguity_arcsec=None, verbose=True, make_plot=True, show_duplicates=True, draw_lines=False):
    """
    Perform a sky coordinate crossmatch between two catalogs.

    Parameters:
        my_labels (dict): Dictionary with keys:
            'file' (str, optional): Path to the first catalog CSV file.
            'df' (pd.DataFrame, optional): DataFrame containing the first catalog.
            'label' (str): Short name for the survey (used in column naming).
            'id' (str or list of str): Name(s) of the unique ID column(s).
            'ra' (str): Name of the Right Ascension column (in degrees).
            'dec' (str): Name of the Declination column (in degrees).
            'extra_columns' (list, optional): Additional columns to include in the output.
            One of 'file' or 'df' must be provided.

        your_labels (dict): Same structure as my_labels, for the second catalog.

        max_sep_arcsec (float): Maximum angular separation for a match (arcseconds).
        ambiguity_arcsec (float, optional): If set, ambiguous matches are discarded.
        verbose (bool): Whether to print match statistics and progress information.
        make_plot (bool): Whether to generate and save diagnostic plots.
        show_duplicates (bool): Whether to display duplicate match information.
        draw_lines (bool): Whether to draw lines between matched sources in the plot.

    Returns:
        final_df (pd.DataFrame): Table of matched sources, including:
            - Unchanged ID, RA, and DEC columns from my_labels
            - Renamed ID, RA, and DEC columns from your_labels (with suffix _{your_label})
            - Any extra columns from either input
            - A 'separation_arcsec' column
            Output is saved as CSV in the same directory as 'file', or in the current working directory
            if DataFrames are used as input.
    """

    # `my_survey` with `my_labels` is the **catalog** you are working on.
    # `your_survey` with `your_labels` is the **source** you want to crossmatch against, typically to retrieve additional properties for `my_survey` entries, such as spectroscopic redshifts or classifications.
    if ambiguity_arcsec is not None and not isinstance(ambiguity_arcsec, (int, float)):
        raise TypeError("ambiguity_arcsec must be a number (float or int).")

    my_label = my_labels['label']
    my_ids = my_labels['id'] if isinstance(my_labels['id'], list) else [my_labels['id']]
    my_ra = my_labels['ra']
    my_dec = my_labels['dec']

    your_label = your_labels['label']
    your_ids = your_labels['id'] if isinstance(your_labels['id'], list) else [your_labels['id']]
    your_ra = your_labels['ra']
    your_dec = your_labels['dec']

    if 'df' in my_labels and 'file' in my_labels:
        raise ValueError("my_labels cannot contain both 'df' and 'file'. Provide only one.")
    elif 'df' in my_labels:
        my_df = my_labels['df']
        my_survey_file = None
    elif 'file' in my_labels:
        my_survey_file = my_labels['file']
        my_df = pd.read_csv(my_survey_file)
    else:
        raise ValueError("my_labels must contain either 'file' or 'df'.")

    if 'df' in your_labels and 'file' in your_labels:
        raise ValueError("your_labels cannot contain both 'df' and 'file'. Provide only one.")
    elif 'df' in your_labels:
        your_df = your_labels['df']
        your_survey_file = None
    elif 'file' in your_labels:
        your_survey_file = your_labels['file']
        your_df = pd.read_csv(your_survey_file)
    else:
        raise ValueError("your_labels must contain either 'file' or 'df'.")

    # to handle possible accidental whitespace in column names from CSVs
    my_df.columns = my_df.columns.str.strip()
    your_df.columns = your_df.columns.str.strip()

    # Check if required columns are present in both DataFrames
    for col in [my_ra, my_dec] + my_ids:
        if col not in my_df.columns:
            raise ValueError(f"Missing column '{col}' in {my_label}")

    for col in [your_ra, your_dec] + your_ids:
        if col not in your_df.columns:
            raise ValueError(f"Missing column '{col}' in {your_label}")

    # your_survey is the source
    ra_source, dec_source = your_df[your_ra].values, your_df[your_dec].values
    # my_survey is the catalog
    ra_cat, dec_cat = my_df[my_ra].values, my_df[my_dec].values

    source = SkyCoord(ra=ra_source * u.deg, dec=dec_source * u.deg)
    cat = SkyCoord(ra=ra_cat * u.deg, dec=dec_cat * u.deg)
    idx, d2d, _ = source.match_to_catalog_sky(cat)

    # for each source entry, find the nearest entry in the catalog
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
    duplicates_my_survey = df_dedu[df_dedu.duplicated(subset=f'{my_label}_idx', keep=False)].sort_values(by=f'{my_label}_idx')


    # define function to identify ambiguous matches
    def get_ambiguous_ids(df, ambiguity_arcsec, my_label):
        grouped = df.groupby(f'{my_label}_idx')
        def is_ambiguous(group):
            if len(group) <= 1:
                return False
            sorted_sep = group['sep'].sort_values().values
            return abs(sorted_sep[1] - sorted_sep[0]) < ambiguity_arcsec
        return [idx for idx, group in grouped if is_ambiguous(group)]

    # optionally filter ambiguous matches
    ambiguous_ids = []
    if ambiguity_arcsec is not None:
        ambiguous_ids = get_ambiguous_ids(df_dedu, ambiguity_arcsec, my_label)
        df_dedu = df_dedu[~df_dedu[f'{my_label}_idx'].isin(ambiguous_ids)]

    # keep only the closest match per my_survey object
    df_unique = df_dedu.sort_values('sep').drop_duplicates(subset=f'{my_label}_idx', keep='first')

    # final matched rows
    # Rename RA and DEC columns to indicate origin
    your_matched = your_df.iloc[df_unique[f'{your_label}_idx']].rename(columns={
        your_ra: f'RA_{your_label}',
        your_dec: f'DEC_{your_label}'
    }).reset_index(drop=True)
    
    # Rename ID columns to indicate origin
    your_id_new = [f'{col}_{your_label}' for col in your_ids]
    your_matched = your_matched.rename(columns=dict(zip(your_ids, your_id_new)))

    my_matched = my_df.iloc[df_unique[f'{my_label}_idx'].values].reset_index(drop=True)

    # Reset index for alignment
    separation = df_unique['sep'].reset_index(drop=True)
    # Concatenate your_survey, my_survey, and separation
    matched_df = pd.concat([your_matched, my_matched, separation.rename('separation_arcsec')], axis=1)

    def safe_columns(df, cols):
        return [c for c in cols if c in df.columns]

    final_df = matched_df[
        my_ids + [my_ra, my_dec] +
        your_id_new + [f'RA_{your_label}', f'DEC_{your_label}', 'separation_arcsec'] +
        safe_columns(matched_df, my_labels.get('extra_columns', [])) +
        safe_columns(matched_df, your_labels.get('extra_columns', []))
    ]
    final_df = final_df.sort_values('separation_arcsec').reset_index(drop=True)

    if my_survey_file is not None:
        output_dir = os.path.dirname(os.path.abspath(my_survey_file))
    else:
        output_dir = os.getcwd()
    output_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}.csv')
    final_df.to_csv(output_path, index=False)
    if final_df.empty:
        print("Warning: No matches found.")

    if verbose:
        print(f"Number of objects in {my_label}: {len(my_df)}")
        print(f"Number of objects in {your_label}: {len(idx)}")
        print(f"Number of matches after imposing the constraint: {len(idx_cat)}")
        if ambiguity_arcsec is not None:
            print(f"Number of matches after filtering ambiguous matches: {len(df_dedu)}")
        print(f"Number of matches after removing duplicates: {len(df_unique)}")
        print(f"Output file saved to: {output_path}")
        print(f"Median separation: {np.median(final_df['separation_arcsec'].values):.2g} arcsec")


    if len(duplicates_my_survey) > 0 and show_duplicates:
        # print(f"Duplicates according to {my_label}:")
        # print(duplicates_my_survey.to_string(index=False))

        is_ambiguous_mask = duplicates_my_survey[f'{my_label}_idx'].isin(ambiguous_ids)
        ambig = duplicates_my_survey[is_ambiguous_mask]
        nonambig = duplicates_my_survey[~is_ambiguous_mask]

        # Group by catalog index and get the two smallest separations for each group
        # Only keep groups with at least two matches
        def compute_diffs(df):
            grouped = df.groupby(f'{my_label}_idx')['sep'].nsmallest(2).reset_index()
            return grouped.groupby(f'{my_label}_idx')['sep'].apply(lambda x: x.iloc[1] - x.iloc[0])

        diffs_ambig = compute_diffs(ambig)
        diffs_nonambig = compute_diffs(nonambig)

        plt.rcParams.update({'font.size': 12})
        plt.rcParams['pdf.fonttype'] = 42  # For editable text in PDFs
        plt.rcParams['pdf.use14corefonts'] = True
        plt.figure(figsize=(6, 3))
        bins = np.linspace(0, max_sep_arcsec, 51)

        if len(diffs_ambig) > 0:
            plt.hist(diffs_ambig, bins=bins, color='red', edgecolor='k', alpha=0.5, label=f'ambiguous ({len(diffs_ambig)})')

        if len(diffs_nonambig) > 0:
            plt.hist(diffs_nonambig, bins=bins, color='blue', edgecolor='k', alpha=0.5, label=f'non-ambiguous ({len(diffs_nonambig)})')

        plt.xlabel('Separation difference (second closest - closest)')
        if ambiguity_arcsec is not None:
            plt.axvline(ambiguity_arcsec, color='red', linestyle='dashed', label=f'Ambiguity threshold: {ambiguity_arcsec} arcsec')
        plt.tight_layout()
        plt.legend()
        plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_duplicates.pdf')
        plt.savefig(plot_path)
        plt.show()

    if make_plot:
        plt.rcParams.update({'font.size': 14})
        plt.rcParams['pdf.fonttype'] = 42  # For editable text in PDFs
        plt.rcParams['pdf.use14corefonts'] = True
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))

        tabplot_galaxies = final_df['separation_arcsec'].values
        axes[0].set_title(f"Median separation: {np.median(tabplot_galaxies):.2g} arcsec — Total matches: {len(tabplot_galaxies)}")
        filtered_plot = d2d[d2d < max_sep * 3].arcsec
        bins = np.logspace(np.log10(filtered_plot.min()), np.log10(filtered_plot.max()), 100)
        axes[0].hist(filtered_plot, bins=bins, log=True)
        axes[0].set_ylabel('Counts')
        axes[0].set_xlabel('Distance (arcsec)')
        axes[0].axvline(max_sep.value, color='red', linestyle='dashed', label=f'Max sep: {max_sep.value:.2g} arcsec')
        axes[0].legend()
        axes[0].set_xscale('log')

        # Plot matched sources and lines between them
        scatter_kwargs = dict(linewidths=0, marker='o')
        axes[1].grid(True)
        
        axes[1].scatter(final_df[my_ra], final_df[my_dec], label=f'{my_label} (matched)',
                        s=8, alpha=1, facecolors='blue', edgecolors='red', zorder=3, **scatter_kwargs)

        axes[1].scatter(final_df[f'RA_{your_label}'], final_df[f'DEC_{your_label}'], label=f'{your_label} (matched)',
                        s=5, alpha=1, facecolors='red', edgecolors='blue', zorder=4, **scatter_kwargs)

        # Draw lines between matched pairs
        if draw_lines:
            for i in range(len(final_df)):
                ra1, dec1 = final_df[my_ra].iloc[i], final_df[my_dec].iloc[i]
                ra2, dec2 = final_df[f'RA_{your_label}'].iloc[i], final_df[f'DEC_{your_label}'].iloc[i]
                axes[1].plot([ra1, ra2], [dec1, dec2], 'k-', lw=0.1, zorder=5)

        # Plot unmatched sources from my_df
        matched_my_indices = final_df[my_ids[0]] if len(my_ids) == 1 else final_df[my_ids].apply(tuple, axis=1)
        if len(my_ids) == 1:
            unmatched_my = my_df[~my_df[my_ids[0]].isin(matched_my_indices)]
        else:
            unmatched_my = my_df[~my_df[my_ids].apply(tuple, axis=1).isin(matched_my_indices)]
        axes[1].scatter(unmatched_my[my_ra], unmatched_my[my_dec], label=f'{my_label} (unmatched)',
                        s=1, alpha=0.6, facecolors='blue', edgecolors='red', zorder=1, **scatter_kwargs)

        # Plot unmatched sources from your_df
        matched_your_indices = final_df[your_id_new[0]] if len(your_id_new) == 1 else final_df[your_id_new].apply(tuple, axis=1)
        if len(your_id_new) == 1:
            unmatched_your = your_df[~your_df[your_ids[0]].isin(matched_your_indices)]
        else:
            unmatched_your = your_df[~your_df[your_ids].apply(tuple, axis=1).isin(matched_your_indices)]
        axes[1].scatter(unmatched_your[your_ra], unmatched_your[your_dec], label=f'{your_label} (unmatched)',
                        s=1, alpha=0.6, facecolors='red', edgecolors='blue', zorder=2, **scatter_kwargs)

        # Final plot settings
        axes[1].set_xlabel('RA [deg]')
        axes[1].set_ylabel('Dec [deg]')
        axes[1].legend(loc='upper center', bbox_to_anchor=(0.5, 1.15), ncol=2, frameon=False)

        # Compute RA and Dec limits from both catalogs
        ra_min = max(my_df[my_ra].min(), your_df[your_ra].min())
        ra_max = min(my_df[my_ra].max(), your_df[your_ra].max())
        dec_min = max(my_df[my_dec].min(), your_df[your_dec].min())
        dec_max = min(my_df[my_dec].max(), your_df[your_dec].max())

        # Apply limits to the plot (RA reversed)
        axes[1].set_xlim(ra_max, ra_min)
        axes[1].set_ylim(dec_min, dec_max)

        plt.tight_layout()
        if len(final_df) > 1e5:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.png')
            plt.savefig(plot_path, dpi=300)
        else:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.pdf')
            plt.savefig(plot_path)
        if verbose:
            print("")
            print(f"Plot saved to: {plot_path}")
        plt.show()

    return final_df