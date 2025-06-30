import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import cartopy.crs as ccrs
import matplotlib.ticker as ticker
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import matplotlib.transforms as mtransforms
from scipy.stats import kstest, chi2
from scipy.stats import multivariate_normal
from scipy.stats import gaussian_kde

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
        median_xm = np.median(final_df['separation_arcsec'].values)
        print(f"Median separation: {median_xm:.2g} arcsec")
        percentile_95 = np.percentile(final_df['separation_arcsec'].values, 95)
        print(f"95th percentile: {percentile_95:.2g} arcsec")


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

        plt.rcParams.update({'font.size': 10})
        plt.rcParams['pdf.fonttype'] = 42  # For editable text in PDFs
        plt.rcParams['pdf.use14corefonts'] = True
        plt.figure(figsize=(6, 2))
        bins = np.linspace(0, max_sep_arcsec, 51)

        if len(diffs_ambig) > 0:
            plt.hist(diffs_ambig, bins=bins, color='red', edgecolor='k', alpha=0.5, label=f'ambiguous ({len(diffs_ambig)} objects)')

        if len(diffs_nonambig) > 0:
            plt.hist(diffs_nonambig, bins=bins, color='blue', edgecolor='k', alpha=0.5, label=f'non-ambiguous ({len(diffs_nonambig)} objects)')

        plt.xlabel('Separation difference (second closest - closest)')
        if ambiguity_arcsec is not None:
            plt.axvline(ambiguity_arcsec, color='green', linewidth=3, linestyle='dashed', label=f'Ambiguity threshold: {ambiguity_arcsec}"')
        plt.tight_layout()
        plt.legend()
        plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_duplicates.pdf')
        plt.savefig(plot_path)
        plt.show()

    if make_plot:

        # unmatched sources from my_df
        matched_my_indices = final_df[my_ids[0]] if len(my_ids) == 1 else final_df[my_ids].apply(tuple, axis=1)
        if len(my_ids) == 1:
            unmatched_my = my_df[~my_df[my_ids[0]].isin(matched_my_indices)]
        else:
            unmatched_my = my_df[~my_df[my_ids].apply(tuple, axis=1).isin(matched_my_indices)]

        # unmatched sources from your_df
        matched_your_indices = final_df[your_id_new[0]] if len(your_id_new) == 1 else final_df[your_id_new].apply(tuple, axis=1)
        if len(your_id_new) == 1:
            unmatched_your = your_df[~your_df[your_ids[0]].isin(matched_your_indices)]
        else:
            unmatched_your = your_df[~your_df[your_ids].apply(tuple, axis=1).isin(matched_your_indices)]

        # Compute the average Dec in radians for tangent-plane projections (gnomonic projection)
        dec_avg_rad = np.deg2rad(final_df[my_dec])
        # Correct ΔRA for declination
        delta_ra = (final_df[f'RA_{your_label}'] - final_df[my_ra]) * np.cos(dec_avg_rad) * 3600  # arcsec
        delta_dec = (final_df[f'DEC_{your_label}'] - final_df[my_dec]) * 3600  # arcsec
        # Compute statistics
        cov = np.cov(delta_ra, delta_dec)
        sigma_ra = np.sqrt(cov[0, 0])
        sigma_dec = np.sqrt(cov[1, 1])
        corr = cov[0, 1] / (sigma_ra * sigma_dec)

        # Stack the data
        samples = np.vstack([delta_ra, delta_dec])
        kde = gaussian_kde(samples)
        # Reference Gaussian
        sigma_rv = (sigma_ra+ sigma_dec) / 2
        rv = multivariate_normal(mean=[0, 0], cov=[[sigma_rv**2, 0], [0, sigma_rv**2]])
        # Sample grid
        xmin, xmax = -percentile_95, percentile_95
        ymin, ymax = -percentile_95, percentile_95
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        # Evaluate densities
        p = kde(positions)
        q = rv.pdf(positions.T)
        # Avoid divide-by-zero or log(0)
        mask = (p > 0) & (q > 0)
        kl_divergence = np.sum(p[mask] * np.log(p[mask] / q[mask])) * ((xmax - xmin) * (ymax - ymin)) / len(p)
        # print(
        #     f"KL divergence = {kl_divergence:.4f}\n"
        #     "Comparing the empirical 2D distribution of (ΔRA cos(Dec), ΔDec)\n"
        #     "to a zero-mean isotropic 2D Gaussian with σ = "
        #     f"{sigma_rv:.3f}\" in both directions"
        # )

        # plots begin
        plt.rcParams.update({'font.size': 12})
        plt.rcParams['pdf.fonttype'] = 42  # For editable text in PDFs
        plt.rcParams['pdf.use14corefonts'] = True
        fig = plt.figure(figsize=(14, 5))
        gs = gridspec.GridSpec(1, 3, wspace=0.25) # width_ratios=[1, 1.5, 1.5]
        axes = [None, None, None]

        axes[0] = fig.add_subplot(gs[0])
        axes[0].set_title(f"Median: {median_xm:.2g}\"  -  95th perc: {percentile_95:.2g}\"", fontsize=12)
        filtered_plot = d2d[d2d < max_sep * 3].arcsec
        bins = np.logspace(np.log10(filtered_plot.min()), np.log10(filtered_plot.max()), 100)
        axes[0].hist(filtered_plot, bins=bins, log=True, label=f"Matches: {len(final_df)}")
        axes[0].set_ylabel('Counts', labelpad=-2, fontsize=12)
        axes[0].set_xlabel('Distance (arcsec)', fontsize=12)
        axes[0].axvline(max_sep.value, color='red', linestyle='dashed', label=f'Max sep: {max_sep.value:.2g}"')
        axes[0].legend(fontsize=10)
        axes[0].set_xscale('log')
        axes[0].set_xlim(2*1e-3, max_sep_arcsec * 2.5)
        # axes[0].set_ylim(1, 500)
        axes[0].tick_params(labelsize=10)
        axes[0].set_aspect(1.0 / axes[0].get_data_ratio(), adjustable='box')

        axes[1] = fig.add_subplot(gs[1])
        axes[1].hist2d(delta_ra, delta_dec, bins=60, cmap='viridis', norm='log')
        axes[1].set_xlabel(r'$\Delta$RA (arcsec)',fontsize=12)
        axes[1].set_ylabel(r'$\Delta$Dec (arcsec)', labelpad=-2, fontsize=12)
        axes[1].set_aspect('equal', adjustable='box')
        axes[1].axhline(0, color='red', linestyle='--', linewidth=1)
        axes[1].axvline(0, color='red', linestyle='--', linewidth=1)
        axes[1].set_title(f"KL divergence = {kl_divergence:.3f}", fontsize=12)
        axes[1].tick_params(labelsize=10)
        circle = patches.Circle((0, 0), max_sep_arcsec, edgecolor='blue', facecolor='none', linestyle='--', linewidth=2)
        axes[1].add_patch(circle)

        stats_text = (
            r'$\sigma_{\mathrm{RA}}$  = ' + f'{sigma_ra:.3f}"' + "    " + r'$\overline{\Delta\mathrm{RA}}$   = ' + f'{np.mean(delta_ra):.3f}"\n'
            r'$\sigma_{\mathrm{Dec}}$ = ' + f'{sigma_dec:.3f}"' + "    " + r'$\overline{\Delta\mathrm{Dec}}$ = ' + f'{np.mean(delta_dec):.3f}"\n'
            r'$\rho$  = ' + f'{corr:.3f}'
        )
        axes[1].text(
            0.02, 0.98, stats_text,
            transform=axes[1].transAxes,
            ha='left', va='top',
            fontsize=10,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.9, edgecolor='none')
        )

        axes[1].set_xlim(-max_sep_arcsec*1.01, max_sep_arcsec*1.01)
        axes[1].set_ylim(-max_sep_arcsec*1.01, max_sep_arcsec*1.01)

        # --- Coordinate conversion functions ---
        def convert_ra_to_long(ravals):
            """Convert RA [0, 360] to longitude [-180, 180], with flipped sign"""
            longvals = np.atleast_1d(ravals).copy()
            longvals[longvals > 180] -= 360
            return -longvals

        def convert_long_to_ra(longvals):
            """Convert longitude [-180, 180] to RA [0, 360], with flipped sign"""
            ravals = np.atleast_1d(-longvals).copy()
            ravals[ravals < 0] += 360
            return ravals

        central_ra = np.median(final_df[my_ra].values)
        central_dec = np.median(final_df[my_dec].values)
        ra_min, ra_max = final_df[my_ra].values.min(), final_df[my_ra].values.max()
        dec_min, dec_max = final_df[my_dec].values.min(), final_df[my_dec].values.max()
        delta_ra = ra_max - ra_min
        delta_dec = dec_max - dec_min
        ra_min -= 0.05 * delta_ra*np.cos(np.deg2rad(central_dec))
        ra_max += 0.05 * delta_ra*np.cos(np.deg2rad(central_dec))
        dec_min -= 0.05 * delta_dec
        dec_max += 0.05 * delta_dec

        # Compute angular spans (RA already corrected by cos(Dec))
        ra_span_proj = (ra_max - ra_min) * np.cos(np.deg2rad(central_dec))
        dec_span = dec_max - dec_min

        # Enforce aspect ratio = 1 by adjusting the smaller span
        if ra_span_proj > dec_span:
            delta = (ra_span_proj - dec_span) / 2
            dec_min -= delta
            dec_max += delta
        else:
            delta = (dec_span - ra_span_proj) / (2 * np.cos(np.deg2rad(central_dec)))
            ra_min -= delta
            ra_max += delta

        axes[2] = fig.add_subplot(gs[2], projection=ccrs.Gnomonic(central_longitude=convert_ra_to_long(central_ra)[0], central_latitude=central_dec))
        axes[2].set_extent([convert_ra_to_long(ra_min), convert_ra_to_long(ra_max), dec_min, dec_max],crs=ccrs.PlateCarree())
        gl = axes[2].gridlines(
            draw_labels=True,
            xformatter=ticker.FuncFormatter(lambda x, pos: f"{convert_long_to_ra(x)[0]:.1f}"),
            yformatter=ticker.StrMethodFormatter("{x:.1f}"),
            x_inline=False,
            y_inline=False
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.bottom_labels = True
        gl.left_labels = True
        gl.xlabel_style = {'rotation': 0, 'va': 'top', 'fontsize': 10}
        gl.ylabel_style = {'rotation': 0, 'ha': 'right', 'fontsize': 10}

        # --- Plotting ---
        axes[2].scatter(convert_ra_to_long(final_df[my_ra]), final_df[my_dec],
                        transform=ccrs.PlateCarree(), zorder=4,
                        s=3, alpha=1, facecolors='blue', edgecolors='red', linewidths=0, marker='o') # label=f'{my_label} (matched)',

        axes[2].scatter(convert_ra_to_long(final_df[f'RA_{your_label}']), final_df[f'DEC_{your_label}'],
                        transform=ccrs.PlateCarree(), zorder=5,
                        s=1, alpha=1, facecolors='red', edgecolors='blue', linewidths=0, marker='o') #label=f'{your_label} (matched)',

        axes[2].scatter(convert_ra_to_long(unmatched_my[my_ra]), unmatched_my[my_dec],
                        transform=ccrs.PlateCarree(), label=f'{my_label}', zorder=2,
                        s=0.5, alpha=0.6, facecolors='blue', edgecolors='red', linewidths=0, marker='o')

        axes[2].scatter(convert_ra_to_long(unmatched_your[your_ra]), unmatched_your[your_dec],
                        transform=ccrs.PlateCarree(), label=f'{your_label}', zorder=3,
                        s=0.5, alpha=0.6, facecolors='red', edgecolors='blue', linewidths=0, marker='o')

        if draw_lines:
            for i in range(len(final_df)):
                ra1 = convert_ra_to_long(final_df[my_ra].iloc[i])
                ra2 = convert_ra_to_long(final_df[f'RA_{your_label}'].iloc[i])
                dec1 = final_df[my_dec].iloc[i]
                dec2 = final_df[f'DEC_{your_label}'].iloc[i]
                axes[2].plot([ra1, ra2], [dec1, dec2], 'k-', lw=0.1, zorder=6, transform=ccrs.PlateCarree())

        renderer = fig.canvas.get_renderer()
        bbox = axes[2].get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted())
        aspect_ratio = bbox.width / bbox.height
        # print(f"Aspect ratio of the sky plot: {aspect_ratio:.4f}")

        axes[2].legend(
            loc='upper center',
            bbox_to_anchor=(0.5, 1+(0.085*aspect_ratio/0.9457)),
            ncol=2,
            frameon=False,
            handlelength=2.5,
            markerscale=10,
            handletextpad=0.5,
            columnspacing=1.0,
            fontsize=12
        )
        axes[2].text(0.5, 0-0.06*aspect_ratio/0.9457, "RA [deg]", transform=axes[2].transAxes,
                    ha='center', va='top', fontsize=12)
        axes[2].text(-0.11, 0.45, "DEC [deg]", transform=axes[2].transAxes,
                    ha='right', va='center', rotation='vertical',
                    fontsize=12)

        # plt.tight_layout() # not supported by cartopy
        fig.subplots_adjust(left=0, right=1, top=1, bottom=0)
        fig.canvas.draw()
        if len(final_df) > 1e5:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.png')
            plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        else:
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.pdf')
            plt.savefig(plot_path, bbox_inches='tight')
            plot_path = os.path.join(output_dir, f'onexmatch_{my_label}_{your_label}_sep_and_skyplot.png')
            plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        if verbose:
            print("")
            print(f"Plot saved to: {plot_path}")
        plt.show()

    return final_df