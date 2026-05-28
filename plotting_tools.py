import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rtm_constants as CONST
import seaborn as sns
import matplotlib.lines as mlines
from matplotlib.lines import Line2D


def _coerce_bool(series):
    if series.dtype == bool:
        return series
    values = series.astype(str).str.strip().str.upper()
    true_set = {'TRUE', 'T', '1', 'YES', 'Y'}
    false_set = {'FALSE', 'F', '0', 'NO', 'N'}
    return values.map(lambda v: True if v in true_set else False if v in false_set else np.nan)


def _baseline_only(df, column):
    return df.loc[_coerce_bool(df[column]) == True].copy()



def analyze_placebo_response(cohort_df,eligibility_min=0,figOut=True):
    # Calculate placebo response statistics
    
    """
    Create a heatmap plotting mSF on the x-axis and PC on the y-axis.
    
    Parameters:
        cohort_df: DataFrame with 'mSF' and 'PC' columns
        bins: number of bins for the 2D histogram
    """

    median_PC = np.nanmean(cohort_df['PC'])
    mean_PC = np.nanmean(cohort_df['PC'])
    median_PCstd = np.nanstd(cohort_df['PC'])
    sem_PC = median_PCstd / np.sqrt(len(cohort_df))

    rr50 = 100*len(cohort_df[cohort_df['PC'] >= 50]) / len(cohort_df)

    subsetLoM = cohort_df[cohort_df['mSF'] < eligibility_min]
    mean_PC_subsetLoM = np.nanmean(subsetLoM['PC']) if len(subsetLoM) > 0 else np.nan
    subsetHiM = cohort_df[cohort_df['mSF'] >= eligibility_min]
    mean_PC_subsetHiM = np.nanmean(subsetHiM['PC']) if len(subsetHiM) > 0 else np.nan

    subsetLoB = cohort_df[cohort_df['diaryBASELINE'] < eligibility_min]
    mean_PC_subsetLoB = np.nanmean(subsetLoB['PC']) if len(subsetLoB) > 0 else np.nan  # Fixed typo: nanman -> nanmean
    subsetHiB = cohort_df[cohort_df['diaryBASELINE'] >= eligibility_min]
    mean_PC_subsetHiB = np.nanmean(subsetHiB['PC']) if len(subsetHiB) > 0 else np.nan
    print(f'{mean_PC:5.3f},{sem_PC:5.3f},{rr50:5.3f},{mean_PC_subsetLoM:5.3f},{len(subsetLoM)/len(cohort_df):5.3f},{mean_PC_subsetHiM:5.3f},{len(subsetHiM)/len(cohort_df):5.3f},{mean_PC_subsetLoB:5.3f},{len(subsetLoB)/len(cohort_df):5.3f},{mean_PC_subsetHiB:5.3f},{len(subsetHiB)/len(cohort_df):5.3f}')
    with open(CONST.OUTPUT_FILENAME, 'a') as f:
        f.write(f'{mean_PC},{sem_PC},{rr50},{mean_PC_subsetLoM},{len(subsetLoM)/len(cohort_df)},{mean_PC_subsetHiM},{len(subsetHiM)/len(cohort_df)},{mean_PC_subsetLoB},{len(subsetLoB)/len(cohort_df)},{mean_PC_subsetHiB},{len(subsetHiB)/len(cohort_df)}\n')
    # Flatten arrays if needed (since mSF and PC may be stored as arrays)
    mSF_values = cohort_df['mSF'].apply(lambda x: x[0] if hasattr(x, '__len__') else x)
    PC_values = cohort_df['PC'].apply(lambda x: x[0] if hasattr(x, '__len__') else x)
    baseline_values = cohort_df['diaryBASELINE'].apply(lambda x: x[0] if hasattr(x, '__len__') else x)
    if figOut:
        plt.figure(figsize=(8, 3))
        plt.subplot(1,2,1)
        plot_heatmap(mSF_values, PC_values, eligibility_min, 'True')
        plt.subplot(1,2,2)
        plot_heatmap(baseline_values, PC_values, eligibility_min,'Observed')
        plt.show()
        

def plot_heatmap(mSF_values, PC_values, eligibility_min, title_prefix):
    bins=CONST.HEATMAP_BINS
    theRange=[[0, 15], [-100, 100]]
    theExtent=[theRange[0][0], theRange[0][1], theRange[1][0], theRange[1][1]]
    y_bins = np.arange(-100, 110, 10)

    
    # Create 2D histogram as heatmap with fixed range
    heatmap, _, _ = np.histogram2d(
        mSF_values, PC_values,
        bins=[bins, y_bins],
        range=theRange, 
        density=True
    )

    # Plot using imshow
    
    plt.imshow(heatmap.T, origin='lower', aspect='auto',
            extent=theExtent,
            cmap='viridis')

    plt.axvline(x=eligibility_min, color='red', linewidth=3)
    plt.colorbar(label='Count')
    plt.xlabel('Monthly Seizure Frequency (mSF)')
    plt.ylabel('Percent Change (PC)')
    plt.title(f'Heatmap: PC vs {title_prefix} mSF')
    
    plt.tight_layout()
    #plt.show()
    
    #return plt.gcf()

def make_FAR_plot():

    # -*- coding: utf-8 -*-
    """
    MPC vs MIN by FAR using output.csv
    - Sensitivity == 1.00 only
    - BaselineTF=TRUE only (eligibility tested during baseline)
    - FAR=0.0 lines highlighted (thicker)
    - Title: 'Effect of changing FAR (baseline eligibility)'
    """


    # --- Load data from CSV ---
    csv_path = CONST.OUTPUT_FILENAME
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip() for c in df.columns}, inplace=True)

    # --- Validate required columns ---
    required_cols = {'Sensitivity', 'FAR', 'BaselineTF', 'MIN'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")  # [1](https://bilh-my.sharepoint.com/personal/dgoldenh_bidmc_harvard_edu/Documents/Microsoft%20Copilot%20Chat%20Files/output.csv)
    mpc_col = 'MPCmean' if 'MPCmean' in df.columns else 'MPC'
    if mpc_col not in df.columns:
        raise ValueError(f"Missing MPC/MPCmean column. Found columns: {list(df.columns)}")

    # --- Filter: Sensitivity == 1.00 ---
    baseline_df = _baseline_only(df, 'BaselineTF')
    sens_mask = pd.to_numeric(baseline_df['Sensitivity'], errors='coerce').eq(1.0)
    true_df = baseline_df.loc[sens_mask, ['FAR', 'MIN', mpc_col]].copy()
    true_df.rename(columns={mpc_col: 'MPC'}, inplace=True)

    # --- Aggregate to one MPC per (FAR, MIN) to avoid duplicate points ---
    true_agg  = (true_df.groupby(['FAR','MIN'], as_index=False)
                .agg(MPC=('MPC','mean')).sort_values(['FAR','MIN']))

    # --- Color map per FAR (same color for TRUE/FALSE of same FAR) ---
    all_fars = sorted(set(true_agg['FAR']))
    cmap = plt.get_cmap('tab10')
    color_map = {far: cmap(i % 10) for i, far in enumerate(all_fars)}

    # --- Line width rules (highlight FAR=0.0) ---
    base_lw = 2.0
    highlight_lw = 4.0

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10.5, 6.2))

    # Baseline-only eligibility: solid + closed markers
    for far_val, grp in true_agg.groupby('FAR'):
        color = color_map[far_val]
        lw = highlight_lw if float(far_val) == 0.0 else base_lw
        ax.plot(grp['MIN'], grp['MPC'],
                linestyle='-', linewidth=lw, color=color,
                marker='o', markersize=5,
                markerfacecolor=color, markeredgecolor=color)

    # --- Axes labels and title ---
    ax.set_title("Effect of changing FAR (baseline eligibility)")
    ax.set_xlabel('MIN')
    ax.set_ylabel('MPC')
    ax.grid(True, linestyle='--', alpha=0.35)

    # --- Legend proxies ---
    true_proxy = [Line2D([0],[0], color=color_map[far], linestyle='-',
                        marker='o', markersize=6,
                        markerfacecolor=color_map[far], markeredgecolor=color_map[far],
                        linewidth=(highlight_lw if float(far)==0.0 else base_lw))
                for far in all_fars if far in set(true_agg['FAR'])]

    true_labels  = [f'FAR={far}' for far in all_fars if far in set(true_agg['FAR'])]

    # Set x-axis to show only integer tick values
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Place legend OUTSIDE to the right
    leg = ax.legend(true_proxy, true_labels,
                    title="Eligibility tested during baseline",
                    frameon=True, facecolor='white', framealpha=0.90,
                    handlelength=3.0,
                    loc='center left', bbox_to_anchor=(1.02, 0.5))

    # Make room on the right for the legend
    plt.subplots_adjust(right=0.75)

    fig.tight_layout()
    plt.show()

def plot_all_combos_of_meanMPC_vs_MIN():
    """Plot MPCmean vs MIN for rows where eligibility is tested during baseline."""
    df = pd.read_csv(CONST.OUTPUT_FILENAME)
    df = _baseline_only(df, 'BaselineTF')

    unique_sensitivities = sorted(df['Sensitivity'].unique())
    unique_fars = sorted(df['FAR'].unique())
    if len(unique_fars) != 4:
        raise ValueError(f"Expected 4 FAR values for a 2x2 plot. Found {len(unique_fars)}: {unique_fars}")

    sens_colors = plt.cm.viridis(np.linspace(0, 0.9, len(unique_sensitivities)))
    sens_color_map = {sens: sens_colors[i] for i, sens in enumerate(unique_sensitivities)}

    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True, sharey=True)
    axes = axes.flat

    for ax, far in zip(axes, unique_fars):
        far_df = df[df['FAR'] == far]
        for sens in unique_sensitivities:
            group = far_df[far_df['Sensitivity'] == sens].sort_values('MIN')
            if group.empty:
                continue
            ax.plot(
                group['MIN'],
                group['MPCmean'],
                color=sens_color_map[sens],
                marker='o',
                linestyle='-',
                markersize=4,
                linewidth=1.4,
                alpha=0.85,
            )
        ax.set_title(f'FAR={far:g}')
        ax.grid(True, alpha=0.3)
        ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    for ax in axes[2:]:
        ax.set_xlabel('MIN', fontsize=11)
    for ax in axes[::2]:
        ax.set_ylabel('MPCmean', fontsize=11)

    fig.suptitle('MPCmean vs MIN by FAR\nBaseline eligibility only', fontsize=13)

    sens_handles = [mlines.Line2D([], [], color=sens_color_map[sens], marker='o', linestyle='-',
                                markersize=6, label=f'Sens={sens:g}') for sens in unique_sensitivities]

    fig.legend(handles=sens_handles, title='Sensitivity', loc='upper right',
               bbox_to_anchor=(0.99, 0.90), fontsize=8, frameon=True)

    fig.tight_layout(rect=[0, 0, 0.86, 0.93])
    plt.show()
    #plt.savefig('/mnt/user-data/outputs/mpc_plot.png', dpi=150, bbox_inches='tight')

def plot_MPCmean_vs_sensitivity_for_MIN(min_value=4):
    """Plot MPCmean vs sensitivity for one MIN value and baseline eligibility only."""
    df = pd.read_csv(CONST.OUTPUT_FILENAME)
    df = _baseline_only(df, 'BaselineTF')
    df = df[df['MIN'] == min_value].copy()
    if df.empty:
        raise ValueError(f"No baseline-eligibility rows found for MIN={min_value}.")

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for far, group in df.groupby('FAR'):
        group = group.sort_values('Sensitivity')
        ax.plot(
            group['Sensitivity'] * 100,
            group['MPCmean'],
            marker='o',
            linewidth=1.8,
            label=f'FAR={far:g}',
        )

    ax.set_xlabel('Sensitivity (%)')
    ax.set_ylabel('MPCmean')
    ax.set_title(f'MPCmean vs Sensitivity\nMIN={min_value}, baseline eligibility only')
    ax.grid(True, alpha=0.3)
    ax.legend(title='FAR', frameon=True)
    fig.tight_layout()
    plt.show()


def plot_MPCmean_vs_FAR_for_MIN(min_value=4):
    """Plot MPCmean vs FAR for one MIN value and baseline eligibility only."""
    df = pd.read_csv(CONST.OUTPUT_FILENAME)
    df = _baseline_only(df, 'BaselineTF')
    df = df[df['MIN'] == min_value].copy()
    if df.empty:
        raise ValueError(f"No baseline-eligibility rows found for MIN={min_value}.")

    unique_sensitivities = sorted(df['Sensitivity'].unique())
    sens_colors = plt.cm.viridis(np.linspace(0, 0.9, len(unique_sensitivities)))
    sens_color_map = {sens: sens_colors[i] for i, sens in enumerate(unique_sensitivities)}

    fig, ax = plt.subplots(figsize=(7, 4.5))
    for sens in unique_sensitivities:
        group = df[df['Sensitivity'] == sens].sort_values('FAR')
        ax.plot(
            group['FAR'],
            group['MPCmean'],
            color=sens_color_map[sens],
            marker='o',
            linewidth=1.5,
            label=f'Sens={sens:g}',
        )

    ax.set_xlabel('FAR')
    ax.set_ylabel('MPCmean')
    ax.set_title(f'MPCmean vs FAR\nMIN={min_value}, baseline eligibility only')
    ax.grid(True, alpha=0.3)
    ax.legend(title='Sensitivity', frameon=True, bbox_to_anchor=(1.02, 1), loc='upper left')
    fig.tight_layout()
    plt.show()

def draw_sens_and_far_vs_RTM(csv_path='rtm_test123_results.csv', correct_far=True):
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {'sensitivity', 'far', 'use_baseline', 'correct_far', 'frac_rtm', 'mpc'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])
    df = df[df['use_baseline'] == True]

    combos = [
        (True, correct_far),
    ]

    line_styles = {True: '-'}
    marker_styles = {True: 'o'}
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3'])
    combo_colors = {combo: color_cycle[i % len(color_cycle)] for i, combo in enumerate(combos)}

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    base_marker_size = 5
    emphasis_marker_size = base_marker_size * 2

    for combo in combos:
        use_baseline, correct_far = combo
        color = combo_colors[combo]
        linestyle = line_styles[use_baseline]
        marker = marker_styles[use_baseline]

        subset_a = df[
            (df['use_baseline'] == use_baseline) &
            (df['correct_far'] == correct_far) &
            (~((df['sensitivity'] == 1) & (df['far'] > 0)))
        ].sort_values('sensitivity')
        if not subset_a.empty:
            axes[0, 0].plot(
                subset_a['sensitivity'],
                subset_a['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            axes[1, 0].plot(
                subset_a['sensitivity'],
                subset_a['mpc'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            if use_baseline:
                sens_mask = subset_a['sensitivity'] == 1.0
                if sens_mask.any():
                    axes[0, 0].scatter(
                        subset_a.loc[sens_mask, 'sensitivity'],
                        subset_a.loc[sens_mask, 'frac_rtm'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
                    axes[1, 0].scatter(
                        subset_a.loc[sens_mask, 'sensitivity'],
                        subset_a.loc[sens_mask, 'mpc'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
        subset_b = df[
            (df['use_baseline'] == use_baseline) &
            (df['correct_far'] == correct_far) &
            (~((df['far'] == 0) & (df['sensitivity'] < 1)))
        ].sort_values('far')
        if not subset_b.empty:
            axes[0, 1].plot(
                subset_b['far'],
                subset_b['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            axes[1, 1].plot(
                subset_b['far'],
                subset_b['mpc'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            if use_baseline:
                far_mask = subset_b['far'] == 0.0
                if far_mask.any():
                    axes[0, 1].scatter(
                        subset_b.loc[far_mask, 'far'],
                        subset_b.loc[far_mask, 'frac_rtm'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
                    axes[1, 1].scatter(
                        subset_b.loc[far_mask, 'far'],
                        subset_b.loc[far_mask, 'mpc'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 70.0)
    axes[0, 0].text(0.5, 1.04, 'A', transform=axes[0, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[0, 1].text(0.5, 1.04, 'B', transform=axes[0, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 0].text(0.5, 1.04, 'C', transform=axes[1, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 1].text(0.5, 1.04, 'D', transform=axes[1, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')

    for ax in axes.flat:
        ax.grid(True, linestyle='--', alpha=0.3)
    axes[0, 0].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 1].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 0].set_xlim(0.0, 1.1)
    axes[1, 0].set_xlim(0.0, 1.1)
    sens_ticks = np.linspace(0.0, 1.0, 6)
    axes[0, 0].set_xticks(sens_ticks)
    axes[1, 0].set_xticks(sens_ticks)
    axes[0, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')
    axes[1, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')

    fig.subplots_adjust(bottom=0.08, wspace=0.15, hspace=0.35)
    fig.tight_layout()
    #plt.show()

def draw_sens_and_far_vs_RTM_v2(csv_path='rtm_test123_results.csv', correct_far=True):
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {'sensitivity', 'far', 'use_baseline', 'correct_far', 'frac_rtm', 'mpc'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])
    df = df[df['use_baseline'] == True]

    correct_far_value = _coerce_bool(pd.Series([correct_far])).iloc[0]
    if pd.isna(correct_far_value):
        raise ValueError(f"Could not interpret correct_far={correct_far!r} as boolean.")

    df_correct_far = df[df['correct_far'] == correct_far_value]
    df_correct_far_true = df[df['correct_far'] == True]

    combos = [
        (True, True),
    ]

    line_styles = {True: '-'}
    marker_styles = {True: 'o'}
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3'])
    combo_colors = {combo: color_cycle[i % len(color_cycle)] for i, combo in enumerate(combos)}

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    base_marker_size = 5
    emphasis_marker_size = base_marker_size * 2

    for combo in combos:
        use_baseline, _ = combo
        color = combo_colors[combo]
        linestyle = line_styles[use_baseline]
        marker = marker_styles[use_baseline]

        subset_a = df_correct_far_true[
            (df_correct_far_true['use_baseline'] == use_baseline) &
            (~((df_correct_far_true['sensitivity'] == 1) & (df_correct_far_true['far'] > 0)))
        ].sort_values('sensitivity')
        if not subset_a.empty:
            axes[0, 0].plot(
                subset_a['sensitivity'],
                subset_a['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            axes[1, 0].plot(
                subset_a['sensitivity'],
                subset_a['mpc'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            if use_baseline:
                sens_mask = subset_a['sensitivity'] == 1.0
                if sens_mask.any():
                    axes[0, 0].scatter(
                        subset_a.loc[sens_mask, 'sensitivity'],
                        subset_a.loc[sens_mask, 'frac_rtm'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
                    axes[1, 0].scatter(
                        subset_a.loc[sens_mask, 'sensitivity'],
                        subset_a.loc[sens_mask, 'mpc'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
        subset_b = df_correct_far[
            (df_correct_far['use_baseline'] == use_baseline) &
            (~((df_correct_far['far'] == 0) & (df_correct_far['sensitivity'] < 1)))
        ].sort_values('far')
        if not subset_b.empty:
            axes[0, 1].plot(
                subset_b['far'],
                subset_b['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            axes[1, 1].plot(
                subset_b['far'],
                subset_b['mpc'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            if use_baseline:
                far_mask = subset_b['far'] == 0.0
                if far_mask.any():
                    axes[0, 1].scatter(
                        subset_b.loc[far_mask, 'far'],
                        subset_b.loc[far_mask, 'frac_rtm'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
                    axes[1, 1].scatter(
                        subset_b.loc[far_mask, 'far'],
                        subset_b.loc[far_mask, 'mpc'],
                        s=emphasis_marker_size ** 2,
                        marker=marker,
                        color=color,
                        edgecolor=color,
                        zorder=3,
                    )
    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 70.0)
    axes[0, 0].text(0.5, 1.04, 'A', transform=axes[0, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[0, 1].text(0.5, 1.04, 'B', transform=axes[0, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 0].text(0.5, 1.04, 'C', transform=axes[1, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 1].text(0.5, 1.04, 'D', transform=axes[1, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')

    for ax in axes.flat:
        ax.grid(True, linestyle='--', alpha=0.3)
    axes[0, 0].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 1].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 0].set_xlim(0.0, 1.1)
    axes[1, 0].set_xlim(0.0, 1.1)
    sens_ticks = np.linspace(0.0, 1.0, 6)
    axes[0, 0].set_xticks(sens_ticks)
    axes[1, 0].set_xticks(sens_ticks)
    axes[0, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')
    axes[1, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')

    fig.subplots_adjust(bottom=0.08, wspace=0.15, hspace=0.35)
    fig.tight_layout()
    #plt.show()

def draw_sens_and_far_vs_RTM_mpc_ci(csv_path='rtm_test123_mpc_ci_results.csv', correct_far=True):
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {
        'sensitivity', 'far', 'use_baseline', 'correct_far', 'frac_rtm',
        'mpc_median', 'mpc_ci_lo', 'mpc_ci_hi'
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])
    df = df[df['use_baseline'] == True]

    correct_far_value = _coerce_bool(pd.Series([correct_far])).iloc[0]
    if pd.isna(correct_far_value):
        raise ValueError(f"Could not interpret correct_far={correct_far!r} as boolean.")

    df_correct_far = df[df['correct_far'] == correct_far_value]
    df_correct_far_true = df[df['correct_far'] == True]

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    color = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0'])[0]
    base_marker_size = 1

    def _mpc_yerr(subset):
        low = subset['mpc_median'] - subset['mpc_ci_lo']
        high = subset['mpc_ci_hi'] - subset['mpc_median']
        return np.vstack([low.clip(lower=0.0), high.clip(lower=0.0)])

    def _add_reference_line(ax, subset, x_col, x_value):
        ref = subset[np.isclose(subset[x_col], x_value)]
        if not ref.empty:
            ax.axhline(
                ref.iloc[0]['mpc_median'],
                color='0.5',
                linestyle='--',
                linewidth=1.0,
                alpha=0.8,
                zorder=0,
            )

    subset_a = df_correct_far_true[
        ~((df_correct_far_true['sensitivity'] == 1) & (df_correct_far_true['far'] > 0))
    ].sort_values('sensitivity')
    if not subset_a.empty:
        axes[0, 0].plot(
            subset_a['sensitivity'],
            subset_a['frac_rtm'],
            linestyle='-',
            marker='o',
            color=color,
            markersize=base_marker_size,
            linewidth=1.8,
        )
        axes[1, 0].errorbar(
            subset_a['sensitivity'],
            subset_a['mpc_median'],
            yerr=_mpc_yerr(subset_a),
            linestyle='-',
            marker='o',
            color=color,
            ecolor=color,
            capsize=3,
            markersize=base_marker_size,
            linewidth=1.8,
        )
        _add_reference_line(axes[1, 0], subset_a, 'sensitivity', 1.0)

    subset_b = df_correct_far[
        ~((df_correct_far['far'] == 0) & (df_correct_far['sensitivity'] < 1))
    ].sort_values('far')
    if not subset_b.empty:
        axes[0, 1].plot(
            subset_b['far'],
            subset_b['frac_rtm'],
            linestyle='-',
            marker='o',
            color=color,
            markersize=base_marker_size,
            linewidth=1.8,
        )
        _add_reference_line(axes[1, 1], subset_b, 'far', 0.0)
        axes[1, 1].errorbar(
            subset_b['far'],
            subset_b['mpc_median'],
            yerr=_mpc_yerr(subset_b),
            linestyle='-',
            marker='o',
            color=color,
            ecolor=color,
            capsize=3,
            markersize=base_marker_size,
            linewidth=1.8,
        )

    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Median placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 70.0)
    axes[0, 0].text(0.5, 1.04, 'A', transform=axes[0, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[0, 1].text(0.5, 1.04, 'B', transform=axes[0, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 0].text(0.5, 1.04, 'C', transform=axes[1, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 1].text(0.5, 1.04, 'D', transform=axes[1, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')

    for ax in axes.flat:
        ax.grid(True, linestyle='--', alpha=0.3)
    axes[0, 0].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 1].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 0].set_xlim(0.0, 1.1)
    axes[1, 0].set_xlim(0.0, 1.1)
    sens_ticks = np.linspace(0.0, 1.0, 6)
    axes[0, 0].set_xticks(sens_ticks)
    axes[1, 0].set_xticks(sens_ticks)
    axes[0, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')
    axes[1, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')

    fig.subplots_adjust(bottom=0.08, wspace=0.15, hspace=0.35)
    fig.tight_layout()
    return fig, axes

def draw_combined_sens_and_far_vs_RTM_mpc_ci(csv_path='rtm_test123_mpc_ci_results.csv'):
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {
        'sensitivity', 'far', 'use_baseline', 'correct_far', 'frac_rtm',
        'mpc_median', 'mpc_ci_lo', 'mpc_ci_hi'
    }
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])
    df = df[df['use_baseline'] == True]

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    style = {
        True: {'color': 'black', 'label': 'FAR corrected'},
        False: {'color': 'red', 'label': 'FAR uncorrected'},
    }
    base_marker_size = 1

    def _mpc_yerr(subset):
        low = subset['mpc_median'] - subset['mpc_ci_lo']
        high = subset['mpc_ci_hi'] - subset['mpc_median']
        return np.vstack([low.clip(lower=0.0), high.clip(lower=0.0)])

    def _add_reference_line(ax, subset, x_col, x_value):
        ref = subset[np.isclose(subset[x_col], x_value)]
        if not ref.empty:
            ax.axhline(
                ref.iloc[0]['mpc_median'],
                color='0.5',
                linestyle='--',
                linewidth=1.0,
                alpha=0.8,
                zorder=0,
            )

    for correct_far_value, attrs in style.items():
        subset = df[df['correct_far'] == correct_far_value]
        subset_a = subset[
            (subset['far'] == 0.0) &
            (~((subset['sensitivity'] == 1) & (subset['far'] > 0)))
        ].sort_values('sensitivity')
        if not subset_a.empty:
            axes[0, 0].plot(
                subset_a['sensitivity'],
                subset_a['frac_rtm'],
                linestyle='-',
                marker='o',
                color=attrs['color'],
                markersize=base_marker_size,
                linewidth=1.8,
                label=attrs['label'],
            )
            axes[1, 0].errorbar(
                subset_a['sensitivity'],
                subset_a['mpc_median'],
                yerr=_mpc_yerr(subset_a),
                linestyle='-',
                marker='o',
                color=attrs['color'],
                ecolor=attrs['color'],
                capsize=3,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            _add_reference_line(axes[1, 0], subset_a, 'sensitivity', 1.0)

        subset_b = subset[subset['sensitivity'] == 1.0].sort_values('far')
        if not subset_b.empty:
            axes[0, 1].plot(
                subset_b['far'],
                subset_b['frac_rtm'],
                linestyle='-',
                marker='o',
                color=attrs['color'],
                markersize=base_marker_size,
                linewidth=1.8,
                label=attrs['label'],
            )
            axes[1, 1].errorbar(
                subset_b['far'],
                subset_b['mpc_median'],
                yerr=_mpc_yerr(subset_b),
                linestyle='-',
                marker='o',
                color=attrs['color'],
                ecolor=attrs['color'],
                capsize=3,
                markersize=base_marker_size,
                linewidth=1.8,
            )
            _add_reference_line(axes[1, 1], subset_b, 'far', 0.0)

    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Median placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 70.0)
    axes[0, 0].text(0.5, 1.04, 'A', transform=axes[0, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[0, 1].text(0.5, 1.04, 'B', transform=axes[0, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 0].text(0.5, 1.04, 'C', transform=axes[1, 0].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')
    axes[1, 1].text(0.5, 1.04, 'D', transform=axes[1, 1].transAxes,
                 fontsize=14, fontweight='bold', va='bottom', ha='center')

    for ax in axes.flat:
        ax.grid(True, linestyle='--', alpha=0.3)
    axes[0, 0].legend(frameon=True, facecolor='white', framealpha=0.90)
    axes[0, 1].legend(frameon=True, facecolor='white', framealpha=0.90)
    axes[0, 0].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 1].yaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}%')
    axes[0, 0].set_xlim(0.0, 1.1)
    axes[1, 0].set_xlim(0.0, 1.1)
    sens_ticks = np.linspace(0.0, 1.0, 6)
    axes[0, 0].set_xticks(sens_ticks)
    axes[1, 0].set_xticks(sens_ticks)
    axes[0, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')
    axes[1, 0].xaxis.set_major_formatter(lambda x, pos: f'{x * 100:.0f}')

    fig.subplots_adjust(bottom=0.08, wspace=0.15, hspace=0.35)
    fig.tight_layout()
    return fig, axes

def plot_efficacy_result(fn='rtm_efficacy_results_100000.csv'):
    df = pd.read_csv(fn)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {'sensitivity', 'far', 'n_for_90m'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")
    if 'usebaselinetf' in df.columns:
        df = df[_coerce_bool(df['usebaselinetf']) == True]
    
    # Filter for 90% power
    df_90 = df[df['power'] == 0.90] if 'power' in df.columns else df
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Subplot A: Sensitivity vs N for MPC
    for far in sorted(df_90['far'].unique()):
        subset = df_90[df_90['far'] == far].sort_values('sensitivity')
        axes[0].plot(subset['sensitivity'], subset['n_for_90m'], marker='o', label=f'FAR={far}')
    
    axes[0].set_xlabel('Sensitivity')
    axes[0].set_ylabel('N')
    axes[0].set_title('MPC: Sensitivity vs N (90% Power, baseline eligibility)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].text(-0.1, 1.05, 'A', transform=axes[0].transAxes,
                    fontsize=14, fontweight='bold', va='bottom', ha='right')
    
    # Subplot B: FAR vs N for NPC
    for sens in sorted(df_90['sensitivity'].unique()):
        subset = df_90[df_90['sensitivity'] == sens].sort_values('far')
        axes[1].plot(subset['far'], subset['n_for_90m'], marker='s', label=f'Sens={sens}')
    
    axes[1].set_xlabel('FAR')
    axes[1].set_ylabel('N')
    axes[1].set_title('MPC: FAR vs N (90% Power, baseline eligibility)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].text(-0.1, 1.05, 'B', transform=axes[1].transAxes,
                    fontsize=14, fontweight='bold', va='bottom', ha='right')
    
    fig.tight_layout()
    plt.show()
