import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import rtm_constants as CONST
import seaborn as sns




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
        plot_heatmap(mSF_values, PC_values, eligibility_min)
        plot_heatmap(baseline_values, PC_values, eligibility_min)

def plot_heatmap(mSF_values, PC_values, eligibility_min):
    bins=CONST.HEATMAP_BINS
    plt.figure(figsize=(3, 2))
    
    # Create 2D histogram as heatmap with fixed range
    heatmap, _, _ = np.histogram2d(
        mSF_values, PC_values,
        bins=bins,
        range=[[0, 15], [-100, 100]], density=True
    )

    # Plot using imshow
    plt.imshow(heatmap.T, origin='lower', aspect='auto',
            extent=[0, 15, -100, 100],
            cmap='viridis')

    plt.axvline(x=eligibility_min, color='red', linewidth=3)
    plt.colorbar(label='Count')
    plt.xlabel('Monthly Seizure Frequency (mSF)')
    plt.ylabel('Percent Change (PC)')
    plt.title('Cohort Heatmap: mSF vs Percent Change')
    
    plt.tight_layout()
    plt.show()
    
    return plt.gcf()

def make_FAR_plot():

    # -*- coding: utf-8 -*-
    """
    MPC vs MIN by FAR using output.csv
    - Sensitivity == 1.00 only
    - BaselineTF=TRUE: solid line + closed markers
    - BaselineTF=FALSE: dashed line + open markers
    - Legend outside (right), two columns: left TRUE, right FALSE
    - FAR=0.0 lines highlighted (thicker)
    - Title: 'Effect of changing FAR'
    """


    # --- Load data from CSV ---
    csv_path = CONST.OUTPUT_FILENAME
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip() for c in df.columns}, inplace=True)

    # --- Validate required columns ---
    required_cols = {'Sensitivity', 'FAR', 'BaselineTF', 'MIN', 'MPC'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")  # [1](https://bilh-my.sharepoint.com/personal/dgoldenh_bidmc_harvard_edu/Documents/Microsoft%20Copilot%20Chat%20Files/output.csv)

    # --- Filter: Sensitivity == 1.00 ---
    sens_mask = pd.to_numeric(df['Sensitivity'], errors='coerce').eq(1.0)

    # --- Baseline masks (support boolean and typical string encodings) ---
    base_series = df['BaselineTF']
    base_str = base_series.astype(str).str.strip().str.upper()
    true_mask  = sens_mask & (base_str.eq('TRUE')  | base_series.eq(True))
    false_mask = sens_mask & (base_str.eq('FALSE') | base_series.eq(False) |
                            base_str.isin({'0','F','NO'}) | base_series.eq(0))

    true_df  = df.loc[true_mask,  ['FAR','MIN','MPC']].copy()
    false_df = df.loc[false_mask, ['FAR','MIN','MPC']].copy()

    # --- Aggregate to one MPC per (FAR, MIN) to avoid duplicate points ---
    true_agg  = (true_df.groupby(['FAR','MIN'], as_index=False)
                .agg(MPC=('MPC','mean')).sort_values(['FAR','MIN']))
    false_agg = (false_df.groupby(['FAR','MIN'], as_index=False)
                .agg(MPC=('MPC','mean')).sort_values(['FAR','MIN']))

    # --- Color map per FAR (same color for TRUE/FALSE of same FAR) ---
    all_fars = sorted(set(true_agg['FAR']).union(set(false_agg['FAR'])))
    cmap = plt.get_cmap('tab10')
    color_map = {far: cmap(i % 10) for i, far in enumerate(all_fars)}

    # --- Line width rules (highlight FAR=0.0) ---
    base_lw = 2.0
    highlight_lw = 4.0

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10.5, 6.2))

    # TRUE: solid + closed markers
    for far_val, grp in true_agg.groupby('FAR'):
        color = color_map[far_val]
        lw = highlight_lw if float(far_val) == 0.0 else base_lw
        ax.plot(grp['MIN'], grp['MPC'],
                linestyle='-', linewidth=lw, color=color,
                marker='o', markersize=5,
                markerfacecolor=color, markeredgecolor=color)

    # FALSE: dashed + open markers
    for far_val, grp in false_agg.groupby('FAR'):
        color = color_map[far_val]
        lw = highlight_lw if float(far_val) == 0.0 else base_lw
        ax.plot(grp['MIN'], grp['MPC'],
                linestyle='--', linewidth=lw, color=color,
                marker='o', markersize=5,
                markerfacecolor='none', markeredgecolor=color)

    # --- Axes labels and title ---
    ax.set_title("Effect of changing FAR")  # <-- updated title
    ax.set_xlabel('MIN')
    ax.set_ylabel('MPC')
    ax.grid(True, linestyle='--', alpha=0.35)

    # --- Legend proxies (two columns: left TRUE, right FALSE) ---
    true_proxy = [Line2D([0],[0], color=color_map[far], linestyle='-',
                        marker='o', markersize=6,
                        markerfacecolor=color_map[far], markeredgecolor=color_map[far],
                        linewidth=(highlight_lw if float(far)==0.0 else base_lw))
                for far in all_fars if far in set(true_agg['FAR'])]

    false_proxy = [Line2D([0],[0], color=color_map[far], linestyle='--',
                        marker='o', markersize=6,
                        markerfacecolor='none', markeredgecolor=color_map[far],
                        linewidth=(highlight_lw if float(far)==0.0 else base_lw))
                for far in all_fars if far in set(false_agg['FAR'])]

    true_labels  = [f'FAR={far}' for far in all_fars if far in set(true_agg['FAR'])]
    false_labels = [f'FAR={far}' for far in all_fars if far in set(false_agg['FAR'])]

    handles_combined = true_proxy + false_proxy
    labels_combined  = true_labels + false_labels

    # Set x-axis to show only integer tick values
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))

    # Place legend OUTSIDE to the right
    leg = ax.legend(handles_combined, labels_combined, ncol=2,
                    title="Legend — left: BaselineTF=TRUE | right: BaselineTF=FALSE",
                    frameon=True, facecolor='white', framealpha=0.90,
                    columnspacing=1.0, handlelength=3.0,
                    loc='center left', bbox_to_anchor=(1.02, 0.5))

    # Make room on the right for the legend
    plt.subplots_adjust(right=0.75)

    fig.tight_layout()
    plt.show()

def plot_all_combos_of_meanMPC_vs_MIN():
    # -*- coding: utf-8 -*-
    """
    MPC vs MIN by FAR using output.csv
    - Sensitivity == 1.00 only
    - BaselineTF=TRUE: solid line + closed markers
    - BaselineTF=FALSE: dashed line + open markers
    - Legend outside (right), two columns: left TRUE, right FALSE
    - FAR=0.0 lines highlighted (thicker)"""


    '''"""
    Build a plot of x=MIN and y=MPCmean. each combination of sensitivity and FAR will be plotted with a different color and marker combination for the line - but if the BaselineTF=TRUE then the line should be solid, and if the BaselineTF=FALSE then the line should be dashed.

    10:26 PM

    View the CSV file to understand its structure
    Show working file
    Now I understand the data structure. Let me create a Python script to generate the plot:


    1 step

    Run the plotting script



    Check the output plot
    Show working file
    The plot is created but the legend is quite crowded due to the many combinations. Let me create an improved version with a cleaner legend that shows the color/marker coding for Sensitivity+FAR separately from the line style for Baseline:


    Hide steps

    Create improved plot with cleaner legend organization
    plot_data_v2.py

    Run the improved plotting script



    Check the improved output plot
    Show working file

    Presented file
    Here's your plot! The visualization shows:

    Colors represent different Sensitivity values (0.25 to 1.0, using a viridis colormap)
    Markers represent different FAR values (circles for 0.0, squares for ~0.033, triangles for ~0.143, diamonds for 1.0)
    Solid lines indicate BaselineTF=True
    Dashed lines indicate BaselineTF=False
    You can see that the Baseline=True conditions (solid lines) generally have positive MPCmean values that increase with MIN, while Baseline=False conditions (dashed lines) cluster around zero or slightly negative values.

    Mpc plot
    Image · PNG 





    Claude is AI and can make mistakes. Please double-check responses.
    Plot data v2 · PY
    Copy
    '''

    # Read the data
    df = pd.read_csv(CONST.OUTPUT_FILENAME)

    # Get unique combinations of Sensitivity and FAR
    unique_combos = df[['Sensitivity', 'FAR']].drop_duplicates().sort_values(['Sensitivity', 'FAR'])

    # Define colors for different sensitivities and markers for different FARs
    unique_sensitivities = sorted(df['Sensitivity'].unique())
    unique_fars = sorted(df['FAR'].unique())

    # Color map for sensitivities
    sens_colors = plt.cm.viridis(np.linspace(0, 0.9, len(unique_sensitivities)))
    sens_color_map = {sens: sens_colors[i] for i, sens in enumerate(unique_sensitivities)}

    # Markers for FAR values
    far_markers = ['o', 's', '^', 'D', 'v', '<', '>', 'p', 'h', '*']
    far_marker_map = {far: far_markers[i % len(far_markers)] for i, far in enumerate(unique_fars)}

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 9))

    # Plot each group
    for (sens, far, baseline), group in df.groupby(['Sensitivity', 'FAR', 'BaselineTF']):
        color = sens_color_map[sens]
        marker = far_marker_map[far]
        linestyle = '-' if baseline else '--'
        
        # Sort by MIN for proper line plotting
        group_sorted = group.sort_values('MIN')
        
        ax.plot(group_sorted['MIN'], group_sorted['MPCmean'], 
                color=color, 
                marker=marker,
                linestyle=linestyle,
                markersize=6,
                linewidth=1.5,
                alpha=0.8)

    # Customize plot
    ax.set_xlabel('MIN', fontsize=12)
    ax.set_ylabel('MPCmean', fontsize=12)
    ax.set_title('MPCmean vs MIN\n(Color = Sensitivity, Marker = FAR, Solid = Baseline True, Dashed = Baseline False)', fontsize=13)
    ax.grid(True, alpha=0.3)

    # Create custom legend
    # First part: Sensitivity colors
    sens_handles = [mlines.Line2D([], [], color=sens_color_map[sens], marker='o', linestyle='-', 
                                markersize=6, label=f'Sens={sens}') for sens in unique_sensitivities]

    # Second part: FAR markers  
    far_handles = [mlines.Line2D([], [], color='gray', marker=far_marker_map[far], linestyle='none',
                                markersize=8, label=f'FAR={far}') for far in unique_fars]

    # Third part: Baseline line styles
    baseline_handles = [
        mlines.Line2D([], [], color='black', linestyle='-', linewidth=2, label='Baseline=True'),
        mlines.Line2D([], [], color='black', linestyle='--', linewidth=2, label='Baseline=False')
    ]

    # Combine legends
    legend1 = ax.legend(handles=sens_handles, title='Sensitivity', loc='upper left', 
                        bbox_to_anchor=(1.01, 1), fontsize=9)
    ax.add_artist(legend1)

    legend2 = ax.legend(handles=far_handles, title='FAR', loc='upper left',
                        bbox_to_anchor=(1.01, 0.55), fontsize=9)
    ax.add_artist(legend2)

    legend3 = ax.legend(handles=baseline_handles, title='Line Style', loc='upper left',
                        bbox_to_anchor=(1.01, 0.25), fontsize=9)

    plt.tight_layout()
    plt.show()
    #plt.savefig('/mnt/user-data/outputs/mpc_plot.png', dpi=150, bbox_inches='tight')

def draw_sens_and_far_vs_RTM(csv_path='rtm_test123_results.csv', correct_far=True):
    df = pd.read_csv(csv_path)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)

    required_cols = {'sensitivity', 'far', 'use_baseline', 'correct_far', 'frac_rtm', 'mpc'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")

    def _coerce_bool(series):
        if series.dtype == bool:
            return series
        values = series.astype(str).str.strip().str.upper()
        true_set = {'TRUE', 'T', '1', 'YES', 'Y'}
        false_set = {'FALSE', 'F', '0', 'NO', 'N'}
        return values.map(lambda v: True if v in true_set else False if v in false_set else np.nan)

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])

    combos = [
        (True, correct_far),
        (False, correct_far),
    ]

    line_styles = {True: '-', False: '--'}
    marker_styles = {True: 'o', False: '^'}
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3'])
    combo_colors = {combo: color_cycle[i % len(color_cycle)] for i, combo in enumerate(combos)}

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    base_marker_size = 5
    emphasis_marker_size = base_marker_size * 2

    handles = []
    labels = []
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
            line = axes[0, 0].plot(
                subset_a['sensitivity'],
                subset_a['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )[0]
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
            handles.append(line)
            label = 'Eligibility from baseline' if use_baseline else 'Eligibility before baseline'
            labels.append(label)

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
            label = 'Eligibility from baseline' if use_baseline else 'Eligibility before baseline'
            if label not in labels:
                handles.append(axes[0, 1].lines[-1])
                labels.append(label)

    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 50.0)
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

    if handles:
        axes[0, 1].legend(handles, labels, loc='upper right',
                       frameon=True, framealpha=0.95)
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

    def _coerce_bool(series):
        if series.dtype == bool:
            return series
        values = series.astype(str).str.strip().str.upper()
        true_set = {'TRUE', 'T', '1', 'YES', 'Y'}
        false_set = {'FALSE', 'F', '0', 'NO', 'N'}
        return values.map(lambda v: True if v in true_set else False if v in false_set else np.nan)

    df['use_baseline'] = _coerce_bool(df['use_baseline'])
    df['correct_far'] = _coerce_bool(df['correct_far'])
    df = df.dropna(subset=['use_baseline', 'correct_far'])

    correct_far_value = _coerce_bool(pd.Series([correct_far])).iloc[0]
    if pd.isna(correct_far_value):
        raise ValueError(f"Could not interpret correct_far={correct_far!r} as boolean.")

    df_correct_far = df[df['correct_far'] == correct_far_value]
    df_correct_far_true = df[df['correct_far'] == True]

    combos = [
        (True, True),
        (False, True),
    ]

    line_styles = {True: '-', False: '--'}
    marker_styles = {True: 'o', False: '^'}
    color_cycle = plt.rcParams['axes.prop_cycle'].by_key().get('color', ['C0', 'C1', 'C2', 'C3'])
    combo_colors = {combo: color_cycle[i % len(color_cycle)] for i, combo in enumerate(combos)}

    fig, axes = plt.subplots(2, 2, figsize=(7, 5), sharex='col', sharey='row')
    base_marker_size = 5
    emphasis_marker_size = base_marker_size * 2

    handles = []
    labels = []
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
            line = axes[0, 0].plot(
                subset_a['sensitivity'],
                subset_a['frac_rtm'],
                linestyle=linestyle,
                marker=marker,
                color=color,
                markersize=base_marker_size,
                linewidth=1.8,
            )[0]
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
            handles.append(line)
            label = 'Eligibility from baseline' if use_baseline else 'Eligibility before baseline'
            labels.append(label)

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
            label = 'Eligibility from baseline' if use_baseline else 'Eligibility before baseline'
            if label not in labels:
                handles.append(axes[0, 1].lines[-1])
                labels.append(label)

    axes[0, 0].set_xlabel('Sensitivity (%)')
    axes[0, 1].set_xlabel('False alarm rate: alarms/day')
    axes[1, 0].set_xlabel('Sensitivity (%)')
    axes[1, 1].set_xlabel('False alarm rate: alarms/day')
    axes[0, 0].set_ylabel('RTM (%)')
    axes[1, 0].set_ylabel('Placebo MPC (%)')
    axes[0, 0].set_ylim(0.0, 1.0)
    axes[1, 0].set_ylim(-10.0, 50.0)
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

    if handles:
        axes[0, 1].legend(handles, labels, loc='upper right',
                       frameon=True, framealpha=0.95)
    fig.subplots_adjust(bottom=0.08, wspace=0.15, hspace=0.35)
    fig.tight_layout()
    #plt.show()

def plot_efficacy_result(fn='efficacy_results_100000.csv'):
    df = pd.read_csv(fn)
    df.rename(columns={c: str(c).strip().lower() for c in df.columns}, inplace=True)
    
    Sensitivity,FAR,CorrectFAR,useBaselineTF,drugEffect,N_for_90M

    required_cols = {'Sensitivity', 'FAR', 'N_for_90M'}
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}. Found columns: {list(df.columns)}")
    
    # Filter for 90% power
    df_90 = df[df['power'] == 0.90] if 'power' in df.columns else df
    
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    
    # Subplot A: Sensitivity vs N for MPC
    for far in sorted(df_90['far'].unique()):
        subset = df_90[df_90['far'] == far].sort_values('sensitivity')
        axes[0].plot(subset['sensitivity'], subset['n'], marker='o', label=f'FAR={far}')
    
    axes[0].set_xlabel('Sensitivity')
    axes[0].set_ylabel('N')
    axes[0].set_title('MPC: Sensitivity vs N (90% Power)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    axes[0].text(-0.1, 1.05, 'A', transform=axes[0].transAxes,
                    fontsize=14, fontweight='bold', va='bottom', ha='right')
    
    # Subplot B: FAR vs N for NPC
    for sens in sorted(df_90['sensitivity'].unique()):
        subset = df_90[df_90['sensitivity'] == sens].sort_values('far')
        axes[1].plot(subset['far'], subset['n'], marker='s', label=f'Sens={sens}')
    
    axes[1].set_xlabel('FAR')
    axes[1].set_ylabel('N')
    axes[1].set_title('NPC: FAR vs N (90% Power)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    axes[1].text(-0.1, 1.05, 'B', transform=axes[1].transAxes,
                    fontsize=14, fontweight='bold', va='bottom', ha='right')
    
    fig.tight_layout()
    plt.show()