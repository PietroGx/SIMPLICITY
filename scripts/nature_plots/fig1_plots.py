import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

# =============================================================================
# One colour language for the whole figure (Okabe-Ito, colourblind-safe).
# Patient data gets its own hue AND a dashed line, so data-vs-model survives
# greyscale printing; the old dark blue sat too close to Control's near-black.
# =============================================================================
SCEN_COLORS = {
    "Control":   "#000000",   # black
    "SOT":       "#0072B2",   # blue
    "HIV":       "#D55E00",   # vermillion
    "Edge Case": "#CC79A7",   # reddish purple
}
DATA_COLOR = "#009E73"        # bluish green -- used by no scenario
DATA_STYLE = (0, (4, 2))

# Panel C keys clinical categories to the scenario they stand in for, so the
# colours mean the same thing there as everywhere else.
CATEGORY_TO_SCENARIO = {
    "SOT": "SOT",
    "HIV": "HIV",
}
LONG_SHEDDER_CUTOFF_DAYS = 300
UNMODELLED_COLOR = "#B0B0B0"

# Fit annotations are for us, not for the journal. Flip to False (or pass
# --no-fit-stats) for the submission version.
SHOW_FIT_STATS = True

LEGEND_ENTRIES = [
    ("Control (standard)", SCEN_COLORS["Control"]),
    ("SOT (solid organ transplant)", SCEN_COLORS["SOT"]),
    ("HIV (HIV/AIDS)", SCEN_COLORS["HIV"]),
    ("Edge case (>%d d shedding)" % LONG_SHEDDER_CUTOFF_DAYS, SCEN_COLORS["Edge Case"]),
]


def add_figure_legend(fig, has_edge_case=True, has_patient_data=True):
    """Single legend for the whole figure, below the panels."""
    from matplotlib.lines import Line2D
    handles = []
    for label, colour in LEGEND_ENTRIES:
        if "Edge case" in label and not has_edge_case:
            continue
        handles.append(Line2D([0], [0], color=colour, lw=2.4, label=label))
    if has_patient_data:
        handles.append(Line2D([0], [0], color=DATA_COLOR, lw=2.4,
                              linestyle=DATA_STYLE, label="Patient data"))
    fig.legend(handles=handles, loc='lower center', ncol=len(handles),
               frameon=False, fontsize=6.4, bbox_to_anchor=(0.5, -0.015),
               handlelength=2.2, columnspacing=1.5, handletextpad=0.5)


def format_clean_axis(ax, remove_ticks=True):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

def plot_fig1_intra_host(ax, df, x_max=800):
    if df.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, "Panel A\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return
    palette = SCEN_COLORS
    sns.lineplot(data=df, x='time', y='p_infectious', hue='cohort', palette=palette, ax=ax, linewidth=1.5, alpha=0.8, legend=False)
    ax.set_xlim(0, x_max)
    ax.set_ylim(bottom=0, top=1.05)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:g}"))
    ax.set_xlabel("Days since infection", fontsize=7)
    ax.set_ylabel("Probability of being infectious", fontsize=7)
    format_clean_axis(ax, remove_ticks=False)

def plot_fig1_violins(ax, df):
    if df.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, "Panel B\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return
    palette = SCEN_COLORS
    order = ["Control", "SOT", "HIV", "Edge Case"]
    existing_order = [c for c in order if c in df['cohort'].unique()]
    sns.violinplot(data=df, x='cohort', y='duration', palette=palette, hue='cohort', order=existing_order, ax=ax, inner="quartile", cut=0, linewidth=0.8, density_norm='width', legend=False)
    ax.set_ylim(0, 600)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:g}"))
    ax.set_xlabel("", fontsize=7)
    ax.set_ylabel("Realized Duration (days)", fontsize=7)
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha('right')
    format_clean_axis(ax, remove_ticks=False)

def plot_infection_duration(ax, df):
    if df.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, "Panel C\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return
    acronym_map = {
        'HIV/AIDS': 'HIV', 'SOT (Solid Organ Transplant)': 'SOT', 'B-Cell neoplasm': 'B-CN',
        'Rheumatological/Autoimmune': 'AuIm', 'PID (Primary Immunodeficiency)': 'PID',
        'aHSCT': 'aHSCT', 'Other hematological diseases': 'OHem',
        'Immunocompromised': 'ImmC', 'Other': 'Other'
    }
    df['clinical_category'] = df['clinical_category'].map(lambda x: acronym_map.get(x, str(x)[:5]))
    avg_durations = df.groupby('clinical_category')['duration'].mean().sort_values()
    sorted_categories = avg_durations.index.tolist()
    
    # Same colours as the modelled scenarios: a clinical category is drawn in
    # the colour of the scenario that stands in for it, and unmodelled
    # categories stay grey. The shared legend names both.
    box_palette = {
        cat: SCEN_COLORS.get(CATEGORY_TO_SCENARIO.get(cat), UNMODELLED_COLOR)
        for cat in sorted_categories
    }
        
    sns.boxplot(data=df, x='clinical_category', y='duration', order=sorted_categories, ax=ax, 
                palette=box_palette, hue='clinical_category', showfliers=False, width=0.5, boxprops=dict(alpha=0.4), legend=False)
    
    def get_pt_color(row):
        # Beyond the cutoff a patient is an edge case whatever their category.
        if row['duration'] > LONG_SHEDDER_CUTOFF_DAYS:
            return SCEN_COLORS["Edge Case"]
        return box_palette[row['clinical_category']]
    df['pt_color'] = df.apply(get_pt_color, axis=1)

    for i, cat in enumerate(sorted_categories):
        cat_data = df[df['clinical_category'] == cat]
        x_jitter = np.random.normal(i, 0.05, size=len(cat_data))
        ax.scatter(x_jitter, cat_data['duration'], c=cat_data['pt_color'], s=12, alpha=0.8, zorder=5)

    ax.set_ylim(0, 600)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f"{y:g}"))
    ax.set_xlabel("", fontsize=7)
    ax.set_ylabel("Literature Duration (days)", fontsize=7)
    ax.set_xticks(range(len(sorted_categories)))
    ax.set_xticklabels(sorted_categories)
    for label in ax.get_xticklabels():
        label.set_rotation(45)
        label.set_ha('right')
    format_clean_axis(ax, remove_ticks=False)

def plot_tempest_regression(ax, df, title, x_col='sampling_date', y_col='hamming_per_overlap', is_date=True, scatter_color='#404040', line_color='black', force_zero_intercept=False, force_intercept_at_min_x=False, x_label="Days since infection"):
    if df.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, f"{title}\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return

    # SIMULATED DATA BRANCH (Panels F & G)
    if 'Sequencing_time' in df.columns and 'Distance_from_root' in df.columns:
        import simplicity.tuning.evolutionary_rate as er
        palette = {"Control": "#333333", "SOT": "#56B4E9", "HIV": "#D55E00", "Edge Case": "#CC79A7"}
        cohorts = df['cohort'].unique() if 'cohort' in df.columns else ['Control']
        
        stats_texts = []
        for cohort in cohorts:
            cdf = df[df['cohort'] == cohort] if 'cohort' in df.columns else df
            if cdf.empty: continue
            
            try:
                fit_cdf = cdf.sample(frac=0.25, random_state=42) if len(cdf) > 500 else cdf
                fitted_model = er.tempest_regression(fit_cdf)
                slope = fitted_model.coef_[0]
                x_vals = fit_cdf['Sequencing_time'].values.reshape(-1, 1)
                y_vals = fit_cdf['Distance_from_root'].values
                r_squared = fitted_model.score(x_vals, y_vals)
            except Exception: continue
            
            c_color = palette.get(cohort, scatter_color)
            
            # Adjusted sampling logic so dots remain visible
            plot_cdf = cdf.sample(frac=0.1, random_state=42) if len(cdf) > 500 else cdf
            
            x_plot_days = plot_cdf['Sequencing_time'].values * 365.25
            y_plot = plot_cdf['Distance_from_root'].values
            ax.scatter(x_plot_days, y_plot, alpha=0.3, s=8, color=c_color, edgecolors='none', zorder=3)
            
            x_line_years = np.array([0, cdf['Sequencing_time'].max()])
            y_line = fitted_model.predict(x_line_years.reshape(-1, 1))
            ax.plot(x_line_years * 365.25, y_line, color=c_color, linewidth=1.5, zorder=4)
            
            stats_texts.append(f"{cohort}: $R^2={r_squared:.2f}$, Rate={slope:.5f} s/s/y")
            
        ax.text(0.05, 0.95, "\n".join(stats_texts), transform=ax.transAxes, fontsize=6, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_ylabel("Hamming D.", fontsize=7)
        ax.set_xlim(left=0, right=400)
        ax.set_ylim(bottom=0, top=0.002)
        ax.text(0.5, 1.05, title, transform=ax.transAxes, ha='center', fontsize=7, fontweight='bold')
        format_clean_axis(ax, remove_ticks=False)
        return

    # REAL DATA BRANCH (Panels D & E)
    df = df.dropna(subset=[x_col, y_col]).copy()
    df = df.sort_values(x_col)
    y = df[y_col].values
    if is_date:
        df[x_col] = pd.to_datetime(df[x_col])
        x_dates = df[x_col]
        x_plot = mdates.date2num(x_dates.dt.to_pydatetime())
    else:
        x_plot = df[x_col].astype(float).values
    if force_zero_intercept:
        slope = np.sum(x_plot * y) / np.sum(x_plot**2)
        r_squared = 1 - (np.sum((y - slope * x_plot)**2) / np.sum((y - np.mean(y))**2))
        ax.plot([0, np.max(x_plot)], [0, slope * np.max(x_plot)], color=line_color, linewidth=1.2, zorder=1)
    elif force_intercept_at_min_x:
        x_min = np.min(x_plot)
        x_shifted = x_plot - x_min
        slope = np.sum(x_shifted * y) / np.sum(x_shifted**2)
        r_squared = 1 - (np.sum((y - slope * x_shifted)**2) / np.sum((y - np.mean(y))**2))
        ax.plot([x_min, np.max(x_plot)], [0, slope * (np.max(x_plot) - x_min)], color=line_color, linewidth=1.2, zorder=1)
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_plot, y)
        r_squared = r_value**2
        sns.regplot(x=x_plot, y=y, ax=ax, scatter=False, color=line_color, ci=95, line_kws={'linewidth': 1.2})
    ax.scatter(x_plot, y, alpha=0.5, s=8, color=scatter_color, edgecolors='none', zorder=3)
    rate_per_year = slope * 365.25 
    stats_text = f"$R^2 = {r_squared:.2f}$\nRate $= {rate_per_year:.5f}$ s/s/y"
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, fontsize=6, va='top', ha='left', bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))
    if is_date:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_ha('right')
        ax.set_xlim(left=np.min(x_plot), right=np.max(x_plot))
    else:
        ax.set_xlabel(x_label, fontsize=7)
        ax.set_xlim(left=0, right=np.max(x_plot) * 1.05)
    ax.set_ylim(bottom=0)
    ax.set_ylabel("Hamming D.", fontsize=7)
    ax.text(0.5, 1.05, title, transform=ax.transAxes, ha='center', fontsize=7, fontweight='bold')
    format_clean_axis(ax, remove_ticks=False)


# =============================================================================
# Clock overlay (panels D and E)
# =============================================================================
# Real data and model on ONE axes, both in days, so the comparison is read off
# a single panel instead of across two. Replaces the old D/E (data) + F/G
# (model) split.
def _through_origin(x, y):
    """Slope of a through-origin fit, plus its R^2. Matches tempest_regression's
    fit_intercept=False so model and data rates are comparable."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    denom = np.sum(x ** 2)
    if denom <= 0:
        return np.nan, np.nan
    slope = np.sum(x * y) / denom
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - np.sum((y - slope * x) ** 2) / ss_tot if ss_tot > 0 else np.nan
    return slope, r2


def _thin(df, n=600, seed=42):
    return df.sample(n=n, random_state=seed) if len(df) > n else df


def plot_clock_overlay(ax, df_real, df_model, real_x, real_y, title,
                       data_label="Patient data", x_label="Days since infection",
                       x_max=400, y_max=None, model_scatter=False):
    """One clock, data and model together.

    Only the REAL cohort is scattered. The model contributes its fitted line
    only: it carries thousands of points against the data's hundreds, so its
    cloud buried the observations it is meant to be compared against. Set
    model_scatter=True to draw it anyway.
    """
    entries = []

    if df_real is not None and not df_real.empty:
        r = df_real.dropna(subset=[real_x, real_y]).copy()
        xr = r[real_x].astype(float).values
        yr = r[real_y].astype(float).values
        xr = xr - np.min(xr)                      # days since first observation
        if len(xr) > 1:
            slope, r2 = _through_origin(xr, yr)
            ax.scatter(xr, yr, s=7, alpha=0.28, color=DATA_COLOR,
                       edgecolors='none', zorder=2, label=data_label)
            xs = np.array([0.0, x_max])           # extend to the panel edge
            ax.plot(xs, slope * xs, color=DATA_COLOR, lw=1.7, zorder=7,
                    linestyle=DATA_STYLE)
            entries.append((data_label, slope * 365.25, r2))

    if df_model is not None and not df_model.empty:
        for cohort in [c for c in ("Control", "SOT", "HIV", "Edge Case")
                       if c in set(df_model.get('cohort', []))]:
            cdf = df_model[df_model['cohort'] == cohort]
            if cdf.empty:
                continue
            xm = cdf['Sequencing_time'].astype(float).values * 365.25
            ym = cdf['Distance_from_root'].astype(float).values
            slope, r2 = _through_origin(xm, ym)
            if not np.isfinite(slope):
                continue
            colour = SCEN_COLORS.get(cohort, "#888888")
            if model_scatter:
                t = _thin(cdf)
                ax.scatter(t['Sequencing_time'].astype(float) * 365.25,
                           t['Distance_from_root'].astype(float),
                           s=7, alpha=0.30, color=colour, edgecolors='none',
                           zorder=3, label=f"Model, {cohort}")
            xs = np.array([0.0, x_max])           # extend to the panel edge
            ax.plot(xs, slope * xs, color=colour, lw=1.7, zorder=6)
            entries.append((f"Model, {cohort}", slope * 365.25, r2))

    if not entries:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, f"{title}\n(Data Missing)", ha='center', va='center',
                transform=ax.transAxes)
        return

    if SHOW_FIT_STATS:
        txt = "\n".join(f"{n}: {rate:.5f} s/s/y ($R^2$={r2:.2f})"
                        for n, rate, r2 in entries)
        ax.text(0.03, 0.97, txt, transform=ax.transAxes, fontsize=5.4,
                va='top', ha='left',
                bbox=dict(facecolor='white', alpha=0.78, edgecolor='none', pad=1.2))

    ax.set_xlim(0, x_max)
    ax.set_ylim(0, y_max)
    ax.set_xlabel(x_label, fontsize=7)
    ax.set_ylabel("Divergence (subs/site)", fontsize=7)
    ax.text(0.5, 1.06, title, transform=ax.transAxes, ha='center',
            fontsize=7, fontweight='bold')
    format_clean_axis(ax, remove_ticks=False)
