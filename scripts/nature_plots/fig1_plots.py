import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
from scipy import stats
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter

def format_clean_axis(ax, remove_ticks=True):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

def plot_fig1_intra_host(ax, df):
    if df.empty:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, "Panel A\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return
    palette = {"Control": "#333333", "SOT": "#56B4E9", "HIV": "#D55E00", "Edge Case": "#CC79A7"}
    sns.lineplot(data=df, x='time', y='p_infectious', hue='cohort', palette=palette, ax=ax, linewidth=1.5, alpha=0.8, legend=False)
    ax.set_xlim(0, 800) 
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
    palette = {"Control": "#333333", "SOT": "#56B4E9", "HIV": "#D55E00", "Edge Case": "#CC79A7"}
    order = ["Control", "SOT", "HIV", "Edge Case"]
    existing_order = [c for c in order if c in df['cohort'].unique()]
    sns.violinplot(data=df, x='cohort', y='duration', palette=palette, hue='cohort', order=existing_order, ax=ax, inner="quartile", cut=0, linewidth=0.8, density_norm='width', legend=False)
    ax.set_ylim(0, 800)
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
    
    box_palette = {}
    for cat in sorted_categories:
        if cat == 'HIV': box_palette[cat] = "#D55E00"
        elif cat == 'SOT': box_palette[cat] = "#56B4E9"
        else: box_palette[cat] = "#B0B0B0"
        
    sns.boxplot(data=df, x='clinical_category', y='duration', order=sorted_categories, ax=ax, 
                palette=box_palette, hue='clinical_category', showfliers=False, width=0.5, boxprops=dict(alpha=0.4), legend=False)
    
    def get_pt_color(row):
        if row['duration'] > 300: return "#CC79A7"
        return box_palette[row['clinical_category']]
    df['pt_color'] = df.apply(get_pt_color, axis=1)

    for i, cat in enumerate(sorted_categories):
        cat_data = df[df['clinical_category'] == cat]
        x_jitter = np.random.normal(i, 0.05, size=len(cat_data))
        ax.scatter(x_jitter, cat_data['duration'], c=cat_data['pt_color'], s=12, alpha=0.8, zorder=5)

    ax.set_ylim(0, 800)
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
