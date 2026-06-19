# nature_figures_plots.py
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from scipy import stats
import numpy as np
import pandas as pd
import seaborn as sns

def format_clean_axis(ax, remove_ticks=True):
    """Applies Nature-standard clean aesthetics to an axis."""
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    if remove_ticks:
        ax.set_xticks([])
        ax.set_yticks([])

def plot_baseline_dynamics(ax, df, title):
    format_clean_axis(ax, remove_ticks=True)
    ax.text(0.5, 0.5, f"{title}\nPlaceholder", ha='center', va='center', transform=ax.transAxes)

def plot_infection_duration(ax, df):
    format_clean_axis(ax, remove_ticks=True)
    ax.text(0.5, 0.5, "Infection Duration\n(Real Data Toy)\nPlaceholder", ha='center', va='center', transform=ax.transAxes)

def plot_tempest_regression(ax, df, title, x_col='sampling_date', y_col='hamming_per_overlap', is_date=True, scatter_color='#404040', line_color='black', force_zero_intercept=False, force_intercept_at_min_x=False):
    """Plots root-to-tip genetic distance vs time with a linear regression line."""
    if df.empty or x_col not in df.columns or y_col not in df.columns:
        format_clean_axis(ax, remove_ticks=True)
        ax.text(0.5, 0.5, f"{title}\n(Data Missing)", ha='center', va='center', transform=ax.transAxes)
        return

    # Drop missing values and sort
    df = df.dropna(subset=[x_col, y_col]).copy()
    df = df.sort_values(x_col)
    
    y = df[y_col].values

    # Handle Date vs Numeric X-axis
    if is_date:
        df[x_col] = pd.to_datetime(df[x_col])
        x_dates = df[x_col]
        x_plot = mdates.date2num(x_dates.dt.to_pydatetime())
    else:
        x_plot = df[x_col].astype(float).values
        x_dates = x_plot

    # --- REGRESSION LOGIC ---
    if force_zero_intercept:
        # 1. Force intercept to exactly 0 (y = m*x)
        slope = np.sum(x_plot * y) / np.sum(x_plot**2)
        y_pred = slope * x_plot
        
        # Calculate R^2 manually
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        # Plot manual line
        x_line = np.array([0, np.max(x_plot)])
        ax.plot(x_line, slope * x_line, color=line_color, linewidth=1.2, zorder=1)
        
    elif force_intercept_at_min_x:
        # 1. Force intercept to 0 at the MINIMUM X value (day of first sample)
        x_min = np.min(x_plot)
        x_shifted = x_plot - x_min  # Shift X so the first sample is day 0
        
        # Calculate zero-intercept slope on the shifted data
        slope = np.sum(x_shifted * y) / np.sum(x_shifted**2)
        y_pred = slope * x_shifted
        
        # Calculate R^2 manually
        ss_res = np.sum((y - y_pred)**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        # Plot manual line (Shift it back to real calendar dates for the graph)
        x_line = np.array([x_min, np.max(x_plot)])
        y_line = slope * (x_line - x_min)
        ax.plot(x_line, y_line, color=line_color, linewidth=1.2, zorder=1)

    else:
        # Standard least-squares regression (includes free intercept)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x_plot, y)
        r_squared = r_value**2
        
        # Plot using Seaborn for the 95% Confidence Interval band
        sns.regplot(x=x_plot, y=y, ax=ax, scatter=False, color=line_color, ci=95, 
                    line_kws={'linewidth': 1.2})

    # --- SCATTER PLOT ---
    ax.scatter(x_plot, y, alpha=0.5, s=8, color=scatter_color, edgecolors='none', zorder=3)

    # --- ANNOTATIONS & FORMATTING ---
    rate_per_year = slope * 365.25 

    stats_text = f"$R^2 = {r_squared:.2f}$\nRate $= {rate_per_year:.5f}$ s/s/y"
    ax.text(0.05, 0.95, stats_text, transform=ax.transAxes, 
            fontsize=6, va='top', ha='left', 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none', pad=1))

    if is_date:
        ax.xaxis.set_major_locator(mdates.YearLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
        for label in ax.get_xticklabels():
            label.set_rotation(45)
            label.set_ha('right')
        ax.set_xlim(left=np.min(x_plot), right=np.max(x_plot))
    else:
        ax.set_xlabel("Days since infection", fontsize=7)
        ax.set_xlim(left=0, right=np.max(x_plot) * 1.05)

    ax.set_ylim(bottom=0)
    ax.set_ylabel("Root-to-tip divergence", fontsize=7)
    ax.text(0.5, 1.05, title, transform=ax.transAxes, ha='center', fontsize=7, fontweight='bold')
    
    format_clean_axis(ax, remove_ticks=False)