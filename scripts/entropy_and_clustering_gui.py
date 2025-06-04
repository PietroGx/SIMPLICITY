# This file is part of SIMPLICITY
# Copyright (C) 2025 Pietro Gerletti
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

import streamlit as st
import matplotlib.pyplot as plt
import simplicity.dir_manager as dm
import simplicity.plots_manager as pm 

# -----------------------------------------------------------------------------
# Import backend functions 
# -----------------------------------------------------------------------------
from temp_backend import (
    filter_seeds,
    build_seed_cache,
    prepare_figure1_data, prepare_figure2_data,
    prepare_figure3_data, prepare_figure4_data,
    prepare_figure5_data,
    plot_figure1, plot_figure2,
    plot_figure3, plot_figure4, plot_figure5,
    summarize_statistical_tests_from_dfs
)

pm.apply_plos_rcparams()
plt.rcParams.update({"font.size": 13})

# -----------------------------------------------------------------------------
# Sidebar Controls
# -----------------------------------------------------------------------------
st.set_page_config(layout="wide")
st.sidebar.header("SIMPLICITY Parameters")

min_time = st.sidebar.slider(
    "Min final time", min_value=0, max_value=1100, value=300, step=10
)
max_time = st.sidebar.slider(
    "Max final time", min_value=0, max_value=1100, value=450, step=10
)
cut_after_max = st.sidebar.checkbox("Cut after max time", value=True)

cluster_threshold = st.sidebar.slider(
    "Cluster threshold", min_value=0, max_value=10, value=5, step=1
)
takeover_threshold = st.sidebar.slider(
    "Takeover threshold", min_value=0.0, max_value=1.0, value=0.50, step=0.05
)

# -----------------------------------------------------------------------------
# Caching Filtered Seeds & Seed Cache (does not depend on seed)
# -----------------------------------------------------------------------------
@st.cache_data
def get_filtered_and_cached(min_time, max_time, cut_after_max, cluster_threshold):
    EXPERIMENT_NAME = "SIMPLICITY_exp_output"
    sim_dirs = dm.get_simulation_output_dirs(EXPERIMENT_NAME)
    sim_out_dir_1, sim_out_dir_2 = sim_dirs

    filtered_seeds_1 = filter_seeds(
        sim_out_dir_1, min_time, max_time, cut_after_max
    )
    filtered_seeds_2 = filter_seeds(
        sim_out_dir_2, min_time, max_time, cut_after_max
    )

    cache_1 = build_seed_cache(
        sim_out_dir_1,
        filtered_seeds_1,
        max_time=max_time,
        cut_after_max=cut_after_max,
        cluster_threshold=cluster_threshold
    )
    cache_2 = build_seed_cache(
        sim_out_dir_2,
        filtered_seeds_2,
        max_time=max_time,
        cut_after_max=cut_after_max,
        cluster_threshold=cluster_threshold
    )
    return (
        sim_out_dir_1,
        sim_out_dir_2,
        filtered_seeds_1,
        filtered_seeds_2,
        cache_1,
        cache_2
    )

(
    sim_dir_1,
    sim_dir_2,
    filtered_seeds_1,
    filtered_seeds_2,
    cache_1,
    cache_2
) = get_filtered_and_cached(
    min_time, max_time, cut_after_max, cluster_threshold
)

# -----------------------------------------------------------------------------
# Seed Dropdown (changes trigger only Figure 1 recompute)
# -----------------------------------------------------------------------------
filtered_seeds = sorted(set(filtered_seeds_1) | set(filtered_seeds_2))
if not filtered_seeds:
    st.warning("No valid seeds in the given (min_time, max_time, cut_after_max) range.")
    st.stop()

selected_seed = st.sidebar.selectbox(
    "Select seed for Figure 1",
    options=filtered_seeds,
    index=0
)

# -----------------------------------------------------------------------------
# Cached preparations for Figures 2–5
# -----------------------------------------------------------------------------
@st.cache_data
def get_figure2_data(
    sim_out_dir_1, sim_out_dir_2,
    filtered_seeds_1, filtered_seeds_2,
    cache_1, cache_2, takeover_threshold
):
    return prepare_figure2_data(
        sim_out_dir_1=sim_out_dir_1,
        sim_out_dir_2=sim_out_dir_2,
        filtered_seeds_1=filtered_seeds_1,
        filtered_seeds_2=filtered_seeds_2,
        cache_1=cache_1,
        cache_2=cache_2,
        takeover_threshold=takeover_threshold
    )

@st.cache_data
def get_figure3_data(
    sim_out_dir_1, sim_out_dir_2,
    filtered_seeds_1, filtered_seeds_2,
    cache_1, cache_2, max_time
):
    return prepare_figure3_data(
        sim_out_dir_1=sim_out_dir_1,
        sim_out_dir_2=sim_out_dir_2,
        filtered_seeds_1=filtered_seeds_1,
        filtered_seeds_2=filtered_seeds_2,
        cache_1=cache_1,
        cache_2=cache_2,
        max_time=max_time
    )

@st.cache_data
def get_figure4_data(
    sim_out_dir_1, sim_out_dir_2,
    filtered_seeds_1, filtered_seeds_2,
    cache_1, cache_2
):
    return prepare_figure4_data(
        sim_out_dir_1=sim_out_dir_1,
        sim_out_dir_2=sim_out_dir_2,
        filtered_seeds_1=filtered_seeds_1,
        filtered_seeds_2=filtered_seeds_2,
        cache_1=cache_1,
        cache_2=cache_2
    )

@st.cache_data
def get_figure5_data(
    sim_out_dir_1, sim_out_dir_2,
    filtered_seeds_1, filtered_seeds_2,
    cache_1, cache_2
):
    return prepare_figure5_data(
        sim_out_dir_1=sim_out_dir_1,
        sim_out_dir_2=sim_out_dir_2,
        filtered_seeds_1=filtered_seeds_1,
        filtered_seeds_2=filtered_seeds_2,
        cache_1=cache_1,
        cache_2=cache_2
    )

fig2_data = get_figure2_data(
    sim_out_dir_1=sim_dir_1,
    sim_out_dir_2=sim_dir_2,
    filtered_seeds_1=filtered_seeds_1,
    filtered_seeds_2=filtered_seeds_2,
    cache_1=cache_1,
    cache_2=cache_2,
    takeover_threshold=takeover_threshold
)

fig3_data = get_figure3_data(
    sim_out_dir_1=sim_dir_1,
    sim_out_dir_2=sim_dir_2,
    filtered_seeds_1=filtered_seeds_1,
    filtered_seeds_2=filtered_seeds_2,
    cache_1=cache_1,
    cache_2=cache_2,
    max_time=max_time
)

fig4_data = get_figure4_data(
    sim_out_dir_1=sim_dir_1,
    sim_out_dir_2=sim_dir_2,
    filtered_seeds_1=filtered_seeds_1,
    filtered_seeds_2=filtered_seeds_2,
    cache_1=cache_1,
    cache_2=cache_2
)

fig5_data = get_figure5_data(
    sim_out_dir_1=sim_dir_1,
    sim_out_dir_2=sim_dir_2,
    filtered_seeds_1=filtered_seeds_1,
    filtered_seeds_2=filtered_seeds_2,
    cache_1=cache_1,
    cache_2=cache_2
)

# -----------------------------------------------------------------------------
# Figure 1 Preparation (not cached by seed)
# -----------------------------------------------------------------------------
fig1_data = prepare_figure1_data(
    sim_out_dir_1=sim_dir_1,
    sim_out_dir_2=sim_dir_2,
    seed=selected_seed,
    cache_1=cache_1,
    cache_2=cache_2
)

# -----------------------------------------------------------------------------
# Compute “Test Table Summary” (Mann–Whitney results for Figs 2, 4, 5)
# -----------------------------------------------------------------------------
stats_table = summarize_statistical_tests_from_dfs(
    fig2_data["summary_df"],
    fig4_data["summary_df"],
    fig5_data["summary_df"]
)

# -----------------------------------------------------------------------------
# Layout with Three Columns
# -----------------------------------------------------------------------------
col1, col2, col3 = st.columns(3)

with col1:
    st.subheader("Simulation Entropy")
    # Render at 20×20, Streamlit will resize it to the column width
    fig1, axs1 = plt.subplots(2, 2, figsize=(5, 5), sharex="col")
    plot_figure1(axs1, fig1_data)
    fig1.subplots_adjust(top=0.95, bottom=0.05)
    st.pyplot(fig1)

    st.subheader("Takeover Events")
    fig2, axs2 = plt.subplots(2, 1, figsize=(5, 5), sharex=False)
    plot_figure2(axs2, fig2_data)
    fig2.subplots_adjust(top=0.95, bottom=0.05)
    st.pyplot(fig2)

with col2:
    st.subheader("Mean Entropy")
    fig3, axs3 = plt.subplots(2, 2, figsize=(5, 5), sharex="col")
    plot_figure3(axs3, fig3_data)
    fig3.subplots_adjust(top=0.95, bottom=0.05)
    st.pyplot(fig3)

    st.subheader("MW Tests Summary")
    st.table(stats_table)

with col3:
    st.subheader("Max Entropy")
    fig4, axs4 = plt.subplots(2, 1, figsize=(5, 5), sharex=False)
    plot_figure4(axs4, fig4_data)
    fig4.subplots_adjust(top=0.95, bottom=0.05)
    st.pyplot(fig4)

    st.subheader("Last-time Entropy")
    fig5, axs5 = plt.subplots(2, 1, figsize=(5, 5), sharex=False)
    plot_figure5(axs5, fig5_data)
    fig5.subplots_adjust(top=0.95, bottom=0.05)
    st.pyplot(fig5)
