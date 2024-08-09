import numpy as np 
from pyteomics import mzml 
import streamlit as st
import requests
import os
from bokeh.plotting import figure 
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral11

def load_MRFA_file(url, MRFA_filename):
    response = requests.get(url)
    if response.status_code == 200:
        with open(MRFA_filename, 'wb') as f:
            f.write(response.content)
    else:
        st.error("Error opening file")

MRFA_file_path = "12Mar2024_MJ_ultramark1621.mzML"
mzml_MRFA_url = "https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_ultramark1621.mzML"

if not os.path.isfile(MRFA_file_path):
    load_MRFA_file(mzml_MRFA_url, MRFA_file_path)

def plot_spectrum(MRFA_file, spectrum_index, isolation_centre, isolation_width):
    reader = mzml.read(MRFA_file, use_index=True)
    spec = reader[spectrum_index]

    mzs = spec["m/z array"]
    intensities = spec["intensity array"]

    tophat_filter = tophat(isolation_centre, isolation_width, mzs)
    filtered_intensities = intensities * tophat_filter

    source_original = ColumnDataSource(data=dict(mzs=mzs, intensities=intensities))
    source_filtered = ColumnDataSource(data=dict(mzs=mzs, filtered_intensities=filtered_intensities))

    p = figure(title="Full MRFA MS1 Scan with Isolation Window",
               x_axis_label='m/z', y_axis_label='Intensity',
               x_range=(1000,1300))
    
    p.line('mzs', 'intensities', source=source_original, legend_label='Original Intensity', color=Spectral11[0])
    p.line('mzs', 'filtered_intensities', source=source_filtered, legend_label='Filtered Intensity', color=Spectral11[1], line_dash='dashed')

    p.legend.location = "top_left"
    p.legend.background_fill_color = 'white'

    return p

def tophat(centre, width, xs):
    ys = np.zeros_like(xs)
    lowerl = centre - (width / 2)
    upperl = centre - (width / 2)
    ys[(xs > lowerl) & (xs < upperl)] = 1 

    return ys

# Streamlit layout
st.title("Isolation Window Visualization")

# Streamlit layout, creation of tabs 
Isolation_tab, Instruction_tab = st.tabs(["Isolation Window", "Instructions"])

with Isolation_tab:
    st.sidebar.title("Isolation Window Visualisation")
    st.sidebar.markdown("Explore the impact of the isolation window.... ")

    isolation_centre = st.sidebar.slider(
        "Select the Centrepoint value (m/z) of the isolation window",
        min_value=1.0,
        max_value=2000.0,
        value=1215.0)

    isolation_width = st.sidebar.slider(
        "Select the width of the isolation window",
        min_value=1.0,
        max_value=100.0,
        value=15.0)

spectrum_index = st.sidebar.slider(
        "Select spectrum index",
        min_value=0,
        max_value=100,
        value=20
    )

bokeh_figure = plot_spectrum(MRFA_file_path, spectrum_index, isolation_centre, isolation_width)
st.bokeh_chart(bokeh_figure)

with Instruction_tab:
    st.header("Instructions")
    st.markdown("Instructions for the use of the Isolation Window")
    st.subheader("Isolation Window")
    st.write("""

""")

