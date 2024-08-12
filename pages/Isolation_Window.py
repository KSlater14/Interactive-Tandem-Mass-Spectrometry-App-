import numpy as np 
from pyteomics import mzml 
import streamlit as st
import requests
import os
from bokeh.plotting import figure 
from bokeh.models import ColumnDataSource
from bokeh.palettes import Spectral11

# Function to download the MRFA file if not already present
def load_MRFA_file(url, MRFA_filename):
    response = requests.get(url)
    if response.status_code == 200:
        with open(MRFA_filename, 'wb') as f:
            f.write(response.content)
    else:
        st.error("Error opening file")

# File path and URL
MRFA_file_path = "12Mar2024_MJ_ultramark1621.mzML"
mzml_MRFA_url = "https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_ultramark1621.mzML"

# Load file if not present
if not os.path.isfile(MRFA_file_path):
    load_MRFA_file(mzml_MRFA_url, MRFA_file_path)

# Function to average spectra
def average_spectra(spectra, bin_width=None, filter_string=None):
    reference_scan = np.unique(spectra[0]['m/z array'])
    if bin_width is None:
        bin_width = np.min(np.diff(reference_scan))

    scan_min = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window lower limit']
    scan_max = spectra[0]['scanList']['scan'][0]['scanWindowList']['scanWindow'][0]['scan window upper limit']

    reference_mz = np.arange(scan_min, scan_max, bin_width)
    merge_intensity = np.zeros_like(reference_mz)

    for scan in spectra:
        tmp_mz = scan['m/z array']
        tmp_intensity = scan['intensity array']
        merge_intensity += np.interp(reference_mz, tmp_mz, tmp_intensity, left=0, right=0)

    merge_intensity = merge_intensity / len(spectra)

    avg_spec = spectra[0].copy()
    avg_spec['m/z array'] = reference_mz
    avg_spec['intensity array'] = merge_intensity
    avg_spec['scanList']['scan'][0]['filter string'] = "AV: {:.2f}-{:.2f}; {}".format(
        spectra[0]['scanList']['scan'][0]['scan start time'], 
        spectra[-1]['scanList']['scan'][0]['scan start time'], 
        filter_string)

    return avg_spec


def plot_spectrum(MRFA_file, isolation_centre, isolation_width):
    reader = mzml.read(MRFA_file, use_index=True)
    spectra = list(reader)

    avg_spectrum = average_spectra(spectra)

    mzs = avg_spectrum["m/z array"]
    intensities = avg_spectrum["intensity array"]

    tophat_filter = tophat(isolation_centre, isolation_width, mzs)
    filtered_intensities = intensities * tophat_filter

    source_original = ColumnDataSource(data=dict(mzs=mzs, intensities=intensities))
    source_filtered = ColumnDataSource(data=dict(mzs=mzs, filtered_intensities=filtered_intensities))

    p = figure(title="Full MRFA MS1 Scan with Isolation Window",
               x_axis_label='m/z', y_axis_label='Intensity',
               x_range=(1000,1300))
    
    p.line('mzs', 'intensities', source=source_original, legend_label='Original Intensity', color=Spectral11[0], alpha=0.3, line_width=2)
    p.line('mzs', 'filtered_intensities', source=source_filtered, legend_label='Filtered Intensity', color=Spectral11[1], line_dash='dashed', line_width=2)

    p.legend.location = "top_left"
    p.legend.background_fill_color = 'white'

    return p


def tophat(centre, width, xs):
    ys = np.zeros_like(xs)
    lowerl = centre - (width / 2)
    upperl = centre + (width / 2)  # Corrected upper limit
    ys[(xs > lowerl) & (xs < upperl)] = 1
    return ys


# Streamlit layout
st.title("Isolation Window Visualization")

# Streamlit layout, creation of tabs 
Isolation_tab, Instruction_tab = st.tabs(["Isolation Window", "Instructions"])

with Instruction_tab:
    st.header("Instructions")
    st.markdown("Instructions for the use of the Isolation Window:")
    st.write("Upon selecting this Isolation Window page, a spectrum plot is generated of the full MS1 scan of MRFA with a default centrepoint of 1215.0 and width of 15.0.")
    st.write("To use this page, the parameters that can be adapted are the centrepoint of the isolation window and the width of the isolation window.")
    st.write("Simply use the addition and subtraction buttons on the input box to select a number for each of these parameters to identify how adjustment of these parameters influence the isolation window.")

with Isolation_tab:
    st.sidebar.title("Isolation Window Visualisation")
    st.sidebar.markdown("Explore the impact of the isolation window through changing the isolation centre and width.")

    isolation_centre = st.sidebar.number_input(
        "Select the centrepoint value (m/z) of the isolation window",
        min_value=1.0,
        max_value=2000.0,
        value=1215.0,
        step=0.1,
        format="%0.1f")

    isolation_width = st.sidebar.number_input(
        "Select the width of the isolation window",
        min_value=1.0,
        max_value=100.0,
        value=15.0,
        step=0.1,
        format="%0.1f")

    bokeh_figure = plot_spectrum(MRFA_file_path, isolation_centre, isolation_width)
    st.bokeh_chart(bokeh_figure)

