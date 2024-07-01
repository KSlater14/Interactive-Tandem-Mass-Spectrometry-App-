import streamlit as st
import numpy as np
import pandas as pd
from bokeh.plotting import figure
from bokeh.models import HoverTool, ColumnDataSource, LabelSet
from pyteomics import mzml
import requests
import io
from scipy import signal

# Define the peak detection function
def peak_detection(spectrum, threshold=5, distance=4, prominence=0.8, width=2, centroid=False):
    relative_threshold = spectrum['intensity array'].max() * (threshold / 100)
    if centroid:
        peaks = np.where(spectrum['intensity array'] > relative_threshold)[0]
        return peaks, None
    else:
        peaks, properties = signal.find_peaks(spectrum['intensity array'], height=relative_threshold, prominence=prominence, width=width, distance=distance)
        return peaks, properties

# Define the function to get centroid values
def get_centroid(spectrum, peaks, properties):
    centroids = np.zeros_like(peaks, dtype='float32')
    mz_array = spectrum['m/z array']
    intensity_array = spectrum['intensity array']
    for i, peak in enumerate(peaks):
        left_index = int(properties['left_ips'][i])
        right_index = int(properties['right_ips'][i])
        peak_range = range(left_index, right_index + 1)
        mz_values = mz_array[peak_range]
        intensity_values = intensity_array[peak_range]
        weighted_mz_sum = np.sum(mz_values * intensity_values)
        total_intensity = np.sum(intensity_values)
        centroids[i] = weighted_mz_sum / total_intensity 
    return centroids

# Define the function to load mzML data
def load_mzml_data(peptide):
    file_map = {
        'MRFA': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_full_scan_enhanced.mzML',
        'Bradykinin': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin.mzML',
        'GRGDS': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_full_ms.mzML',
        'SDGRG': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_full_ms.mzML'
    }
    file_url = file_map.get(peptide)
    if file_url:
        response = requests.get(file_url)
        mzml_data = io.BytesIO(response.content)
        spectra = mzml.read(mzml_data)
        return spectra
    else:
        st.error("No mzML data available for the selected peptide.")
        return None
    
# Streamlit layout
Spectrum_tab, Instruction_tab = st.tabs(["Spectrum", "Instructions"])


with Instruction_tab:  

    image_peptide_selection = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/peptide%20selection.png'
    image_full_peptide_plot = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/plotting%20of%20peptide.png'
    image_plot_expansion = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/View%20fullscreen.png'
    image_zoom_function =  'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Zoom%20function.png'
    image_hover_function = 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Instruction%20images/Hover%20function.png'

    st.header("Instructions")
    st.markdown("Instructions for the use of the Peptide Explorer")
    st.subheader("Peptide Selection")
    st.write("""
- Within this interactive peptide explorer, the entire spectra for a single peptide can be explored. 
- The user can select which peptide they would like to investigate. 
- Once selected the spectra is automatically plotted.
             """)

    col1, col2 = st.columns(2)
    with col1: 
        st.image(image_peptide_selection, caption='Data Selection: Using data already available', width=300)
    with col2:
        st.image(image_full_peptide_plot, caption='Automatic plotting of the full spectrum specific to each peptide', width=900)    

    st.subheader("Plot interactivity")
    st.write("""
    - Once plotted, the only settings that can be changed is whether the user wants to 'Show m/z Labels', 
             which can be selected via the checkbox. 
    - The plot has various interactive features that allows the user to undertake a thorough exploration. These features include:
   
    - The expansion of the plot to a full screen.""")
    st.image(image_plot_expansion, caption='Button for spectrum plot expansion', width=800)
    ("""
        - The ability to drag the cursor of a section of the plot to zoom in for exploration.""")
    st.image(image_zoom_function, caption='The cursor drag zoom function', width=600)
    ("""
        - A hover tool allows the cursor, when hovering over a peak, to display the m/z, intensity and centroid data regarding each peak.
             """)
    st.image(image_hover_function, caption='The hover function in use', width=600)

    
with Spectrum_tab:
    st.markdown("Explore the parameters influencing the spectra, over a series of scans. Select instructions for help.") 

# Define the function to plot the spectrum

def plot_spectrum(spectrum, show_labels):
    mz_values = spectrum['m/z array']
    intensity_values = spectrum['intensity array']

    peaks, properties = peak_detection(spectrum, threshold=5, centroid=False)
    peak_centroids = get_centroid(spectrum, peaks, properties)

    source = ColumnDataSource(data=dict(
        x=mz_values[peaks],
        y=intensity_values[peaks],
        cent=["%.2f" % x for x in peak_centroids]
    ))

    TOOLTIPS = [
                ("m/z", "@x{0.00}"),
                ("intensity", "@y{0.0}"),
                ("centroid", "@cent{0.00}")
    ]

    p = figure(title="", x_axis_label="m/z", y_axis_label="Intensity", tooltips=None,
            tools='pan,box_zoom,xbox_zoom, reset,save', active_drag='xbox_zoom')
        
    p.line(mz_values, intensity_values, line_width=1, color='black')
    r = p.circle('x', 'y', size=5, source=source, color='red')

    if show_labels:
        labels = LabelSet(x='x', y='y', text='cent', source=source, text_font_size='8pt', text_color='black')
        p.add_layout(labels)

    p.x_range.start = min(mz_values)
    p.x_range.end = max(mz_values)
    p.y_range.start = 0
    p.y_range.end = max(intensity_values) * 1.1

    hover = HoverTool(tooltips=TOOLTIPS, renderers=[r])
    p.add_tools(hover)
    st.bokeh_chart(p, use_container_width=True)


st.sidebar.title("Interactive Peptide Explorer")
st.sidebar.markdown("Select a peptide to explore its mass spectrometry data.")

peptide_options = ["MRFA", "Bradykinin", "GRGDS", "SDGRG"]
selected_peptide = st.sidebar.selectbox("Select Peptide", peptide_options)

show_labels = st.sidebar.checkbox("Show m/z Labels", value=False)

with Spectrum_tab:
    spectra = load_mzml_data(selected_peptide)
    if spectra is not None:
        first_spectrum = next(spectra, None)  # Get the first spectrum or None if no spectra
        if first_spectrum and 'm/z array' in first_spectrum and 'intensity array' in first_spectrum:
            plot_spectrum(first_spectrum, show_labels)
            mz_values = first_spectrum['m/z array']
            intensity_values = first_spectrum['intensity array']
        else:
            st.error("The mzML data for the selected peptide is not in the expected format.")
        