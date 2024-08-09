## FRAGMENT VISUALISATION 2D ##
"""
This script allows for the 2D visualisation of the fragments found in each peptide.

"""
import streamlit as st
import streamlit.components.v1 as components
from bokeh.models import ColumnDataSource, HoverTool, Legend, LegendItem
from bokeh.plotting import figure
from bokeh.palettes import Category10
import numpy as np
from pyteomics import mzml, parser, mass
import requests
import io
import math
from scipy import signal
import pandas as pd

## FUNCTIONS ##

# Peak detection function
def peak_detection(spectrum, threshold=5, distance=4, prominence=0.8, width=2, centroid=False):
    relative_threshold = spectrum['intensity array'].max() * (threshold / 100)
    if centroid:
        peaks = np.where(spectrum['intensity array'] > relative_threshold)[0]
        return peaks
    else:
        peaks, properties = signal.find_peaks(
            spectrum['intensity array'], 
            height=relative_threshold, 
            prominence=prominence, 
            width=width, 
            distance=distance)
        return peaks, properties

# Centroid calculation function
def return_centroid(spectrum, peaks, properties):
    centroids = np.zeros_like(peaks, dtype='float32')
    for i, peak in enumerate(peaks):
        left_ip = int(properties['left_ips'][i])
        right_ip = int(properties['right_ips'][i])
        peak_range = range(left_ip, right_ip + 1)
        mz_range = spectrum['m/z array'][peak_range]
        intensity_range = spectrum['intensity array'][peak_range]
        centroids[i] = np.sum(mz_range * intensity_range) / np.sum(intensity_range)
    return centroids

aa_mass = mass.std_aa_mass
aa_mass['p'] = 79.966331  # phosphorylation (STY)
aa_mass['ox'] = 15.994915  # oxidation (MW)
aa_mass['d'] = 0.984016  # deamidation (NQ)
aa_mass['am'] = -0.984016  # amidation (C-term)

# Fragment annotation function with isolation window
def get_fragments(sequence, selected_charge_state, peaks_data, isolation_window, ion_types=('b', 'y')):
    fragments = []
    _sequence = parser.parse(sequence)
    pep_length = len(_sequence)
    print(f"Peptide length: {pep_length}")

    neutral_losses = {
        '': 0, 
        '-H2O': -18.01528, 
        '-NH3': -17.02655
    }

    def is_in_peaks_data(mass):
        tolerance = 0.2
        lower_bound = isolation_window[0] - tolerance
        upper_bound = isolation_window[1] + tolerance
        return any(lower_bound <= peak <= upper_bound for peak in peaks_data)

    for pos in range(1, pep_length):
        for ion_type in ion_types:
            if ion_type[0] in ('a', 'b', 'c'):
                seq = ''.join(_sequence[:pos])
            elif ion_type[0] in ('x', 'y', 'z'):
                seq = ''.join(_sequence[-pos:])
            
            for charge in range(1, selected_charge_state + 1):
                for loss, mass_diff in neutral_losses.items():
                    _mass = mass.fast_mass2(seq, ion_type=ion_type, charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
                    ion_label = ion_type + str(pos) + loss + "+"*charge

                    if is_in_peaks_data(_mass):
                        fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type, 'start': _mass - 0.5, 'end': _mass + 0.5})
                        print(f"Annotated fragment: {ion_label}, m/z: {_mass}")


    # Precursor ion annotation
    for charge in range(1, selected_charge_state + 1):
        for loss, mass_diff in neutral_losses.items():
            seq = ''.join(_sequence)
            _mass = mass.fast_mass2(seq, ion_type="M", charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
            ion_label = "M" + loss + "+"*charge

            if is_in_peaks_data(_mass):
                fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': "M", 'start': _mass - 0.5, 'end': _mass + 0.5})
                print(f"Annotated fragment: {ion_label}, m/z: {_mass}")
    
    return fragments

# Function to load mzML data
def load_mzml_data(peptide):
    file_map = {
        'MRFA': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/MRFA/12Mar2024_MJ_MRFA_full_scan_enhanced.mzML',
        'Bradykinin': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/Bradykinin/13Mar2024_MJ_bradykinin.mzML',
        'GRGDS': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/GRGDS/21Mar2024_MJ_GRGDS_full_ms.mzML',
        'SDGRG': 'https://raw.githubusercontent.com/KSlater14/Interactive-Tandem-Mass-Spectrometry-App-/main/Data/SDGRG/21Mar2024_MJ_SDGRG_full_ms.mzML'
    }
    file_url = file_map.get(peptide)
    if file_url:
        try:
            response = requests.get(file_url)
            if response.status_code == 200:
                mzml_data = io.BytesIO(response.content)
                spectra = list(mzml.read(mzml_data))
                return spectra
            else:
                print(f"Failed to fetch data. HTTP Status Code: {response.status_code}")
                return None
        except Exception as e:
            print(f"Error loading mzML data: {e}")
            return None
    else:
        print("No mzML data available for the selected peptide.")
        return None



# Plot fragments function
def plot_fragments(fragments):
    if not fragments:
        st.write("No fragments found to plot.")
        return
    
    # Extract m/z values, ion types, and ion labels from the fragments
    mz_values = [fragment['m/z'] for fragment in fragments]
    ion_labels = [fragment['ion'] for fragment in fragments]
    ion_types = [fragment['type'] for fragment in fragments]
    
    ion_start = [fragment.get('start', fragment['m/z'] - 0.5) for fragment in fragments]
    ion_end = [fragment.get('end', fragment['m/z'] + 0.5) for fragment in fragments]

    # Define colors for different ion types
    colors = {'b': 'blue', 'y': 'red', 'M': 'orange', '': 'black'}
    color_values = [colors.get(ion_type, 'black') for ion_type in ion_types]

    # Create a ColumnDataSource from the fragment data
    fragment_data = ColumnDataSource(data=dict(
        mz_values=mz_values,
        ion_labels=ion_labels,
        ion_start=ion_start,
        ion_end=ion_end,
        ion_types=ion_types,
        color_values=color_values  # Add color values to the data source
    ))
    
    # Create a Bokeh figure for the plot
    p = figure(
        title="Fragment Ion Visualisation",
        x_axis_label='m/z',
        y_axis_label='Ion',
        y_range=fragment_data.data['ion_labels'],
        tools='pan,box_zoom,xbox_zoom,reset,save',
        active_drag='xbox_zoom'
    )

    p.hbar( 
        y='ion_labels',
        left='ion_start',
        right='ion_end',
        height=0.7,
        color='color_values',
        source=fragment_data
    )

    if mz_values:
        max_mz = math.ceil(max(mz_values))  # Round up to nearest integer
        x_ticks = list(range(0, max_mz + 1, 10))  # Create ticks at every 10 units
        p.xaxis.ticker = x_ticks
        p.xaxis.major_label_overrides = {tick: str(tick) for tick in x_ticks}
    else:
        p.xaxis.ticker = [0]

    # Add HoverTool for displaying tooltips
    hover = HoverTool(tooltips=[
        ("m/z", "@mz_values"),
        ("Ion", "@ion_labels"),
        ("Type", "@ion_types")
    ])
    p.add_tools(hover)

    legend_items = []
    unique_ion_types = set(ion_types)  # Get unique ion types to avoid duplicate legend items
    for ion_type in unique_ion_types:
        color = colors.get(ion_type, 'black')
        
        # Create a dummy plot for the legend item
        legend_source = ColumnDataSource(data=dict(x=[0], y=[0]))
        
        # Create a glyph for the legend item
        legend_glyph = p.scatter(
            x='x',
            y='y',
            color=color,
            size=8,
            source=legend_source,
            legend_label=ion_type
        )
        
        legend_items.append(LegendItem(label=f"{ion_type}", renderers=[legend_glyph]))

    # Configure the appearance of the legend
    p.legend.title = 'Ion Type'
    p.legend.location = 'right'
 
    st.bokeh_chart(p, use_container_width=True)

# Function to display the Streamlit app interface 
def show():
    peptide_options = {
        'MRFA': 'MRFA',
        'Bradykinin': 'RPPGFSPFR',
        'GRGDS': 'GRGDS',
        'SDGRG': 'SDGRG'
    }

    selected_peptide_name = st.sidebar.selectbox(
        "Select a peptide sequence", 
        list(peptide_options.keys())
    )

    peptide_sequence = peptide_options[selected_peptide_name]
    st.write(f"Selected Peptide Sequence: {peptide_sequence}")

    # User-selectable slider for charge state
    selected_charge_state = st.slider(
        "Select Charge State", 
        min_value=1, 
        max_value=3, 
        value=1)

    # User-selectable slider for isolation window
    isolation_window = (0.0, 2000.0)
    
    # Load mzML data
    mzml_data = load_mzml_data(selected_peptide_name)
    if mzml_data is not None:
        spectrum = mzml_data[0]
        peaks, properties = peak_detection(spectrum)


        # Get centroid m/z values for the peaks
        peaks_data = return_centroid(spectrum, peaks, properties)
        

        # Annotate fragments based on peaks data and isolation window
        fragments = get_fragments(peptide_sequence, selected_charge_state, peaks_data, isolation_window)
        plot_fragments(fragments)
        
        # Converts fragments data to a DataFrame for display 
        df_fragments = pd.DataFrame(fragments)
        # Display the fragment in a table, with columns as specified 
        st.table(df_fragments[['type', 'm/z', 'ion']])

       
if __name__ == "__main__":
    show()
