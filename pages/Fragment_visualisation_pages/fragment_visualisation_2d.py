import streamlit as st
import streamlit.components.v1 as components
import plotly.graph_objects as go
import numpy as np
from pyteomics import mzml, parser, mass
import requests
import io
from scipy import signal
import pandas as pd
import py3Dmol

aa_mass = mass.std_aa_mass
aa_mass['p'] = 79.966331  # phosphorylation (STY)
aa_mass['ox'] = 15.994915  # oxidation (MW)
aa_mass['d'] = 0.984016  # deamidation (NQ)
aa_mass['am'] = -0.984016  # amidation (C-term)

# Peak detection function
def peak_detection(spectrum, threshold=5, distance=4, prominence=0.8, width=2, centroid=False):
    relative_threshold = spectrum['intensity array'].max() * (threshold / 100)
    if centroid:
        peaks = np.where(spectrum['intensity array'] > relative_threshold)[0]
        return peaks
    else:
        peaks, properties = signal.find_peaks(spectrum['intensity array'], height=relative_threshold, prominence=prominence, width=width, distance=distance)
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
                        fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': ion_type})
                        print(f"Annotated fragment: {ion_label}, m/z: {_mass}")

    # precursor ion annotation
    for charge in range(1, selected_charge_state + 1):
        for loss, mass_diff in neutral_losses.items():
            seq = ''.join(_sequence)
            _mass = mass.fast_mass2(seq, ion_type="M", charge=charge, aa_mass=aa_mass) + (mass_diff / charge)
            ion_label = "M" + loss + "+"*charge

            if is_in_peaks_data(_mass):
                fragments.append({'seq': seq, 'ion': ion_label, 'm/z': _mass, 'type': "M"})
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
    # Extract m/z values, ion types, and ion labels
    mz_values = [fragment['m/z'] for fragment in fragments]
    ion_labels = [fragment['ion'] for fragment in fragments]
    ion_types = [fragment['type'] for fragment in fragments]
   
 

    # Define colors for different ion types
    colors = {'b': 'blue', 'y': 'red', 'M': 'yellow', '': 'black'}
    unique_ion_types = set(ion_types)

    # Create a Plotly figure
    fig = go.Figure()

    # Plot each fragment with different colors for different ion types
    for ion_type in unique_ion_types:
        indices = [i for i, t in enumerate(ion_types) if t == ion_type]
        fig.add_trace(go.Scatter(
            x=[mz_values[i] for i in indices],
            y=[ion_labels[i] for i in indices],
            mode='markers',
            marker=dict(color=colors.get(ion_type, 'black')),
            name=f'Ion Type {ion_type}'
        ))

    # Update layout
    fig.update_layout(
        title='Annotated Fragments',
        xaxis_title='m/z',
        yaxis_title='Ion',
        legend_title='Ion Type'
    )

    # Display the figure in Streamlit
    st.plotly_chart(fig)


def show():
    peptide_options = {
        'MRFA': 'MRFA',
        'Bradykinin': 'RPPGFSPFR',
        'GRGDS': 'GRGDS',
        'SDGRG': 'SDGRG'
    }

    selected_peptide_name = st.sidebar.selectbox(
        "Select a peptide sequence", list(peptide_options.keys())
    )

    peptide_sequence = peptide_options[selected_peptide_name]
    st.write(f"Selected Peptide Sequence: {peptide_sequence}")

    # User-selectable charge state
    selected_charge_state = st.slider("Select Charge State", min_value=1, max_value=3, value=1)

    # User-selectable isolation window
    isolation_window = st.slider(
        "Select Isolation Window (m/z)", 
        min_value=0.0, 
        max_value=2000.0, 
        value=(0.0, 2000.0)
    )
    st.write(f"Isolation Window: {isolation_window}")

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
        
        df_fragments = pd.DataFrame(fragments)
        st.table(df_fragments[['type', 'm/z', 'ion']])
    
       
        


       
if __name__ == "__main__":
    show()
