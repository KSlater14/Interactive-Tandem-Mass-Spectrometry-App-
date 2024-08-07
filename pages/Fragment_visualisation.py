import streamlit as st
from streamlit_option_menu import option_menu
import plotly.graph_objects as go
import numpy as np

# Navigation menu
with st.sidebar:
    selected = option_menu(
        "Menu", 
        ["Home", "2D Fragment Visualization", "3D Fragment Visualization"], icons=["house", "scatter-plot", "graph-up"], 
        menu_icon="cast", 
        default_index=0)

# Navigation logic
if selected == "Home":
    from pages.Fragment_visualisation_pages import home
    home.show()
elif selected == "2D Fragment Visualization":
    from pages.Fragment_visualisation_pages import fragment_visualisation_2d as fragment_visualisation_2d
    fragment_visualisation_2d.show()
elif selected == "3D Fragment Visualization":
    from pages.Fragment_visualisation_pages import fragment_visualisation_3d as fragment_visualisation_3d
    fragment_visualisation_3d.show()

