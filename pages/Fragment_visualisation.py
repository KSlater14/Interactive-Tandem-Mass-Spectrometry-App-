## FRAGMENT VISUALISATION MENU PAGE ##
"""
This script produces a Fragment visualisation menu to access different types of fragment 
visualisation to meet the users requirements. 
"""
import streamlit as st
from streamlit_option_menu import option_menu
import plotly.graph_objects as go
import numpy as np

# Sidebar navigation menu setup
with st.sidebar:
    selected = option_menu(
        "Menu", 
        ["Home", "2D Fragment Visualization", "3D Fragment Visualization"], 
        icons=["house", "scatter-plot", "graph-up"], 
        menu_icon="cast", 
        default_index=0)

# Navigation logic based on selected menu option and imports necessary page when selected
if selected == "Home":
    from pages.Fragment_visualisation_pages import home
    home.show()
elif selected == "2D Fragment Visualization":
    from pages.Fragment_visualisation_pages import fragment_visualisation_2d as fragment_visualisation_2d
    fragment_visualisation_2d.show()
elif selected == "3D Fragment Visualization":
    from pages.Fragment_visualisation_pages import fragment_visualisation_3d as fragment_visualisation_3d
    fragment_visualisation_3d.show()

