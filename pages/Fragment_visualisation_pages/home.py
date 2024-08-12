## FRAGMENT VISUALISATION HOME ##
"""
This script displays necessary information and instructions to help user utilise the fragment visualisation 
menu and features. 
"""
import streamlit as st 

# Displays introductory text and instructions for the fragment visualisation features within the app
def show(): 
    st.title("Fragment Visualisation")
    st.write("""
             
             Explore the fragments through visualisation. 
            
              Use the sidebar to navigate between visualisations: 
                - **2D Fragment Visualisation**: View the fragmentation pattern in a 2D Scatter plot.
                - **Peptide and Protein Visualisation**: View Bradykinin in peptide structure and the full protein structure of Ubiquitin.
                
                ### How to Use:
                1. Select the type of visualisation you would like to see from the sidebar, i.e. 2D or 3D. 
                2. Select the peptide you would like to view.
                3. Interact with the plots to explore the fragments, peptides and final protein structure.
                4. Use the available controls to customise visualisation. 
                
              For the best view of each plot, expand to full screen. 
                """)