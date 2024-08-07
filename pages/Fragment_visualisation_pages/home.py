import streamlit as st 

def show(): 
    st.title("Fragment Visualisation")
    st.write("""
             
             Explore the fragments through visualisation. 
            
              Use the sidebar to navigate between visualisations: 
                - **2D Fragment Visualisation**: View the fragmentation pattern in a 2D Scatter plot.
                - **3D Fragment Visualisation**: View the fragmentation pattern in a 3D scatter plot. 
                
                ### How to Use:
                1. Select the type of visualisation you would like to see from the sidebar. 
                2. Interact with the plots to explore the fragments.
                3. Use the available controls to customise visualisation. 
                
                """)