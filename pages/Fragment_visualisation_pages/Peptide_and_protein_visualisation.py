import streamlit as st
import streamlit.components.v1 as components 

# Function to generate the HTML for JSmol viewer
def generate_jsmol_html(pdb_id):
    return f"""
    <script type="text/javascript" src="https://chemapps.stolaf.edu/jmol/jsmol/JSmol.min.js"></script>
    <script type="text/javascript">
        Jmol.getApplet("jsmolApplet0", {{
            width: 500,
            height: 700,
            use: "HTML5",
            j2sPath: "https://chemapps.stolaf.edu/jmol/jsmol/j2s",
            script: "load https://files.rcsb.org/view/{pdb_id}.pdb; cartoons only; color structure;"
        }});
    </script>
    <div id="jsmolApplet0" style="margin:auto;"></div>
    """

# Streamlit app to visualize a protein structure in 3D using JSmol
def show():
    st.title("3D Protein Structure Visualization")
    st.markdown("Explore the 3D protein structure of Bradykinin as a peptide and the full protein structure of Ubiquitin.") 
    st.write("With structures directly taken from RCSB Protein Data Bank to provide a visual image of both a peptide structure and a protein structure.")

    peptide_options = {"Ubiquitin", "Bradykinin"}
    selected_peptide = st.sidebar.selectbox("Select Peptide", peptide_options)
                                        
    peptide_pdb = {
        "Ubiquitin": '1UBQ',
        "Bradykinin": '6F3V', 
        
    }
    # PDB ID for 6ZNS is hardcoded here
    pdb_id = peptide_pdb[selected_peptide]

    # Generate the JSmol HTML
    jsmol_html = generate_jsmol_html(pdb_id)

    # Display the HTML with JSmol embedded
    st.components.v1.html(jsmol_html, height=650, width=550)

# Example of how to use this in your Streamlit navigation
if __name__ == "__main__":
    show()


