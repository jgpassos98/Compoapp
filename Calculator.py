import streamlit as st
import itertools
import pandas as pd
import numpy as np

# Read the Excel file containing element-VEC data
df_VEC = pd.read_excel("VEC.xlsx")

# Read the Excel file containing element-Bulk Modulus data
bulk_modulus_values = pd.read_excel("bulk_modulus.xlsx")

# Read the Excel file containing element-Radius data
atomic_radius_values = pd.read_excel("radius.xlsx")

# Extract the "Element" and "VEC" columns into separate lists
elements_VEC = df_VEC["Element"].tolist()
vec_values = df_VEC["VEC"].tolist()

# Read the Excel file containing element-Tm data
df_Tm = pd.read_excel("Tm.xlsx")

# Extract the "Element" and "Tm (°C)" columns into separate lists
elements_Tm = df_Tm["Element"].tolist()
Tm_values = df_Tm["Tm (Celsius)"].tolist()

# Read the Excel file containing element-w0el data
df_w0el = pd.read_excel("w0el.xlsx", index_col=0)

# Read the Excel file containing element-VEC data
df_w1el = pd.read_excel("w1el.xlsx", index_col=0)

# Read the Excel file containing element-VEC data
df_w2el = pd.read_excel("w2el.xlsx", index_col=0)

# Read the Excel file containing element-VEC data
df_w3el = pd.read_excel("w3el.xlsx", index_col=0)

def calculate_atomic_fraction_product(df_filtered_DHel2):
    # Extract elements and compositions from the DataFrame
    elements = df_filtered_DHel2.columns.tolist()

    # Initialize a list to store the atomic fraction product matrices for each alloy
    atomic_fraction_product_matrices = []

    for _, row in df_filtered_DHel2.iterrows():
        compositions = row.to_dict()

        # Calculate the atomic fraction product matrix for the current alloy
        n = len(elements)
        atomic_fraction_product_matrix = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                if i != j:
                    atomic_fraction_product_matrix[i][j] = compositions[elements[i]] * compositions[elements[j]] / 10000

        atomic_fraction_product_matrices.append(atomic_fraction_product_matrix)

    return atomic_fraction_product_matrices
# Now, 'transpose_matrixes' contains the transpose of each composition matrix
def calculate_relative_composition(df_filtered_DHel2):
    # Extract elements and compositions from the DataFrame

    elements = df_filtered_DHel2.columns.tolist()

    # Initialize a list to store relative composition matrices for each alloy
    relative_composition_matrices = []

    for _, row in df_filtered_DHel2.iterrows():
        compositions = row.to_dict()

        # Calculate the relative composition matrix for the current alloy
        n = len(elements)
        relative_composition_matrix = [[0.0] * n for _ in range(n)]

        for i in range(n):
            for j in range(n):
                if i != j:
                    relative_composition_matrix[i][j] = compositions[elements[i]] / (compositions[elements[i]] + compositions[elements[j]])

        relative_composition_matrices.append(relative_composition_matrix)

    return relative_composition_matrices

def calculate_mixing_enthalpy(df_filtered_DHel2, df_w0el, df_w1el, df_w2el, df_w3el, relative_composition_matrices, atomic_fraction_product_matrices):
    # Initialize a list to store the mixing enthalpy matrices for each alloy
    total_enthalpies = []

    # Iterate over each alloy (each row in df_elements)
    for idx, (index, alloy_row) in enumerate(df_filtered_DHel2.iterrows()):
        # Extract the elements and their compositions for the current alloy
        elements = alloy_row.index.tolist()
        compositions = alloy_row.tolist()

        # Initialize a matrix to store mixing enthalpy for the current alloy
        n = len(elements)
        mixing_enthalpy_matrix = np.zeros((n, n))

        # Iterate over all possible pairs of elements
        total_enthalpy = 0.0
        for i in range(n):
            for j in range(n):
                if i != j:
                    # Get the names of the current element pair
                    element_i = elements[i]
                    element_j = elements[j]

                    # Fetch the values for the current element pair from the element matrices
                    el_values = {
                        'w0el': df_w0el.loc[element_i, element_j],
                        'w1el': df_w1el.loc[element_i, element_j],
                        'w2el': df_w2el.loc[element_i, element_j],
                        'w3el': df_w3el.loc[element_i, element_j]
                    }

                    # Calculate the mixing enthalpy using the provided formula
                     # Use the relative composition matrix for the current alloy
                    relative_composition_matrix = relative_composition_matrices[idx]
                    
                     # Get the atomic fraction product for the current element pair
                    atomic_fraction_product = atomic_fraction_product_matrices[idx][i][j]

                    delta_C = relative_composition_matrix[i][j] - relative_composition_matrix[j][i]
                    H_mix = 4 * (el_values['w0el'] + el_values['w1el'] * delta_C + el_values['w2el'] * delta_C ** 2 + el_values['w3el'] * delta_C ** 3) * atomic_fraction_product

                    # Replace NaN values with 0
                    H_mix = 0 if np.isnan(H_mix) else H_mix

                    # Add the enthalpy to the total enthalpy for the alloy
                    total_enthalpy += H_mix
        # Append the total enthalpy for the current alloy to the list
        total_enthalpies.append(total_enthalpy)
        
    # After the outer loop, outside both loops, assign DHel_values to the DataFrame column
    df_filtered_DHel2['Enthalpy of mixing (DHmix)'] = total_enthalpies

    # Finally, return the DataFrame
    return df_filtered_DHel2

def calculate_elastic_strain_energy(df_filtered_Tm2, bulk_modulus_values, atomic_radius_values):
    num_alloys = len(df_filtered_Tm2)
    DHel_values = []

    # Step 1: Calculate ci.Bi for each element
    ci_Bi_values = []
    for index, row in df_filtered_Tm2.iterrows():
        alloy_composition = row / 100  # Convert compositions to fractions
        ci_Bi = []

        for element, composition in alloy_composition.items():
            if composition > 0:  # Exclude elements with zero composition
                # Find the index of the element in the bulk modulus values dataframe
                element_index = bulk_modulus_values.index[bulk_modulus_values['Element'] == element].tolist()[0]

                # Retrieve the bulk modulus value for the element
                Bi = bulk_modulus_values.loc[element_index, 'Bulk Modulus']

                # Calculate ci.Bi for the element
                ci_Bi_element = Bi * composition * 1e9
                ci_Bi.append(ci_Bi_element)

        ci_Bi_values.append(ci_Bi)

    # Initialize vi_values vector
    vi_values = []

    # Step 3: Calculate ci.Bi.Vi for each element
    ci_Bi_Vi_values = []
    for index, row in df_filtered_Tm2.iterrows():
        alloy_composition = row / 100  # Convert compositions to fractions
        ci_Bi_Vi = []
        vi = []  # Initialize vi for the alloy composition
        for element, composition in alloy_composition.items():
            if composition > 0:  # Exclude elements with zero composition
                # Find the index of the element in the atomic radius values dataframe
                element_index = atomic_radius_values.index[atomic_radius_values['Element'] == element].tolist()[0]

                # Retrieve the radius value for the element
                ri = atomic_radius_values.loc[element_index, 'radius']

                # Calculate the corresponding Vi value using the retrieved radius
                vi_element = (4 / 3) * np.pi * (ri ** 3) * (10 ** -27) * 6.02e23
                
                # Append the Vi value for the alloy composition
                vi.append(vi_element)
                
        vi_values.append(vi)

        # Calculate ci.Bi.Vi for the alloy composition
        for ci_Bi_element, vi_element in zip(ci_Bi_values[index], vi_values[index]):  # Using index to access ci_Bi for the current alloy composition
            ci_Bi_Vi_element = ci_Bi_element * vi_element  # Calculate ci.Bi.Vi for the element
            ci_Bi_Vi.append(ci_Bi_Vi_element)

        # Append ci_Bi_Vi for the current alloy composition
        ci_Bi_Vi_values.append(ci_Bi_Vi)

    V_values = []
    for ci_Bi, ci_Bi_Vi in zip(ci_Bi_values, ci_Bi_Vi_values):
        sum_ci_Bi = np.sum(ci_Bi)
        sum_ci_Bi_Vi = np.sum(ci_Bi_Vi)

        # Avoid division by zero
        if sum_ci_Bi != 0:
            V_value = sum_ci_Bi_Vi / sum_ci_Bi
        else:
            V_value = 0

        V_values.append(V_value)

    # Now, you can use vi_values in the calculation of DHel_element
    
    DHel_values = []
    for ci_Bi, V, vi in zip(ci_Bi_values, V_values, vi_values):
        # Initialize total_DHel as a scalar
        total_DHel = 0
        for ci_Bi_element, vi_element in zip(ci_Bi, vi):
            DHel_element = ci_Bi_element * (((vi_element - V) ** 2) / (2 * vi_element))
            total_DHel += DHel_element  # Accumulate the DHel_element into total_DHel

        # Append the total_DHel divided by 1000 after the inner loop
        DHel_values.append(total_DHel / 1000)
        

    # After the outer loop, outside both loops, assign DHel_values to the DataFrame column
    df_filtered_Tm2['Elastic Strain Energy (DHel)'] = DHel_values

    # Finally, return the DataFrame
    return df_filtered_Tm2

    
def calculate_configurational_entropy(df):
    num_alloys = len(df)
    entropies = []

    for index, row in df.iterrows():
        mole_fractions = row / 100  # Convert compositions to mole fractions
        mole_fractions = mole_fractions[mole_fractions > 0]  # Exclude elements with zero composition
        entropy = np.sum(mole_fractions * np.log(mole_fractions))
        entropies.append(entropy)

    df['Configurational Entropy (J/mol·K)'] = entropies
    return df
    
    #df is the original df
    
def calculate_average_vec(df_filtered2, df_VEC):
    num_alloys = len(df_filtered2)
    average_vecs = []

    for index, row in df_filtered2.iterrows():
        alloy_composition = row / 100  # Convert compositions to fractions
        alloy_vec = 0
        total_composition = 0

        for element, composition in alloy_composition.items():
            if composition > 0:  # Exclude elements with zero composition
                element_vec = df_VEC.loc[df_VEC['Element'] == element, 'VEC'].iloc[0]
                alloy_vec += element_vec * composition
                total_composition += composition

        average_vec = alloy_vec / total_composition if total_composition > 0 else np.nan
        average_vecs.append(average_vec)

    df_filtered2['Average VEC'] = average_vecs
    return df_filtered2
    
    #df_filtered2 is the original df filtered for configurational entropy without the entropy tab
    
def calculate_weighted_average_Tm(df_filtered_VEC2, df_Tm):
    num_alloys = len(df_filtered_VEC2)
    average_Tms = []

    for index, row in df_filtered_VEC2.iterrows():
        alloy_composition = row / 100  # Convert compositions to fractions
        alloy_Tm = 0
        total_composition = 0

        for element, composition in alloy_composition.items():
            if composition > 0:  # Exclude elements with zero composition
                # Find the melting temperature for the element from df_Tm
                element_Tm = df_Tm.loc[df_Tm['Element'] == element, 'Tm (Celsius)'].iloc[0]
                alloy_Tm += element_Tm * composition
                total_composition += composition

        average_Tm = alloy_Tm / total_composition if total_composition > 0 else np.nan
        average_Tms.append(average_Tm)

    df_filtered_VEC2['Average Tm'] = average_Tms
    return df_filtered_VEC2
    
    #Same structure with the df filtered for VEC without VEC column

def generate_compositions_recursive(num_elements, elements, lower_limits, upper_limits, step, composition, index, compositions):
    # If all elements except the last one have been assigned a value, add the composition to the list
    if index == num_elements - 1:
        # Calculate the remaining composition for the last element
        remaining_composition = 100 - sum(composition)
        # Check if the remaining composition is within the range for the last element
        if lower_limits[num_elements - 1] <= remaining_composition <= upper_limits[num_elements - 1]:
            # Append the composition to the list of compositions
            compositions.append(composition + [remaining_composition])
        return

    # Generate compositions for the current element within the specified range
    for value in np.arange(lower_limits[index], upper_limits[index] + step, step):
        # Assign the value to the current element in the composition
        composition[index] = value
        # Recur to the next element
        generate_compositions_recursive(num_elements, elements, lower_limits, upper_limits, step, composition, index + 1, compositions)


def generate_alloy_compositions(num_elements, elements, lower_limits, upper_limits, step):
    compositions = []
    # Initialize an empty composition array
    composition = [0] * (num_elements - 1)
    # Generate compositions recursively
    generate_compositions_recursive(num_elements, elements, lower_limits, upper_limits, step, composition, 0, compositions)
    # Create a DataFrame from the compositions
    df = pd.DataFrame(compositions, columns=elements)
    return df
    
def df_to_excel(df):
    # Convert DataFrame to Excel file content in memory
    excel_buffer = pd.ExcelWriter(path=None, engine="openpyxl")
    df.to_excel(excel_buffer, index=False, header=True)
    excel_buffer.seek(0)
    excel_content = excel_buffer.read()
    excel_buffer.close()
    return excel_content
        
def main():

    st.title("Alloy Composition App")

    # Create tabs for different functionalities
    tabs = ["Generate Compositions", "Filter Uploaded Compositions"]
    selected_tab = st.sidebar.selectbox("Select Tab", tabs)

    if selected_tab == "Generate Compositions":
        generate_compositions_tab()
    elif selected_tab == "Filter Uploaded Compositions":
        view_uploaded_compositions_tab()
        
def generate_compositions_tab():
        st.title("Generate Compositions")
    
        num_elements = st.number_input("Enter the number of elements in the alloy:", min_value=1, step=1)
        elements = [st.text_input(f"Element {i+1}:") for i in range(num_elements)]
        
        st.write("Enter lower and upper limits for each element:")
        limits_col1, limits_col2 = st.columns(2)
        lower_limits = []
        upper_limits = []
        for i in range(num_elements):
            with limits_col1:
                lower_limits.append(st.number_input(f"Lower limit for {elements[i]}:", step=0.01))
            with limits_col2:
                upper_limits.append(st.number_input(f"Upper limit for {elements[i]}:", step=0.01))
        
        step = st.number_input("Step size:", step=0.01)

        if st.button("Generate Alloy Compositions"):
            df = generate_alloy_compositions(num_elements, elements, lower_limits, upper_limits, step)
            st.write(df)
        
        if st.button("Save to Excel"):
            df = generate_alloy_compositions(num_elements, elements, lower_limits, upper_limits, step)
            excel_file = df_to_excel(df)
            st.download_button(label="Download Excel", data=excel_file, file_name="alloy_compositions.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", key=None)
            st.success("Click the button above to download the Excel file.")
            

def view_uploaded_compositions_tab():
    st.title("Apply filters (follow the steps every time a value is changed)")

    # Upload Excel file
    uploaded_file = st.file_uploader("Upload Excel file", type=["xlsx", "xls"])

    # Initialize session state if it doesn't exist
    if "df" not in st.session_state:  # Added this line
        st.session_state.df = None  # Added this line

    if uploaded_file is not None:
        # Read Excel file as DataFrame
        if st.session_state.df is None:  # Added this line
            st.session_state.df = pd.read_excel(uploaded_file)  # Changed this line
        else:
            df = pd.read_excel(uploaded_file)
            st.session_state.df = df  # Changed this line

        # Input fields for configurational entropy threshold
        min_entropy, max_entropy = st.columns(2)
        min_entropy_value = min_entropy.number_input("Minimum Configurational Entropy (J/mol·K):", key="min_entropy", step=0.1)
        max_entropy_value = max_entropy.number_input("Maximum Configurational Entropy (J/mol·K):", key="max_entropy", step=0.1)
        
        # Input fields for VEC threshold
        min_VEC, max_VEC = st.columns(2)
        min_vec_value = min_VEC.number_input("Minimum VEC:", key="min_VEC", step=0.1)
        max_vec_value = max_VEC.number_input("Maximum VEC:", key="max_VEC", step=0.1)
        
        # Input fields for Tm threshold
        min_Tm, max_Tm = st.columns(2)
        min_Tm_value = min_Tm.number_input("Minimum Tm (°C):", key="min_Tm", step=0.1)
        max_Tm_value = max_Tm.number_input("Maximum Tm (°C):", key="max_Tm", step=0.1)
        
        # Input fields for DHel threshold
        min_DHel, max_DHel = st.columns(2)
        min_DHel_value = min_DHel.number_input("Minimum DHel (kJ/mol):", key="min_DHel", step=0.1)
        max_DHel_value = max_DHel.number_input("Maximum DHel (kJ/mol):", key="max_DHel", step=0.1)
        
        # Input fields for DHmix threshold
        min_DHmix, max_DHmix = st.columns(2)
        min_DHmix_value = min_DHmix.number_input("Minimum DHmix (kJ/mol):", key="min_DHmix", step=0.1)
        max_DHmix_value = max_DHmix.number_input("Maximum DHmix (kJ/mol):", key="max_DHmix", step=0.1)

    filter_button_col, _ = st.columns([3, 1])
    df = st.session_state.df
    if filter_button_col.button("Step 1: Apply Configurational entropy filter"):
            entropies = calculate_configurational_entropy(df)
            entropies_array = np.array(entropies)
            
            # Filter the dataframe based on configurational entropy threshold
            df_filtered = df[(entropies_array >= min_entropy_value) & (entropies_array <= max_entropy_value)]
            
            # Drop the configurational entropy column
            df_filtered2 = df_filtered.drop(columns=["Configurational Entropy (J/mol·K)"], errors="ignore")
            
            # Overwrite the original DataFrame with the filtered DataFrame
            st.session_state.df = df_filtered2
            
            st.session_state.df_filtered2 = df_filtered2
            
            # Save the filtered dataframe back to the original Excel file
            df_filtered2.to_excel(uploaded_file.name, index=False, engine='openpyxl')
            st.success("Filtered data saved to the uploaded file.")
                
    # Display button for VEC filter
    if st.button("Step 2: Apply VEC Filter"):
        if st.session_state.df_filtered2 is not None:

            # Calculate average VEC for the filtered compositions
            df_filtered2 = calculate_average_vec(st.session_state.df_filtered2, df_VEC)
            st.session_state.df_filtered2 = df_filtered2

            # Filter based on VEC values
            df_filtered_VEC = st.session_state.df_filtered2[(st.session_state.df_filtered2['Average VEC'] >= min_vec_value) & 
                                                            (st.session_state.df_filtered2['Average VEC'] <= max_vec_value)]

            # Drop the 'Average VEC' column
            df_filtered_VEC2 = df_filtered_VEC.drop(columns=["Average VEC"], errors="ignore")
            
            st.session_state.df_filtered_VEC2 = df_filtered_VEC2

            # Save the filtered dataframe back to the original Excel file
            df_filtered_VEC2.to_excel(uploaded_file.name, index=False, engine='openpyxl')
            st.success("Filtered data saved to the uploaded file.")
            
    # Display button for Tm filter
    if st.button("Step 3: Apply Tm Filter"):
        if st.session_state.df_filtered_VEC2 is not None:

            # Calculate average Tm for the filtered compositions
            df_filtered_VEC2 = calculate_weighted_average_Tm(st.session_state.df_filtered_VEC2, df_Tm)
            st.session_state.df_filtered_VEC2 = df_filtered_VEC2

            # Filter based on VEC values
            df_filtered_Tm = st.session_state.df_filtered_VEC2[(st.session_state.df_filtered_VEC2['Average Tm'] >= min_Tm_value) & 
                                                               (st.session_state.df_filtered_VEC2['Average Tm'] <= max_Tm_value)]

            # Drop the 'Average VEC' column
            df_filtered_Tm2 = df_filtered_Tm.drop(columns=["Average Tm"], errors="ignore")

            st.session_state.df_filtered_Tm2 = df_filtered_Tm2

            # Save the filtered dataframe back to the original Excel file
            df_filtered_Tm2.to_excel(uploaded_file.name, index=False, engine='openpyxl')
            st.success("Filtered data saved to the uploaded file.")
            
    # Display button for DHel filter
    if st.button("Step 3: Apply DHel Filter"):
        if st.session_state.df_filtered_Tm2 is not None:

            # Calculate average DHel for the filtered compositions
            df_filtered_Tm2 = calculate_elastic_strain_energy(st.session_state.df_filtered_Tm2, bulk_modulus_values, atomic_radius_values)
            st.session_state.df_filtered_Tm2 = df_filtered_Tm2

            # Filter based on VEC values
            df_filtered_DHel = st.session_state.df_filtered_Tm2[(st.session_state.df_filtered_Tm2['Elastic Strain Energy (DHel)'] >= min_DHel_value) & 
                                                              (st.session_state.df_filtered_Tm2['Elastic Strain Energy (DHel)'] <= max_DHel_value)]

            # Drop the 'Average VEC' column
            df_filtered_DHel2 = df_filtered_DHel.drop(columns=["Elastic Strain Energy (DHel)"], errors="ignore")

            st.session_state.df_filtered_DHel2 = df_filtered_DHel2

            st.write("Filtered DataFrame with Tm Filter:")
            st.write(df_filtered_DHel2)

            # Save the filtered dataframe back to the original Excel file
            df_filtered_DHel2.to_excel(uploaded_file.name, index=False, engine='openpyxl')
            st.success("Filtered data saved to the uploaded file.")
            
        # Display button for DHmix filter
    if st.button("Step 4: Apply DHmix Filter"):
        if st.session_state.df_filtered_DHel2 is not None:

            # Calculate the relative composition matrices
            relative_composition_matrices = calculate_relative_composition(st.session_state.df_filtered_DHel2)
            
            # Calculate the atomic fraction product matrices
            atomic_fraction_product_matrices = calculate_atomic_fraction_product(st.session_state.df_filtered_DHel2)

            # Calculate average DHel for the filtered compositions
            df_filtered_DHel2 = calculate_mixing_enthalpy(st.session_state.df_filtered_DHel2, df_w0el, df_w1el, df_w2el, df_w3el, relative_composition_matrices, atomic_fraction_product_matrices)
            
            st.session_state.df_filtered_DHel2 = df_filtered_DHel2  

            # Print the DataFrame with average Tm to check its content
            st.write("DataFrame with average Tm:4")
            st.write(df_filtered_DHel2)

            # Filter based on DHmix values
            df_filtered_DHmix = st.session_state.df_filtered_DHel2[(st.session_state.df_filtered_DHel2['Enthalpy of mixing (DHmix)'] >= min_DHmix_value) & 
                                                              (st.session_state.df_filtered_DHel2['Enthalpy of mixing (DHmix)'] <= max_DHmix_value)]

            # Print the filtered DataFrame to check its content
            st.write("Filtered DataFrame with Tm Filter:")
            st.write(df_filtered_DHmix)

            # Drop the 'DHmix' column
            df_filtered_DHmix2 = df_filtered_DHmix.drop(columns=["Enthalpy of mixing (DHmix)"], errors="ignore")

            # Print the DataFrame without 'DHmix' column to check its content
            st.write("Filtered DataFrame with DHmix Filter and column removed:")
            st.write(df_filtered_DHmix2)
                        
            st.session_state.df_filtered_DHmix2 = df_filtered_DHmix2

            # Save the filtered dataframe back to the original Excel file
            df_filtered_DHmix2.to_excel(uploaded_file.name, index=False, engine='openpyxl')
            st.success("Filtered data saved to the uploaded file.")
                
if __name__ == "__main__":
    main()
