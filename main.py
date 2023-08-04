# import necessary libraries
import streamlit as st
from Bio import Entrez, SeqIO
from collections import Counter
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import multinomial
import numpy as np
st.set_page_config(page_title="MOV mutations analysis", page_icon="🧬", layout="wide")

spectra = [
    {"name":"BA.1",
     "url": "https://raw.githubusercontent.com/theosanderson/molnupiravir/main/mutational_spectra/BA.1_SBS_spectrum_Ruis.csv"
    },
    {"name":"High G-to-A",
        "url": "https://raw.githubusercontent.com/theosanderson/molnupiravir/main/mutational_spectra/long_phylogenetic_branches/long_branch_spectrum_rescaled.csv"    
    },
]

# function to fetch the reference genome
def fetch_reference_genome(accession):
    Entrez.email = "theo@theo.io"  # always tell NCBI who you are
    handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)

def count_all_trinucleotide_contexts(genome):
    counts = Counter()
    for i in range(len(genome)-2):
        context = genome[i:i+3]
        counts[context] += 1
    return counts

# function to get trinucleotide context
def get_context(genome_seq, mutation):
    pos = int(mutation[1:-1]) - 1  # -1 as Python uses 0-based indexing
    context = genome_seq[pos-1:pos+2]  # get the base before and after
    return context

# function to normalize mutation data
def normalize_data(mutation_data):
    # dummy function, replace with actual normalization code
    return mutation_data

# function to download and parse spectra data
def get_spectra_data(spectra_urls):
    spectra_data = {}
    for spectrum in spectra_urls:
        df = pd.read_csv(spectrum["url"])
        # we get mutation type from the 3rd-5th letters
        df['type'] = df['Substitution'].apply(lambda x: x[2:5])
        # we get context from 1st, 3rd and 5th letters
        df['context'] = df['Substitution'].apply(lambda x: x[0] + x[2] + x[6])
        df['value'] = df['Number_of_mutations']
        spectra_data[spectrum['name']] = df

    return spectra_data

def get_mut_type(mut_string):
    return mut_string[0] + '>' + mut_string[-1]

def plot_comparison(spectrum1, spectrum2, name1, name2, type_of_interest):
    # we want to filter to the type of interest, then plot.
    spectrum1 = spectrum1[spectrum1['type'] == type_of_interest]
    spectrum2 = spectrum2[spectrum2['type'] == type_of_interest]
    starting_nuc = type_of_interest[0] 
    # if either spectrum is empty, return NaN
    if len(spectrum1) == 0 or len(spectrum2) == 0:
        return np.nan
    # filter so context[1] is the starting nucleotide
    spectrum1 = spectrum1[spectrum1['context'].apply(lambda x: x[1] == starting_nuc)]
    spectrum2 = spectrum2[spectrum2['context'].apply(lambda x: x[1] == starting_nuc)]
    # We should join and assume any missing values are 0, then plot a scatter plot of each "value" col
    joined = spectrum1.merge(spectrum2, on='context', how='outer', suffixes=('_' + name1, '_' + name2))
    joined = joined.fillna(0)
    #st.write(joined)
    # normalize both values to sum to 1
    joined['value_' + name1] = joined['value_' + name1]/joined['value_' + name1].sum()
    joined['value_' + name2] = joined['value_' + name2]/joined['value_' + name2].sum()
    #st.write(spectrum1)
    #st.write(spectrum2)
 
    
    plt.scatter(joined['value_' + name1], joined['value_' + name2])
    plt.xlabel(name1)
    plt.ylabel(name2)
    plt.title(type_of_interest)
    #st.pyplot()
       
            
#   calculate multinomial log likelihood
    counts = joined['counts'].values
  
    p = (joined['Number_of_mutations']/joined['Number_of_mutations'].sum() ).values
  
    log_likelihood =float(multinomial.logpmf(counts, n=np.sum(counts), p=p))

    #st.write('Log likelihood: ', log_likelihood)
                                        
    return log_likelihood




def compare_and_report(mutation_info, spectra_data, comparison_list, spectra_list, mutation_name):
    # Ensure both lists have the same length
    assert len(comparison_list) == len(spectra_list)
   
    log_liks = []

    for i in range(len(comparison_list)):
        comp = plot_comparison(mutation_info, spectra_data[spectra_list[i]], 'My muts', spectra_list[i], mutation_name)
        log_liks.append(comp)

    # calculate the log likelihood ratio
    log_lik_ratio = log_liks[0] - log_liks[1]
    st.write('Log likelihood increment towards MOV for ', mutation_name, ': ', np.round(log_lik_ratio, 2))

    return log_lik_ratio

# Streamlit app
def main():
    st.title("Analysis of branch mutation spectra for MOV signature")
    # window title
 
    
                        
    # window title
    

    # input mutation data
    mutation_data = st.text_input("Enter mutations (comma-separated, e.g. 'G123A,T5343A'): ", "C2595T, T3607C, C4464T, C6525T, G9092A, G11272A, G12067A, C12784T, G14430A, C14605T, G18589A, G18712A, G20839A, A22124G, C23170T, T25735C, A26169G, C28344T, C28697T, C29098T")
    mutation_data = mutation_data.replace("nt:","")
    mutation_data = mutation_data.replace(" ","")
    mutations = mutation_data.split(',')

    # fetch reference genome
    genome_seq = fetch_reference_genome('NC_045512.2')

    # get trinucleotide context and calculate normalized values
    mutation_info = {}
    for mutation in mutations:
        context = get_context(genome_seq, mutation)
        mutation_info[mutation] = {'context': context, 'type': get_mut_type(mutation)}
        
    # make a bar chart of mutation_types
    counts_types = Counter([mutation_info[mutation]['type'] for mutation in mutations])

    # hr 
    st.markdown("***")
    st.subheader("Illustration of mutation types")


    st.bar_chart(pd.DataFrame.from_dict(counts_types, orient='index'))
    st.markdown("***")
    st.subheader("Assessment of MOV-likeness of mutation context")
    st.markdown("*This section is wholly separate from the mutation classes, and considers only the contexts in which transition mutations occur.*")
    
    genome_counts = count_all_trinucleotide_contexts(genome_seq)
    # make df
    mutation_info = pd.DataFrame.from_dict(mutation_info, orient='index')
    
    # sum by context and type
    mutation_info = mutation_info.groupby(['context','type']).size().reset_index(name='counts')

    # normalize to genome_counts
    mutation_info['value'] = mutation_info.apply(lambda x: x['counts']/genome_counts[x['context']], axis=1)

   
    # get spectra data
    spectra_data = get_spectra_data(spectra)
    #st.write(spectra_data['BA.1'])
    #st.write(spectra_data['High G-to-A'])



    # Define the list of comparisons and spectra to use
    comparison_list = ['High G-to-A', 'BA.1']
    spectra_list = ['High G-to-A', 'BA.1']
    # count transition muts
    transition_muts = mutation_info[mutation_info['type'].isin(['G>A', 'C>T', 'A>G', 'T>C'])]
    st.write("You entered ", len(transition_muts), " transition mutations. ", "**This may be too few too have any confidence in the results.**" if len(transition_muts) < 5 else "")
    
    # Make comparisons
    ga = compare_and_report(mutation_info, spectra_data, comparison_list, spectra_list, 'G>A')
    ct = compare_and_report(mutation_info, spectra_data, comparison_list, spectra_list, 'C>T')
    ag = compare_and_report(mutation_info, spectra_data, comparison_list, spectra_list, 'A>G')
    tc = compare_and_report(mutation_info, spectra_data, comparison_list, spectra_list, 'T>C')
    # calculate mean weighted by number of mutations of each class in counts_types
    # account for the fact that any of these could be nan, in which case it is excluded
    
    summed = np.nansum([ga, ct, ag, tc])

    st.write('Summed log likelihood ratio: ', np.round(summed, 3))

    # calculate probability of MOV-like contexts
    prob = 1/(1+np.exp(-summed))
    # below should be bold
    st.write('**Estimated un-normalised probability of MOV-like contexts:** ', np.round(prob, 3))
    # horizontal rule
    st.markdown('<hr>', unsafe_allow_html=True)
    st.write("The probability above does not account for the fact that the prior for MOV is low, and therefore in some sense may be overestimated. On the other hand, the sequence only looks at contexts and not at mutation types, and so if the sequence was identified on the basis of a high G-to-A signature and so the prior is already high, the probability may not be underestimated.")

    # in red, add a disclaimer
    st.write("""Disclaimer: Do not rely primarily on this application for interpreting molnupiravir\'s role 
             in a spectific sequence. This application should only be used by experts who can interpret the
                results in the context of other evidence. """)
    
    st.markdown('<hr>', unsafe_allow_html=True)
    st.write(mutation_info)




  

if __name__ == '__main__':
    main()
