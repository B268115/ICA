# Programme to analyze protein from NCBI
# Written by s2762031
# Version 5, 12 Dec 2024

#---0. Preparation------------------------------------------------------------------------------------------------------------#
# Import the required modules
import os, subprocess
import re
import sys
import shutil

# Make an Output folder and the directory
# Our_Output = os.mkdir("./Our_Output")
our_output_directory = "./Our_Output"

# Welcome message
print("-----------------------WELCOME TO TINY PROTEIN ANALYZER----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")
print("glucose-6-phosphatase Aves")

#---1. Search protein from NCBI------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
## Add input from user: the name of the protein, and organism
#print("please write down your answer in lowercase and put '-' instead of space")
name_protein = input("What is the name of the protein?\n").strip()
name_organism = input("What is the name of the Organism? (aves or human)\n").strip()
### Error handling: spaces. Using 'raise ValueError to stop the programm'
Check01 = re.findall(r" ",name_protein)
Check02 = re.findall(r" ",name_organism)
if len(Check01) !=0  or len(Check02) !=0:
   raise ValueError("Your input are invalid, there are spaces, please change into '-' symbol")
### Error handling: empty or return
if not name_protein or not name_organism:
    raise ValueError("Input cannot be empty. Please provide a valid answer.")
### Error handling: input are too short
if len(name_protein) <3 or len(name_organism) <3:
    raise ValueError("Your input are too short. Please provie a real answer")
### Error handling: input in all capslock
if name_protein.isupper() or name_organism.isupper():
    raise ValueError("Your input is in caps lock, please only capitalize the first letter")

## 1.A Run the searching command
request_NCBI= f'esearch -db protein -query "{name_protein}[Protein] AND {name_organism}[Organism]" | efetch -format gb'
run01 = subprocess.run(request_NCBI,shell=True, capture_output=True,text=True)
### Save the output in variable
run01_output = run01.stdout
### Notify wheter the run01 command run successfully or not
print(f'Searching protein database in NCBI for {name_protein} in {name_organism}')
if run01.returncode == 0:
   print("Successfully searching your request in NCBI")
else:
   print("Something wrong")
   print("Error message:\n", process.stderr)

## 1.B Notify how many sequences you get
run02 = re.findall(r".*LOCUS*",run01_output)
print(f'Number of sequences: {len(run02)}')
### Error Handling: Notify if the number of sequences too big to be handled in this script.
if len(run02) >1000:
   raise ValueError("I am sorry, the number of sequences are too big too handle. Our limit is 1000 sequences. Please change the name of protein or the organism")

## 1.C Notify how many Organism you get
run03 = re.findall(r"ORGANISM\s+(.+)", run01_output, re.MULTILINE)
### Count only different species
unique_species =set(run03)
name_unique_species = sorted(unique_species)
print(f'Your sequences are coming from {len(unique_species)} different species')
#print(f'This is the name of the species:\n{name_unique_species}\n')

## 1.D Notify length of the organism
### get the length of all sequences, by retrieving the digit number in feature>source part
run04 = re.findall(r"\s+source\s+\d+\.\.(\d+)",run01_output,re.MULTILINE)
#print(run04)
#### chrosscheck num of length. It should be the same with the results of run02
#print(len(run04))
#### Loop through list run04 and average them
total_length = 0
for number in run04:
   total_length += int(number)
average_length = round(total_length/len(unique_species),2)
print(f'Average length of your sequences:{average_length}')
#### Error Handling: Notify if the length of aa too long to be handled in this script.
if average_length >= 1000:
   raise ValueError("I am sorry, the average num of your aa are too big for us to handle. Our limit is 1000 aa. Please search different protein name or organism.")

#---2. Get the Protein Sequences-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
## 2.A define a function to download .fasta
def download_protein_sequences(name_protein, name_organism, our_output_directory):
    fasta_file_name = "protein_sequence.fasta"
    fasta_file_path = os.path.join(our_output_directory, fasta_file_name)
    request_NCBI_2 = f'esearch -db protein -query "{name_protein}[Protein] AND {name_organism}[Organism]" | efetch -format fasta'

    with open(fasta_file_path, "w") as file:
        subprocess.run(request_NCBI_2, shell=True, text=True, stdout=file)

    print(f"We successfully downloaded your protein fasta. Check the Output folder: {fasta_file_path}\n\n")

## 2.B Ask the user again, whether they want to keep continue or not
if len(unique_species)==0:
    print("No data in NCBI database. Please input different name of organism or protein. Thank you")
    sys.exit()
elif len(unique_species) > 1:
    while True:  # Loop for error trapping
        confirmation = input("Your data seems to have more than one species. Would you like to continue? (yes/no):\n").lower()
        if confirmation == "yes":
            print("You choose 'yes'. We will download your protein sequences.")
            download_protein_sequences(name_protein, name_organism, our_output_directory)
            break # moving on from the while loop
        elif confirmation == "no":
            print("You choose 'no'. Thank you for using the Tiny Protein Analyzer. Have a nice day!")
            sys.exit()
        else:
            print("Invalid input. Please type 'yes' or 'no'.")
else:
   print("Your protein seq is only from one organism. We will download your protein sequences")
   download_protein_sequences(name_protein, name_organism, our_output_directory)


#---3. Filter Sequences Prior to Analyze the Level of Conservation----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
print("---------------Pre-Process: Filter Out non complete Protein Seq---------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

## 3.A Path to the pullseq .exe
pullseq_path = "/localdisk/data/BPSM/OptionalPythonICA/pullseq"

## 3.B Run pullseq command
input_pullseq = "Our_Output/protein_sequence.fasta"
output_pullseq = "Our_Output/filtered_protein_sequences.fasta"
### before filtering the sequences, define the max and min length of sequences
max_length = round(average_length + 200)
min_length = round(average_length - 100)
print(f'The max length of your input: {max_length} aa')
print(f'The min length of your input: {min_length} aa')
### get the header line of excluded sequences
header_excluded = "Our_Output/header_excluded_pullseq.txt"
with open(input_pullseq, "r") as input, open(header_excluded,"w") as output:
   for line in input:
      if line.startswith(">"):
         header = line[1:].strip()
         if "partial" in header:
            output.write(header + "\n")
### define the command for pullseq. In this script we only choose complete seq, length 300-500
pullseq_command = (f'{pullseq_path} -i {input_pullseq} -v -m {min_length} -a {max_length} -n {header_excluded} -e > {output_pullseq}')
#run05 = subprocess.run(pullseq_command,shell=True, check= True)
try:
    run05 = subprocess.run(pullseq_command, shell=True, check=True)
    print(f"Pullseq completed successfully! Output saved to {output_pullseq}\n")
except subprocess.CalledProcessError as e:
    print(f"Error running pullseq: {e}")

## 3.C Rename header of filtered fasta and make the sequence of each sequence as one long string
### currently the header of each sequences contain space. We will replace the space with _ and only grab the ID code and species name inside the brackets.
### define the input and output
input_reformat = output_pullseq
output_reformat = "Our_Output/filtered_protein_sequences_reformat.fasta"
### read the content of input files, output from this code is list
# Read the input file and open the output file
with open(input_reformat, 'r') as input_file, open(output_reformat, 'w') as output_file:
    # read the content, turn it into list
    input_file_content = input_file.readlines()
    sequence = ""  # Make an  empty string to store the sequence
    # loop through the list content
    for line in input_file_content:
        line = line.strip()  # Remove whitespace or newline in each element list
        if line.startswith(">"):
            # If sequence varible not empty, write it before the new header
            if sequence:
                output_file.write(sequence + '\n')  # Write the accumulated sequence with a newline
            # Rename the header
            ## find the id in header
            header_id = line.split(" ", 1)
            ID = header_id[0]
            ## find the organism name in header (inside the bracket pattern). We will find the position of the bracket
            header_organism_start = line.find("[") + 1
            header_organism_stop = line.find("]")
            organism = line[header_organism_start:header_organism_stop]
            ## assembly the new header
            new_header = f'{ID}_{organism}\n'.replace(" ", "_")
            ## write the new header to our output file
            output_file.write(new_header)
            # Reset the sequence for the next protein
            sequence = ""
        else:
            # To write the last sequence
            sequence += line  # The sequence lines are concatenated together
    # At the end, write the last sequence (if any)
    if sequence:
        output_file.write(sequence + '\n')
print(f'The header in fasta file: {input_reformat} has been updated to: {output_reformat}\n\n')


#---4. Multiple Sequence Alignment: Clustal-----------------------------------------------------------------------------------------------------------------------------------#
print("---------------Process_01: Multiple Sequence Allignment using Clustal-----------------------------")
## 4.A Path to clustalo .exe
### since the .exe alreade in /usr/bin, we dont have to specify it again.

## 4.B Run the clustalo command
### define the input and output location
input_clustalo = output_reformat
print(f'this is the input for clustalo:\n{input_clustalo}')
output_clustalo = "Our_Output/filtered_protein_sequences_aligned.msf"
### define the command for clustalo
clustalo_command = (f'clustalo -i {input_clustalo} -t protein --force --threads=4 --outfmt=msf -o {output_clustalo}')
#run06 = subprocess.run(clustalo_command,shell=True,check= True)
try:
    run06 = subprocess.run(clustalo_command, shell=True, check=True)
    print(f"Clustalo completed successfully! Output saved to {output_clustalo}\n\n")
except subprocess.CalledProcessError as e:
    print(f"Error running clustalo: {e}")

## 4.C Run Showalign (Emboss) to preview the results
### define the input and output
input_showalign = output_clustalo
print(f'this is the input for showing the result of clustalo:\n{input_showalign}')
output_showalign = "Our_Output/filtered_protein_sequences_aligned_pretty_format.txt"
### define the command for showalign
### we will preview the disimilarities
showalign_command = (f'showalign -show=d {input_showalign} -outfile {output_showalign}')
#run07 = subprocess.run(showalign_command,shell=True,check= True)
try:
    run07 = subprocess.run(showalign_command, shell=True, check=True)
    run08 = subprocess.run(f'head -20  {output_showalign}', shell=True, check=True)
    print(f"Preview clustal omega result is completed. To see more, open: {output_showalign}\n\n")
except subprocess.CalledProcessError as e:
    print(f"Error running preview output clustalo: {e}")


#---5. Plot conservation of sequence alignment: Plotcon----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
print("---------------Process_02: Plot conserve area of the sequences  using Plotcon----------------------------------------------------------------------------------------------")
## 5.A Path to plotcon
### since the .exe alreade in /usr/bin, we dont have to specify it again.

## 5.B Run the plotcon
### define the input and output
input_plotcon = "Our_Output/filtered_protein_sequences_aligned.msf"
output_plotcon_directory = "Our_Output/Output_Plotcon"
#### To prevent upcoming error, when repeat runing nano. We need to make>
if os.path.exists(output_plotcon_directory):
   shutil.rmtree(output_plotcon_directory)
os.mkdir(output_plotcon_directory)
### define the command for plotcon
plotcon_command = (f'plotcon -sformat msf {input_plotcon} -winsize 4 -graph png -gdirectory {output_plotcon_directory} -verbose -stdout')
#run07 = subprocess.run(showalign_command,shell=True,check= True)
try:
    run09 = subprocess.run(plotcon_command, shell=True, check=True)
    run10 = subprocess.run(f'cd "Our_Output/Output_Plotcon" && display *png', shell=True, check=True)
    print(f"Plot the sequence conservation is completed. The output is saved in directory {output_plotcon_directory}/n/n")
except subprocess.CalledProcessError as e:
    print(f"Error running preview output plotcon: {e}")

#---6. Find Percent Identity among our sequences using BLASTP -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------#
print("---------------Process_03: Analysis Percent Identity using Blastp-------------------------------------------------------------------------------------")

## 6.A Make database Blastp
## define input and output blast
input_blast = output_reformat
output_blast = "Our_Output/output_blastp.txt"
output_blast_2 = "Our_Output/output_blastp_high_identity.txt"
## make database blast
database_blast = "our_database"
makeblastdb_command = (f'makeblastdb -dbtype prot -in {input_blast} -out {database_blast}')
try:
   run11 = subprocess.run(makeblastdb_command, shell=True,check=True)
   print(f'Successfully makedatabase for blast')
except subprocess.CalledProcessError as e:
   print(f'Error running makeblasdb: {e}')

## 6.B Run the Blastp
### we set the output as tabular
blastp_command = (f'blastp -query {input_blast} -db {database_blast} -out {output_blast} -outfmt 6')
try:
   run12 = subprocess.run(blastp_command, shell=True,check=True)
   # show results with high  %identity
   run13 = subprocess.run(f"awk '$3 >95 && $3 <100 {{ print $0 }}' {output_blast} > { output_blast_2}", shell=True,check=True)
   # show the result of blastp on the screen
   run14 = subprocess.run(f'head -20 {output_blast_2}', shell=True, check=True)
   print(f'Successfully run blastp, you can see complete result of blastp in {output_blast} and {output_blast_2}')
except subprocess.CalledProcessError as e:
   print(f'Error running blastp: {e}')


#---7. Find the interesting domain: Using Blast search entire motif in Prosite  ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print("---------------Process_04: Searching Motif in Protein Sequences using Blast against Prosite database---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------")

## 7.A Download prosite database
### Define the location for database
prosite_database = "Our_Output/prosite_database.dat"
### Define the command
prosite_command = (f'wget https://ftp.expasy.org/databases/prosite/prosite.dat -O {prosite_database}')
### Run the command
try:
   run15 = subprocess.run(prosite_command, shell=True, check=True)
except subprocess.CalledProcessError as e:
   print(f'Error downloading prosite database: {e}')

## 7.B Running Motif searching using Blast
## must parse the content of prosite database first, into dictionary
