# Programme to analyze protein from NCBI
# Written by s2762031
# Version 1, 22 Nov 2024

#---0. Preparation------------------------------------------------------------------------------------------------------------#
# Import the required modules
import os, subprocess
import re
import sys

# Make an Output folder and the directory
# Our_Output = os.mkdir("./Our_Output")
our_output_directory = "./Our_Output"

# Welcome message
print("-----------------------WELCOME TO TINY PROTEIN ANALYZER------------------------------------------")
print("glucose-6-phosphatase Aves")

#---1. Search protein from NCBI--------------------------------------------------------------------------------------------------#
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

#---2. Get the Protein Sequences-------------------------------------------------------------------------------------------------------------------#
## 2.A define a function to download .fasta
def download_protein_sequences(name_protein, name_organism, our_output_directory):
    fasta_file_name = "protein_sequence.fasta"
    fasta_file_path = os.path.join(our_output_directory, fasta_file_name)
    request_NCBI_2 = f'esearch -db protein -query "{name_protein}[Protein] AND {name_organism}[Organism]" | efetch -format fasta'

    with open(fasta_file_path, "w") as file:
        subprocess.run(request_NCBI_2, shell=True, text=True, stdout=file)

    print(f"We successfully downloaded your protein fasta. Check the Output folder: {fasta_file_path}")

## 2.B Ask the user again, whether they want to keep continue or not
if len(unique_species)==0:
    print("No data in NCBI database. Please input different name of organism or protein. Thank you")
    sys.exit()
elif len(unique_species) > 1:
    while True:  # Loop for error trapping
        confirmation = input("Your data seems to have more than one species. Would you like to continue? (yes/no): ").lower()
        if confirmation == "yes":
            print("You chose 'yes'. We will download your protein sequences.")
            download_protein_sequences(name_protein, name_organism, our_output_directory)
            break # moving on from the while loop
        elif confirmation == "no":
            print("You chose 'no'. Thank you for using the Tiny Protein Analyzer. Have a nice day!")
            sys.exit()
        else:
            print("Invalid input. Please type 'yes' or 'no'.")
else:
   print("Your protein seq is only from one organism. We will download your protein sequences")
   download_protein_sequences(name_protein, name_organism, our_output_directory)


#---3. Filter Sequences Prior to Analyze the Level of Conservation-------------------------------------------------------------------------------------------------------------------------#
print(f'Starting to analyze conservation level between protein sequences')
