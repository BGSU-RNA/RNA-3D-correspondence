import csv
from datetime import datetime
import gzip
import os
import re
import requests
import sys

from RPFAM_configuration import BASE_DIRECTORY
from RPFAM_configuration import LOOP_DIRECTORY
from RPFAM_configuration import RFAM_ALIGNMENT_DIRECTORY
from RPFAM_configuration import PFAM_ALIGNMENT_DIRECTORY
from RPFAM_configuration import OUTPUT_DIRECTORY

get_organisms = False
if get_organisms:
    from bs4 import BeautifulSoup

# manage imports according to Python version present
if sys.version_info[0] < 3:
    # for experimental sequence mapping to resolved sequences
    import pdbx    # currently only really working for Python 2.7
    from pdbx.reader.PdbxReader import PdbxReader
    from urllib import urlretrieve as urlretrieve
    read_mode = 'rb'
    write_mode = 'w'
else:
    from urllib.request import urlretrieve as urlretrieve
    read_mode = 'rt'
    write_mode = 'wt'   # write as text

# from pfam_chain_mapping import Pfam_family_mapping

ANCHOR_PDB_STRUCTURES = []
PRESENCE_OF_PROTEINS = False
CIF_FILE_LOCATION = "cif_files/"

# lists of Rfam families to include in cross-domain alignments
# align according to the covariance model for the first family
joint_alignments = {}
joint_alignments['SSU'] = ['RF01960','RF01959','RF00177','RF02545','RF02542']
joint_alignments['LSU'] = ['RF02543','RF02540','RF02541','RF02546']



class PhasedMatch:
    def __init__(self, ncbi_taxid=None, rfam_seq_id=None, pfam_id=None, pfam_sequence=None, rfam_sequence=None, organism=None, fasta_title_remainder_rna=None, fasta_title_remainder_protein=None, rfam_aligned_sequence=None, pfam_aligned_sequence=None, phased_id_RNA=None, phased_id_protein=None):
        self.ncbi_taxid = ncbi_taxid
        self.rfam_seq_id = rfam_seq_id
        self.pfam_id = pfam_id
        self.pfam_sequence = pfam_sequence
        self.rfam_sequence = rfam_sequence
        self.organism = organism
        self.fasta_title_remainder_rna = fasta_title_remainder_rna
        self.fasta_title_remainder_protein = fasta_title_remainder_protein
        self.rfam_aligned_sequence = rfam_aligned_sequence
        self.pfam_aligned_sequence = pfam_aligned_sequence
        self.phased_id_RNA = phased_id_RNA
        self.phased_id_protein = phased_id_protein

        # This establishes a class allowing for all alignment information to be stored
        # except for anchor sequences which are stored as lists and passed to functions


    def __str__(self):
        return "ncbi_taxid is %s, rfam_seq_id is %s, pfam_id is %s, organism is %s" % (self.ncbi_taxid, self.rfam_seq_id, self.pfam_id, self.organism)


#!!!How can we still allow PDB sequences to be used, new method of family identification
    #only works for the Author chain naming method
    #would it be possible to

aa_translation = {'ALA':'A', 'CYS':'C', 'ASP':'D', 'GLU':'E', 'PHE':'F', 'GLY':'G', 'HIS':'H','ILE':'I', 'LYS':'K',
'LEU':'L', 'MET':'M', 'ASN':'N', 'PRO':'P', 'GLN':'Q', 'ARG':'R', 'SER':'S', 'THR':'T', 'TRP':'W', 'TYR':'Y', 'VAL':'V'}


def read_rfam_to_pdb_mappings(download=False):
    rfam_family_to_chains = {}
    rfam_chain_to_family = {}
    pdb_chain_to_rfam_family = {}

    filename = os.path.join(BASE_DIRECTORY,'pdb_chain_to_best_rfam.txt')

    if download or not os.path.exists(filename):
        url = 'https://rna.bgsu.edu/data/pdb_chain_to_best_rfam.txt'
        urlretrieve(url,filename)

    lines = open(filename,read_mode).read().split("\n")

    for line in lines:
        #print(line)
        if line.startswith('RF'):
            fields = line.split()

            # check for Rfam positions out of order; it happens and we don't know why
            # it seems that this is a match of Infernal to the opposite strand of the DNA
            # but that doesn't make sense ... so skip these
            if int(fields[3]) > int(fields[4]):
                #print("Out of order: %s" % line.replace("\n",""))
                continue

            family = fields[0]
            if not family in rfam_family_to_chains:
                rfam_family_to_chains[family] = []

            chain_id = '_'.join([fields[1].upper(), fields[2], fields[3], fields[4]])
            rfam_family_to_chains[family].append(chain_id)
            rfam_chain_to_family[chain_id] = family
            pdb_chain = fields[1].upper() + "|1|" + fields[2]
            pdb_chain_to_rfam_family[pdb_chain] = family

    # for joint alignments like SSU and LSU, add together chain lists across all families
    # in the joint alignment
    for molecule in joint_alignments.keys():
        rfam_family_to_chains[molecule] = []
        for family in joint_alignments[molecule]:
            rfam_family_to_chains[molecule] += rfam_family_to_chains[family]

    return rfam_family_to_chains, rfam_chain_to_family, pdb_chain_to_rfam_family

def split_unit_id_string(unit_id_string):

    return re.split(',:',unit_id_string)

def detect_presence_of_proteins(all_unit_ids):
    #this function determines whether any of the unit IDs in the input have proteins

    ###global PRESENCE_OF_PROTEINS
    protein_found = False

    for unit_id in all_unit_ids:
        for key in aa_translation.keys():
            if key in unit_id:
                protein_found = True

    return protein_found

def get_loop_id_to_unit_id(PDB_ids):
    # load or download list of loops and their nucleotides

    loop_id_to_unit_id = {}
    loop_id_to_breaks = {}

    for PDB_id in PDB_ids:
        local_file = os.path.join(LOOP_DIRECTORY,"%s_loop_list.csv" % PDB_id)

        if not os.path.exists(local_file):
            try:
                url = "https://rna.bgsu.edu/rna3dhub/loops/download/%s" % PDB_id
                url = "https://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/%s" % PDB_id
                print("Downloading loop lists for %s using %s" % (PDB_id,url))
                urlretrieve(url,local_file)
            except:
                print("Unable to download loop list for %s" % PDB_id)

        if os.path.exists(local_file):
            with open(local_file,read_mode) as fid:
                loop_list_text = fid.read()
            for line in loop_list_text.splitlines():
                fields = line.split('"')
                if len(fields) >= 5:
                    loop_id = fields[1]
                    loop_id_to_unit_id[loop_id] = fields[3]  # text string of unit ids separated by commas
                    loop_id_to_breaks[loop_id] = fields[5]   # text string of 1,0,0,1 to indicate start and end of strands

    return loop_id_to_unit_id, loop_id_to_breaks


def convert_loop_id_to_unit_ids(loop_id,unit_id_string,border_string):
    """
    convert an individual loop id to ranges of unit ids

    loop_id is a BGSU RNA loop id like HL_1S72_001 or IL_1S72_003
    unit_id_list is the list of unit ids for that loop, already retrieved

    For HL, return a list with one string with unit ids separated by colons.
    For IL, J3, etc., return a list with one string for each strand.
    """

    if loop_id.startswith("HL"):
        unit_ids = [unit_id_string.replace(",",":")]  # list of one text string separated by :

    else:
        unit_ids = []                # list of unit ids or ranges of unit ids

        unit_id_list = unit_id_string.split(",")
        border_list = [int(b) for b in border_string.split(",")]

        unit_id = unit_id_list[0]    # current string being built
        border_total = 1             # how many border nucleotides have been hit

        for i in range(1,len(unit_id_list)):

            if border_total % 2 == 1:
                unit_id += ":" + unit_id_list[i]      # connect with : to indicate a range
            else:
                unit_ids.append(unit_id)      # end current string and append to list
                unit_id = unit_id_list[i]     # non-consecutive, start new string

            border_total += border_list[i]

        unit_ids.append(unit_id)              # append last id or range

    return unit_ids


def convert_loop_ids_to_unit_ids(input_ids):
    # determines if the input IDs are hairpin, internal, J3, ... loop IDs from the BGSU RNA site
    # as list of loop ids and their nucleotides
    # converts HL to one range of ids, converts IL to two ranges, J3 to three ranges, etc.
    # input_ids could be [["IL_4V9F_004"]] or [["IL_4V9F_004"],["IL_4V9F_007"]]
    # or [["4V9F|1|0|...:4V9F|1|0|"],["IL_4V9F_005"]]
    # Also replaces comma-separated lists of unit ids with a list of unit id strings

    loop_types = ["HL","IL","J3","J4","J5","J6","J7","J8","J9"]

    # determine which loop lists are needed
    PDB_list = set([])

    for id_list in input_ids:               # id_list is a list of strings
        for id_range in id_list:            # ids possibly separated by :
            for id in id_range.split(","):  # individual ids might have been separated by commas
                for loop_type in loop_types:
                    if id.startswith(loop_type):
                        PDB_list.add(id.split("_")[1].upper())  # extract PDB ID

    PDB_list = list(PDB_list)

    # load or download map of loop id to unit ids
    loop_id_to_unit_id, loop_id_to_breaks = get_loop_id_to_unit_id(PDB_list)

    # convert loop ids to lists and ranges of unit ids
    input_unit_ids = []      # list of lists of unit ids possibly connected with :

    for id_list in input_ids:
        new_list = []
        for id_string in id_list:                # text string of unit ids or possibly loop id or ids, possibly range
            for id in id_string.split(","):      # in case of comma-separated ids like "IL_4V9F_002,4V9F|1|0|G|25"
                if id in loop_id_to_unit_id:  # this is a loop id
                    unit_ids = convert_loop_id_to_unit_ids(id,loop_id_to_unit_id[id],loop_id_to_breaks[id]) # list of unit ids or ranges
                    new_list += unit_ids      # HL add one range to list, IL add two, etc.
                elif "|" in id:
                    new_list.append(id)       # this is a unit id
        input_unit_ids.append(new_list)

    all_unit_ids = set([])   # set of all unit ids
    for id_list in input_unit_ids:
        for id_range in id_list:
            if ":" in id_range:
                for id in id_range.split(":"):
                    all_unit_ids.add(id)
            elif "," in id_range:
                for id in id_range.split(","):
                    all_unit_ids.add(id)
            else:
                all_unit_ids.add(id_range)

    return input_unit_ids, list(all_unit_ids)


def extract_PDB_and_chains(all_unit_ids):

    all_PDB_ids = []
    all_chain_ids = []

    for unit_id in all_unit_ids:
        fields = unit_id.split("|")
        all_PDB_ids.append(fields[0])
        all_chain_ids.append("|".join(fields[0:3]))

    return all_PDB_ids, all_chain_ids


def identify_anchor_PDB_structures(input_unit_ids):
    # input_unit_ids is a list of lists
    # develops a list "anchor_PDB_structures" for later use from input ids in author format
    # one anchor PDB per list in input_unit_ids

    for unit_id_list in input_unit_ids:
        pdb_structure, other = ("|".join(unit_id_list)).split("|", 1)
        ANCHOR_PDB_STRUCTURES.append(pdb_structure.lower())

    return ANCHOR_PDB_STRUCTURES


def map_seqres_chains_to_sequences():
    """
    create a dictionary that maps chain ids (found in headers) to sequences
    Header format:
    >100d_A mol:na length:10  DNA/RNA (5'-R(*CP*)-D(*CP*GP*GP*CP*GP*CP*CP*GP*)-R(*G)-3')
    """

    lines = gzip.open('pdb_seqres.txt.gz',read_mode).read().split("\n")

    chain_to_sequence = {}

    # every pair of lines is a header followed by a sequence
    for i in range(0, len(lines), 2):
        if i+1 < len(lines):
            header = lines[i]
            chain = header.split(" ")[0].replace(">","")
            sequence = lines[i+1]
            chain_to_sequence[chain] = sequence

    return chain_to_sequence


def obtain_anchor_PDB_sequences(input_unit_ids, seqres_lines):
    # determines anchor sequences

    global PRESENCE_OF_PROTEINS

    if PRESENCE_OF_PROTEINS == True:
        protein_anchor_sequences = []
        protein_struc_header = []
        #these are only necessary if proteins are present in input ids (not for internal loop queries)

    rna_anchor_sequences = []

    headers, sequences = process_lines(seqres_lines)
    for series in input_unit_ids:
        for unit_id in series:
            struc, model, chain, monomer, position = unit_id.split("|")
            entry_identifier = struc.lower() + "_" + chain + " "
            # the end space helps eliminate "4v4q_A" or "4v4q_AA" from both being returned
            # when 4v4q_A is the actual entry_identifier

            if PRESENCE_OF_PROTEINS == True:
                for value in aa_translation.keys():
                    if value == monomer: #this logic checks if it is an amino acid
                        for i, header in enumerate(headers):
                            if entry_identifier in header:
                                temp_header = struc + " " + sequences[i]
                                protein_anchor_sequences.append(sequences[i])
                                protein_struc_header.append(temp_header)

                    else:
                        for i, header in enumerate(headers):
                            if entry_identifier in header:
                                rna_anchor_sequences.append(sequences[i])
            else:
                for i, header in enumerate(headers):
                    if entry_identifier in header:
                        rna_anchor_sequences.append(sequences[i])

    if PRESENCE_OF_PROTEINS == True:
        return protein_anchor_sequences, protein_struc_header, rna_anchor_sequences
    else:
        return rna_anchor_sequences


def map_unit_id_to_rfam(input_unit_ids, rfam_family_to_chains, unit_id_to_sequence_id={}, use_joint = False):
    """
    input_unit_ids is a list of lists of unit ids possibly separated by :
    rfam_family_to_chains is a dictionary mapping Pfam and Rfam families to chains and ranges
    use_joint = False prevents matches to joint alignments
    """

    unit_id_string_to_family_and_chain = {}  # dictionary mapping a unit id string to Pfam/Rfam family

    for unit_id_list in input_unit_ids:
        for unit_id_string in unit_id_list:
            chains = set([])   # make sure the range here only includes one PDB ID and chain

            if "," in unit_id_string:
                unit_ids = unit_id_string.split(",")
            else:
                unit_ids = unit_id_string.split(":")

            for unit_id in unit_ids:
                #struc, model, chain, monomer, position = unit_id.split("|") #4y4o|1|A|U|118 (example)
                fields = unit_id.split("|")     #4y4o|1|A|U|118 (example)
                chains.add(fields[0].upper() + "_" + fields[2])  # 1ABC_W for example

            if len(chains) > 1:
                print("Error: unit id string %s has more than one chain" % unit_id_string)
                continue

            # now it's clear there is only one chain in unit_id_string

            chain = list(chains)[0]

            # loop over Pfam/Rfam families and chains, record where the chain is listed
            possible_matches = []
            for family, corresponding_ids in rfam_family_to_chains.items():
                if not use_joint and family in joint_alignments.keys():
                    continue
                else:
                    for chain_and_range in corresponding_ids:
                        if chain in chain_and_range:
                            possible_matches.append((family,chain_and_range))

            print("RPFAM: possible matching Rfam families for %s are %s" % (unit_id_string,possible_matches))

            F = set([])   # identified families for this unit id string
            C = set([])   # identified chain and range for this unit id string

            for unit_id in unit_ids:
                #struc, model, chain, monomer, position = unit_id.split("|") #4y4o|1|A|U|118 (example)
                fields = unit_id.split("|") #4y4o|1|A|U|118 (example)
                if len(fields) < 5:
                    # not using a unit id, just trust it
                    for family, chain_and_range in possible_matches:
                        F.add(family)
                        C.add(chain_and_range)

                else:
                    # using a unit id, check that its sequence position is in the right range
                    if unit_id in unit_id_to_sequence_id:
                        sequence_id = unit_id_to_sequence_id[unit_id]
                        position = int(sequence_id.split("|")[4])
                    else:

                        # print('RPFAM: did not find sequence id for %s' % unit_id)

                        # use nucleotide number, which will not work in all cases; nucleotide numbers are fragile
                        try:
                            position = int(fields[4])
                        except:
                            position = -1  # just in case

                    # check that residue number (position) is within the range detected by Rfam; reduce multiple matches
                    for family, chain_and_range in possible_matches:
                        _struc, _chain, mon_begin, mon_end = chain_and_range.split("_")

                        if int(mon_begin) <= position and position <= int(mon_end):
                            F.add(family)
                            C.add(chain_and_range)

            F = list(F)
            C = list(C)

            if len(F) == 1:
                unit_id_string_to_family_and_chain[unit_id_string] = (F[0],C[0]) # use the one family found
            elif len(F) == 0:
                print("No Rfam family found for these unit ids: %s" % unit_id_string)
                # return no match, deal with the lack of a match in later functions
                # unit_id_string_to_family_and_chain[unit_id_string] = (None,None)
            elif len(F) > 1:
                print("More than one family found for these unit ids!")
                print(unit_id_list)
                print("Found matches to:")
                print(F)
                print(C)
                print("Using:")
                print(F[0])
                print(C[0])
                unit_id_string_to_family_and_chain[unit_id_string] = (F[0],C[0]) # use the first family found

    return unit_id_string_to_family_and_chain


def determine_xfam_families(input_unit_ids, unit_id_to_sequence_id = {}):
    # for proteins and separately for RNA, map chains to Rfam/Pfam families, return lists

    rfam_family_to_chains, rfam_chain_to_family, pdb_chain_to_rfam_family = read_rfam_to_pdb_mappings()
    #rfam_family_to_chains.update = Pfam_family_mapping

    unit_id_string_to_family_and_chain = map_unit_id_to_rfam(input_unit_ids, rfam_family_to_chains, unit_id_to_sequence_id)

    family_chain = set([])

    for family,chain in unit_id_string_to_family_and_chain.values():
        family_chain.add((family,chain))

    return unit_id_string_to_family_and_chain, family_chain


def determine_families(input_unit_ids):
    # for proteins and separately for RNA, map chains to Rfam/Pfam families, return lists

    global PRESENCE_OF_PROTEINS
    if PRESENCE_OF_PROTEINS == True:
        Pfam_families, protein_chains = map_unit_id_to_rfam(input_unit_ids, Pfam_family_mapping)
    else:
        Pfam_families = [None for unit_id_list in input_unit_ids]
        protein_chains = [None for unit_id_list in input_unit_ids]

    # rfam_family_to_chains is already imported
    Rfam_families, RNA_chains = map_unit_id_to_rfam(input_unit_ids, rfam_family_to_chains)

    #import pdb; pdb.set_trace()

    return Rfam_families, RNA_chains, Pfam_families, protein_chains


def open_family_files(RF, PF):
    rna_info = {}
    with open(RF + "ncbi_acc.csv", "r") as r_file:
        for line in r_file:
            rfam_seq_id, ncbi_taxid = line.split(",", 1)
            rna_info[ncbi_taxid] = rfam_seq_id
            match = PhasedMatch(ncbi_taxid, rfam_seq_id)
            phased_alignment.append(match)

    with open(PF + "ncbi_acc.csv", "r") as p_file:
        for line in p_file:
            pfam_id, ncbi_taxid = line.split(",", 1)
            if ncbi_taxid in rna_info.keys():
                rfam_seq_id = rna_info[ncbi_taxid]
                match = PhasedMatch(ncbi_taxid.strip(), rfam_seq_id, pfam_id)
            else:
                match = PhasedMatch(ncbi_taxid.strip(), None, pfam_id)
            phased_alignment.append(match)


def get_equivalence_classes(resolution):
    """
    resolution is a text string like '2.5A' or 'all'
    Download and read all equivalence classes at the current resolution.
    """

    ife_to_equivalence_class = {}
    equivalence_class_to_ifes = {}
    chain_to_equivalence_class_and_position = {}

    local_file = os.path.join(BASE_DIRECTORY,'alignments','equivalence_class_%s.txt' % resolution)

    # equivalence classes change every week, try to download to update
    try:
        url = 'https://rna.bgsu.edu/rna3dhub/nrlist/download/current/%s/csv' % resolution
        print("Getting data from %s" % url)
        urlretrieve(url,local_file)
    except:
        print('Unable to download equivalence classes from %s' % url)

    if os.path.exists(local_file):
        f = open(local_file,read_mode)
        lines = f.readlines()
        f.close()

        for line in lines:
            # print(line)
            fields = line.split('"')
            class_name = fields[1]
            ifes = fields[5].split(",")

            equivalence_class_to_ifes[class_name] = ifes

            for i in range(0,len(ifes)):
                ife = ifes[i]
                ife_to_equivalence_class[ife] = class_name
                for chain in ife.split("+"):
                    chain_to_equivalence_class_and_position[chain] = (class_name,i)

    else:
        print('No local file %s of equivalence classes found' % local_file)

    return ife_to_equivalence_class, equivalence_class_to_ifes, chain_to_equivalence_class_and_position

def write_file_protein_sequences_with_anchors(protein_struc_header, PF):
    #if UPDATE:
        PDB_structure_orgs = []
        #print protein_struc_header

        with open(PF + "_sequences_with_anchors.fa", "w") as out:
            for i, spliced_seq in enumerate(protein_struc_header):
                structure, sequence = spliced_seq.split(" ")
                out.write(">" + structure + " [PDB]\n" + sequence + "\n")
            for i, match in enumerate(phased_alignment):
                if match.pfam_id:
                    out.write(">" + str(match.pfam_id) + " " + str(match.organism) + str(match.pfam_sequence) + "\n")


        print("hmmer input file created")
    # Creates readable UniProt_sequences_file.txt for when update = false
    # This file contains UniProt IDs, organisms, and protein sequences


def hmmbuild(PF):
    # Builds HMM for the PFAM family
    ## Check the span of the HMMER model Pfam uses
    os.system("hmmbuild -O " + PF + ".seed.sto " + PF + ".hmm " + PF + "_seed_from_PFAM.sto")


def hmmalign_sequences(PF):
    #Aligns sequences to HMM from hmmbuild
    os.system('hmmalign ' + PF +'.hmm ' + PF + "_sequences_with_anchors.fa > " + PF + "_alignments_with_anchors.sto")

def merge_alignments_protein(PF):
    # Merges alignments based on Reference Annotation
    if PF + ".seed.sto" == None:
        os.system('wget -O ' + PF+ '.seed.sto https://pfam.org/family/' + PF + '/alignment')
    os.system('./esl-alimerge ' + PF + "_alignments_with_anchors.sto " + PF + ".seed.sto > " + PF + "merged.sto")


def append_hmmaligned_sequences(protein_anchor_sequences, PF):
    #at the bottom, I detail the procedure for obtaining this alignment file
    pdb_anchors_aligned_protein = []
    protein_full_sequence_records = AlignIO.read(PF +"merged.sto", "stockholm")
    #print(protein_full_sequence_records)
    for i, structure in enumerate(ANCHOR_PDB_STRUCTURES):
        for record in protein_full_sequence_records:
            if structure in record.name:
                #print record.name
                pdb_anchors_aligned_protein.append(record.name + "|" + record.seq)


    for record in protein_full_sequence_records:
        if "|" in record.name:
            number, name_id = str(record.id).split("|", 1)
            if "|" not in name_id:
                for i, match in enumerate(phased_alignment):
                        match.pfam_aligned_sequence = record.seq
#!!!temporary patch here, removes every other line with "if "|" not in name_id:" but we still duplicate entries!

    #print(pdb_anchors_aligned)
    print("hmmaligned sequences added to classes")
    #import pdb; pdb.set_trace()
    return pdb_anchors_aligned_protein
    #Implements aligned sequences to classes for output purposes













#-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------












def pfam_seed_id_and_taxid(PF):
    pfam_info = csv.reader(open(PF +"SQLquery.csv", 'r'))
    for pfam_id, ncbi_taxid in pfam_info:

        if pfam_id != 'pfamseq_id':
            match = PhasedMatch(ncbi_taxid, None, pfam_id)
            phased_alignment.append(match)


    # At the bottom I detail the SQL Pfam search query that obtains this result file.
    # This determines the ncbi_taxid and PFAM ID from the SQL query results

def obtain_pfam_sequences(PF):
    pfam_alignment = AlignIO.read(PF + "_seed.txt", "stockholm")
    for i, match in enumerate(phased_alignment):
        for record in pfam_alignment:
            if record.name == match.pfam_id:
                match.pfam_sequence = record.seq


def write_file_rnacentral_sequences_with_anchors(rna_anchor_sequences, RF):
    #if UPDATE:
        with open(RF + "_sequences_with_anchors.fa", "w") as out:
            for j, value in enumerate(rna_anchor_sequences):
                print(j)
                #print ANCHOR_PDB_STRUCTURES[j]
                #print rna_anchor_sequences[j]
                #out.write(">" + ANCHOR_PDB_STRUCTURES[j] + "\r\n" + rna_anchor_sequences[j] + "\r\n")
            import pdb; pdb.set_trace()
            for i, match in enumerate(phased_alignment):
                out.write(">" + str(match.rfam_seq_id) + "\r\n" + str(match.rfam_sequence) + "\r\n")
        print("Infernal input file created")

    #else:
    #    print("RNACentral_sequences_with_anchors.fa exists\nInfernal input file is prepared")
    # Creates RNACentral_sequences_with_anchors.fa for >PDB_structure\r\nrna_anchor_sequence\r\n
    # And >URS\r\nrnacentral_sequence\r\n
    # For Infernal input file

def cmalign_sequences(RF):
    print("Running cmalign")
    #os.system("pwd")
    os.system('cmalign ' + RF + '.cm ' + RF + '_sequences_with_anchors.fa > RFAM_' + RF + '_alignments_with_anchors.sto')

    # Aligns to existing RFAM covariance model

def merge_alignments_RNA(RF):
    #os.system('wget -O ' + RF + '.seed.sto https://rfam.org/family/' + RF + '/alignment')
    #os.system("cd infernal-1.1.3/miniapps/")
    #os.system("pwd")
    print("Running esl-alimerge")
    os.system("./esl-alimerge RFAM_" + RF +"_alignments_with_anchors.sto " + RF + ".seed.sto > " + RF + "merged.sto")
    print("esl-alimerge RFAM_" + RF +"_alignments_with_anchors.sto " + RF + ".seed.sto > " + RF + "merged.sto")
    # esl-alimerge RFAM_RF00177_alignments_with_anchors.sto RF00177.seed.sto > RF00177merged.sto
    # Merges alignments

def append_cmaligned_sequences(rna_anchor_sequences, RF):
    """
    not being used
    """


    #at the bottom, I detail the procedure for obtaining this alignment file
    #RNA_full_sequence_records = AlignIO.read(RF + "merged.sto", "stockholm")
    alignment_path_file = "alignments/rfam/" + RF + "_full_anchored.sto.gz"
    RNA_full_sequence_records = AlignIO.read(gzip.open(alignment_path_file), "stockholm")

    print(rna_anchor_sequences)

    #for record in RNA_full_sequence_records:
    #    print("%s | %s" % (record.name,record.seq))

    pdb_anchors_aligned_rna = []
    for i, structure in enumerate(ANCHOR_PDB_STRUCTURES):
        for record in RNA_full_sequence_records:
            if structure in record.name:
                pdb_anchors_aligned_rna.append(record.name + "|" + record.seq)

    for i, match in enumerate(phased_alignment):
        for record in RNA_full_sequence_records:
            if "|" in record.name:
                n, record_id = str(record.name).split("|")

            if match.rfam_seq_id == record.name:
                match.rfam_aligned_sequence = record.seq

#    for record in RNA_full_sequence_records:
 #       for j, value in enumerate(RNA_id_of_interest):
  #          if record.name ==  RNA_id_of_interest[j]:
   #             continue
    #    n, structure_id = str(record.name).split("|")
     #   RNA_ecoli_sequence_of_interest = record.seq
      #  return RNA_ecoli_sequence_of_interest
    print("cmalinged sequences re-added")
    return pdb_anchors_aligned_rna




#------------------------------------------------------------------------------------------



def create_Phased_id():
    for i, match in enumerate(phased_alignment):
        if match.ncbi_taxid:
            if match.rfam_seq_id and match.uniprot_id:
                if match.uniprot_id != "None":
                    match.phased_id_RNA = match.ncbi_taxid + "/" + match.rfam_seq_id
                    match.phased_id_protein = match.ncbi_taxid + "/" + match.uniprot_id

            elif match.URS and match.pfam_id:
                match.phased_id_RNA = match.ncbi_taxid + "/" + match.URS
                match.phased_id_protein = match.ncbi_taxid + "/" + match.pfam_id


#----------------------------------------------------------------------------------------------------------------



def stockholm_individual_headers(pdb_anchors_aligned, outfile):
    for i, match in enumerate(pdb_anchors_aligned):
        anchor_num, anchor_id, anchor_sequence = pdb_anchors_aligned[i].split("|")
        sequence_GS_id = str("#=GS " + anchor_id).ljust(45, ' ')
        description = str("DE " + (str("PDB STRUCTURE " + anchor_id))).ljust(70, ' ')
        outfile.write(sequence_GS_id + description + "\n")
    return outfile



def stockholm_individual_body(pdb_anchors_aligned, outfile):
    for i, aligned_anchor in enumerate(pdb_anchors_aligned):
        anchor_num, anchor_id, anchor_sequence = aligned_anchor.split("|")
        length = len(str(anchor_sequence))
        single_alignment = str(anchor_id).ljust(39, ' ') + (str(anchor_sequence).ljust(length + 25, ' '))
        outfile.write(single_alignment + "\n")
    return outfile

def output_organism_information(phased_id, organism, ncbi_taxid, outfile):
    sequence_GS_id = str("#=GS " + phased_id).ljust(45, ' ')
    if "None" not in organism:
        description = str("DE " + str(organism.strip())).ljust(70, ' ')
    else:
        org_url = ("https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=" + ncbi_taxid + "&lvl1")
        org_get = requests.get(org_url)

        if get_organisms:
            soup = BeautifulSoup(org_get.text, 'html.parser')

        junk, organism_and_parenth = str(soup.find("title")).split("(")
        organism_info, junk = organism_and_parenth.split(")")
        description = str("DE " + str(organism_info)).ljust(70, ' ')
    outfile.write(sequence_GS_id + description + "\n")
    return outfile

def forming_individual_stockholm_files(PF, RF, pdb_anchors_aligned_protein, pdb_anchors_aligned_rna):
    if OUTPUT_TYPE == "coupled":
        with open(RF + "_coupled_full.sto", 'w') as outfile_RNA:
            outfile_RNA.write("# STOCKHOLM 1.0" + "\n\n")

            with open(PF + "_coupled_full.sto", 'w') as outfile_protein:
                outfile_protein.write("# STOCKHOLM 1.0" + "\n\n")
                outfile_protein = stockholm_individual_headers(pdb_anchors_aligned_protein, outfile_protein)
                outfile_RNA = stockholm_individual_headers(pdb_anchors_aligned_rna, outfile_RNA)

                ncbi_list = []

                for i, match in enumerate(phased_alignment):
                    if match.ncbi_taxid not in ncbi_list:
                        ncbi_list.append(match.ncbi_taxid)

                        outfile_protein = output_organism_information(match.phased_id_protein, match.organism, match.ncbi_taxid, outfile_protein)
                        outfile_RNA = output_organism_information(match.phased_id_RNA, match.organism, match.ncbi_taxid, outfile_RNA)

                outfile_protein.write("\n")
                outfile_RNA.write("\n")
                #writes out the GS ids of non-anchoring sequences

                outfile_protein = stockholm_individual_body(pdb_anchors_aligned_protein, outfile_protein)
                outfile_RNA = stockholm_individual_body(pdb_anchors_aligned_rna, outfile_RNA)


                ncbi_list = []

                for i, match in enumerate(phased_alignment):
                    if match.ncbi_taxid not in ncbi_list:
                        if match.rnacentral_aligned_sequence and match.uniprot_aligned_sequence:

                            ncbi_list.append(match.ncbi_taxid)

                            if match.uniprot_aligned_sequence:
                                length = len(str(match.uniprot_aligned_sequence))
                                single_alignment = (str(match.phased_id_protein).ljust(39, ' ')) + (str(match.uniprot_aligned_sequence).ljust(length + 25, ' '))
                                outfile_protein.write(single_alignment + "\n")

                            elif match.pfam_aligned_sequence:
                                length = len(str(match.pfam_aligned_sequence))
                                single_alignment = (str(match.phased_id_protein).ljust(39, ' ')) + (str(match.pfam_aligned_sequence).ljust(length + 25, ' '))
                                outfile_protein.write(single_alignment + "\n")

                            length = len(str(match.rnacentral_aligned_sequence))
                            single_alignment = (str(match.phased_id_RNA).ljust(39, ' ')) + (str(match.rnacentral_aligned_sequence).ljust(length + 25, ' '))
                            outfile_RNA.write(single_alignment + "\n")
                #writes out the output alignments
                end_notation = "//"
                outfile_protein.write(end_notation)
                outfile_protein.close()
                print("Protein coupled alignment file is complete")

            outfile_RNA.write(end_notation)
            outfile_RNA.close()
            print("RNA coupled alignment file is complete")
            #End of output files


def forming_RNA_only_query_stockholm_file(RF, pdb_anchors_aligned_rna):
    #!!! this is used for the internal loop query - current code is filler
    with open(RF + "_IL_ouput.sto", "w") as outfile_RNA:
        for line in RF_full_family_alignment:
            print(line)

def map_unit_ids_to_sequence_positions(ANCHOR_PDB_STRUCTURES):
    """
    Not needed for RNA structures; just use the API as in map_unit_ids_to_sequence_positions_chains
    loop over all required PDB structures and map unit ids to sequence positions and vice versa
    """

    unit_id_to_sequence_id = {}
    sequence_id_to_unit_id = {}

    for structure_to_map in set(ANCHOR_PDB_STRUCTURES):

        cif_path_file = CIF_FILE_LOCATION + structure_to_map + ".cif"

        # download .cif file if not already downloaded
        if not os.path.exists(cif_path_file):
            pdb_mmcif_download_url = ("https://files.rcsb.org/view/%s.cif" % structure_to_map)
            mmcif_get = requests.get(pdb_mmcif_download_url)
            mmcif_content = mmcif_get.text
            with open(cif_path_file, "w") as f:
                f.write(mmcif_content)
            f.close()
            print("Downloaded %s mmcif file for mapping" % structure_to_map)

        # load mmcif data file
        ifh = open(cif_path_file)
        pRd = PdbxReader(ifh)

        # read mmcif data
        data = []
        attempt = pRd.read(data)

        # propagate list data with one or more DataContainer objects
        block = data[0]

        # obtain block of interest
        mapping = block.getObj("pdbx_poly_seq_scheme") #_pdbx_poly_seq_scheme

        # loop over lines of data in the block
        for i, line in enumerate(mapping):
            pdb_chain_id = mapping.getValue("asym_id", i) #pdb-given chain id for each residue
            auth_chain_id = mapping.getValue("pdb_strand_id", i) #author chain id for each residue
            resolved_seq_num = mapping.getValue("seq_id", i) #numbering found in PDB unit id
            auth_seq_num = mapping.getValue("auth_seq_num", i) #author sequence number/letters combo
            ndb_seq_num = mapping.getValue("ndb_seq_num", i) #experimental sequence numbering
            auth_mon_id = mapping.getValue("auth_mon_id", i)
            pdb_mon_id = mapping.getValue("pdb_mon_id", i) #monomer: 3 letter if protein, 1 letter nucleotide
            #print(len(pdb_mon_id))
            #if auth_chain_id == '1a' and auth_seq_num == '1240' and auth_mon_id == 'U':
            #    print resolved_seq_num, auth_seq_num, ndb_seq_num

            # assume that we want model 1; if not, we'll get a key error
            unit_id = structure_to_map.upper() + "|1|" + auth_chain_id + "|" + auth_mon_id + "|" + auth_seq_num
            sequence_id = structure_to_map.upper() + "|sequence|" + auth_chain_id + "|" + pdb_mon_id + "|" + resolved_seq_num

            # make mappings both directions for all unit ids and sequence ids
            unit_id_to_sequence_id[unit_id] = sequence_id
            sequence_id_to_unit_id[sequence_id] = unit_id

    return unit_id_to_sequence_id, sequence_id_to_unit_id

def map_unit_ids_to_sequence_positions_chains(all_chain_ids,unit_id_to_sequence_id = {},sequence_id_to_unit_id = {}):
    """
    loop over all required PDB structures and map unit ids to sequence positions and vice versa

    We're developing this on a PC and using OneDrive and Dropbox to share files.
    Unfortunately those systems do not allow case-dependent filenames.
    So if a windows machine is detected, the local filenames will put a hyphen
    in front of lowercase letters in the chain identifier.
    """

    distinct_chain_ids = list(set(all_chain_ids))

    if 'win32' in sys.platform:
        convert_lowercase = True
    else:
        convert_lowercase = False

    for chain_id in distinct_chain_ids:
        fields = chain_id.split("|")

        chain = fields[2]

        if convert_lowercase:
            temp = chain
            chain = ''
            for c in temp:
                if not c == c.upper():
                    chain += "-"
                chain += c

        directory = os.path.join(BASE_DIRECTORY,'alignments','sequence_to_unit_id')
        if not os.path.exists(directory):
            os.mkdir(directory)

        local_file = os.path.join(directory,fields[0]+"_"+chain+".txt")

        if not os.path.exists(local_file):

            try:
                url = 'https://rna.bgsu.edu/rna3dhub/rest/SeqtoUnitMapping?ife=%s' % chain_id
                #urlretrieve(url,local_file)
                url_data = requests.get(url)
            except:
                print("Unable to download sequence to unit id mapping for %s" % chain_id)

            f = open(local_file,write_mode)
            f.write(url_data.text.replace('</br>','\n'))
            f.close()

        f = open(local_file,read_mode)
        lines = f.readlines()
        f.close()

        for line in lines:
            correspondence = line.rstrip('\n').rstrip('\r').split(" ")
            if len(correspondence) == 3:
                sequence_id = correspondence[0]
                unit_id = correspondence[2]
                unit_id_to_sequence_id[unit_id] = sequence_id

                # multiple unit ids can map to the same sequence id, choose the best
                if sequence_id in sequence_id_to_unit_id:
                    if not unit_id == sequence_id_to_unit_id[sequence_id]:
                        if len(unit_id) < len(sequence_id_to_unit_id[sequence_id]):
                            print("RPFAM: 1 Prefer unit id %s over %s" % (unit_id,sequence_id_to_unit_id[sequence_id]))
                            sequence_id_to_unit_id[sequence_id] = unit_id
                        elif int(unit_id.split("|")[1]) < int(sequence_id_to_unit_id[sequence_id].split("|")[1]):
                            print("RPFAM: 2 Prefer unit id %s over %s" % (unit_id,sequence_id_to_unit_id[sequence_id]))
                            sequence_id_to_unit_id[sequence_id] = unit_id
                else:
                    sequence_id_to_unit_id[sequence_id] = unit_id

    return unit_id_to_sequence_id, sequence_id_to_unit_id

def look_up_sequence_ids(input_unit_ids, unit_id_to_sequence_id):
    # loop over list of lists of unit ids and convert to sequence ids
    # unit ids may be separated by :, then sequence ids will be also

    sequence_ids = []          # list of lists of sequence ids

    for unit_id_list in input_unit_ids:
        new_list = []
        for unit_id in unit_id_list:  # loop over unit ids, possibly ranges
            colon_list = []
            for id in unit_id.split(":"):  # split by colon if present
                colon_list.append(unit_id_to_sequence_id[id])
            new_list.append(":".join(colon_list))
        sequence_ids.append(new_list)

    return sequence_ids


def map_sequence_ids_to_columns(alignment,rfam_chain_id,start_index=1):
    """
    Map positions in the sequence from the sequence ids to the
    corresponding column in the alignment.
    Start at start_index in the sequence because some PDB chain
    sequences are longer than what matches the Rfam family.
    Only the range that Rfam says matches the family is put into
    the alignment, as of 2022-09-30.
    This function makes up the sequence ids from the position
    in the sequence and the base observed there.
    Modified nucleotides will be converted to parent nucleotides.
    Not sure what happens to bases like I that have no parent.
    """

    #print("Looking for chain %s in the alignment" % rfam_chain_id)

    # identify the row of the alignment corresponding to chain
    row_found = False
    j = 0
    while not row_found and j < len(alignment["accession"]):
        if rfam_chain_id in alignment["accession"][j]:
            row_found = True
            sequence = alignment["sequence"][j]
        else:
            j += 1

    if not row_found:
        print("Found and checked %d rows in the alignment" % j)
        print("No row was found in the alignment for chain %s" % rfam_chain_id)
        return {},{}

    sequence_id_to_column = {}
    column_to_sequence_id = {}

    if "_" in rfam_chain_id:
        fields = rfam_chain_id.split("_")   # split out PDB ID
        chain_id = fields[0].upper() + "|sequence|" + fields[1]
    else:
        fields = rfam_chain_id.split("|")   # split out PDB ID
        chain_id = fields[0].upper() + "|sequence|" + fields[2]

    sequence_position = start_index
    column = 0

    #print("Chain %s has alignment row" % chain_id)
    #print(sequence)

    for column in range(len(sequence)):
        monomer = sequence[column]
        if monomer != "-" and monomer != '.':
            sequence_id = chain_id + "|" + monomer.upper() + "|" + str(sequence_position)
            sequence_id_to_column[sequence_id] = column
            column_to_sequence_id[column] = sequence_id
            sequence_position += 1

    return sequence_id_to_column, column_to_sequence_id


def map_sequence_ids_to_alignment(alignment,rfam_chain_id,start_index=1):
    """
    This function is just like map_sequence_ids_to_columns but here alignment is
    of type <class 'Bio.Align.MultipleSeqAlignment'>

    Map positions in the sequence from the sequence ids to the
    corresponding column in the alignment.
    Start at start_index in the sequence because some PDB chain
    sequences are longer than what matches the Rfam family.
    Only the range that Rfam says matches the family is put into
    the alignment, as of 2022-09-30.
    This function makes up the sequence ids from the position
    in the sequence and the base observed there.
    Modified nucleotides will be converted to parent nucleotides.
    Not sure what happens to bases like I that have no parent.
    """

    print("Looking for chain %s in the alignment" % rfam_chain_id)

    # identify the row of the alignment corresponding to chain
    row_found = False

    for record in alignment:
        if rfam_chain_id in record.description:
            row_found = True
            sequence = record.seq
            #print(record)
            #print(record.id)
            #print(record.description)
            break

    if not row_found:
        print("Found and checked every row in the alignment")
        print("No row was found in the alignment for chain %s" % rfam_chain_id)
        return {}

    sequence_id_to_column = {}

    fields = rfam_chain_id.split("_")   # split out PDB ID
    chain_id = fields[0].upper() + "|sequence|" + fields[1]

    sequence_position = start_index
    column = 0

    #print("Chain %s has alignment row" % chain_id)
    #print(sequence)

    for column in range(len(sequence)):
        monomer = sequence[column]
        if monomer != "-":
            sequence_id = chain_id + "|" + monomer.upper() + "|" + str(sequence_position)
            sequence_id_to_column[sequence_id] = column
            sequence_position += 1

    return sequence_id_to_column


def developing_output_columns(residue_output_column, sequence_ids, column_mapping_struc, output_data):

    for i, sequence_id in enumerate(sequence_ids):
        counter = 0
        for letter in sequence_id:
            if letter == "|":
                counter += 1
        #print(counter)
        if counter == 5:
            struc, model, chain, monomer, position, sub_exp_position = sequence_id.split("|")
        elif counter == 4:
            struc, model, chain, monomer, position = sequence_id.split("|")


        count = column_mapping_struc[position]
        residue_columns = []
        for record in output_data:
            residue_columns.append(record.seq[count])
        residue_output_column[str(sequence_id)] = residue_columns ##add chains!! or entire unit ids
        #print sequence_id

    return residue_output_column


def load_alignment(family,filename=[]):
    """
    Given a Pfam or Rfam family name, load the alignment and extract species, NCBI taxid, and description
    This should only need to be done once per alignment.
    """

    read_RF00005_fasta = False

    if not filename:
        if family == 'RF00005':
            filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,"%s_full_anchored.fa" % family)
            read_RF00005_fasta = True
        else:
            filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,"%s_full_anchored.sto" % family)

    filename_gz = filename + ".gz"

    accession_list = []
    NCBI_taxid_list = []
    domain_list = []
    description_list = []
    sequence_list = []
    codon_list = []           # only non-empty for RF00005

    print("Time %s. Loading alignment for %s from %s.  " % (str(datetime.now()),family,filename_gz))

    if family.startswith("RF") or family in joint_alignments.keys():
        # see if the compressed alignment file is available locally
        if not os.path.exists(filename_gz):
            # if not available locally, download it
            # Example: https://rna.bgsu.edu/data/alignments/rfam/RF00002_full_anchored.sto.gz
            url = "https://rna.bgsu.edu/data/" + filename_gz
            url = url.replace("\\","/")
            print("Downloading alignment for %s from %s" % (family,url))

            urlretrieve(url,filename_gz)

        """
        if family == 'RF00005':

            f = gzip.open(filename_gz,read_mode)
            counter = 0
            while True:
                header = f.readline()  # omit > character, remove \n if present
                counter += 1
                if not header:
                    f.close()
                    return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list

#                    raise Exception('Got to 1094')

                header = header[1:].rstrip("\n")
                sequence = f.readline().rstrip("\n")

                fields = header.split()

                if len(fields) >= 5:
                    accession = fields[0]
                    codon = fields[1]
                    NCBI_taxid = fields[2]
                    domain = fields[3]
                    description = " ".join(fields[4:])
                    accession_list.append(accession)
                    NCBI_taxid_list.append(NCBI_taxid)
                    domain_list.append(domain)
                    description_list.append(description)
                    sequence_list.append(sequence)
                    codon_list.append(codon)
                else:
                    print(header)
                    print(fields)

                    accession = fields[0]
                    accession_list.append(accession)
                    sequence_list.append(sequence)
                    description_list.append(accession)  # current way to tell that this is a PDB sequence

                    NCBI_taxid_list.append('')
                    domain_list.append('')
                    codon_list.append('')
        """

        # read the compressed alignment file
        try:
            if read_RF00005_fasta:

                f = gzip.open(filename_gz,read_mode)
                counter = 0
                while True:
                    header = f.readline()  # omit > character, remove \n if present
                    counter += 1
                    if not header:
                        f.close()
                        return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list

                    header = header[1:].rstrip("\n")
                    sequence = f.readline().rstrip("\n")

                    fields = header.split()

                    if len(fields) >= 5:
                        accession = fields[0]
                        codon = fields[1]
                        NCBI_taxid = fields[2]
                        domain = fields[3]
                        description = " ".join(fields[4:])
                        accession_list.append(accession)
                        NCBI_taxid_list.append(NCBI_taxid)
                        domain_list.append(domain)
                        description_list.append(description)
                        sequence_list.append(sequence)
                        codon_list.append(codon)
                    else:
                        accession = fields[0]
                        accession_list.append(accession)
                        sequence_list.append(sequence)
                        description_list.append(accession)  # current way to tell that this is a PDB sequence

                        NCBI_taxid_list.append('')
                        domain_list.append('')
                        codon_list.append('')

            else:
                f = gzip.open(filename_gz,read_mode)     # need this format for Flask server
                alignment = AlignIO.read(f,"stockholm")
                f.close()
        except:
            print("Unable to read alignment file %s" % filename_gz)
            return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list

    else:
        if not filename:
            filename = os.path.join(PFAM_ALIGNMENT_DIRECTORY,"%s_full_anchored.sto" % family)
        alignment = AlignIO.read(filename, "stockholm")

    counter = 0
    for record in alignment:
        counter += 1
        if counter % 30000 == 0:
            print("Time %s. Read %6d alignment rows. " % (str(datetime.now()),counter))

        name = record.id           # accession with sequence range

        # split taxid and domain out of the description
        fields = record.description.split()

        if len(fields) >= 3:
            NCBI_taxid = fields[0]
            domain = fields[1]
            description = " ".join(fields[2:])
        else:
            NCBI_taxid = 'None'
            domain = 'None'
            description = record.description

        accession_list.append(name)
        NCBI_taxid_list.append(NCBI_taxid)
        domain_list.append(domain)
        description_list.append(description)
        sequence_list.append(str(record.seq))

    return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list


def load_fasta_alignment(family,filename=[]):
    """
    Given an Rfam family name and maybe a filename,
    load the alignment and extract species, NCBI taxid, and description
    This should only need to be done once per alignment.
    """

    if not filename:
        filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,"%s_combined.fa" % family)

    filename_gz = filename + ".gz"

    accession_list = []
    NCBI_taxid_list = []
    domain_list = []
    description_list = []
    sequence_list = []
    codon_list = []           # only non-empty for RF00005

    #print("Time %s. Loading alignment for %s from %s.  " % (str(datetime.now()),family,filename_gz))
    # print("Loading alignment for %s from %s.  " % (family,filename_gz))

    if family.startswith("RF") or family in joint_alignments.keys():
        # see if the compressed alignment file is available locally
        if not os.path.exists(filename_gz):
            # if not available locally, download it
            # Example: https://rna.bgsu.edu/data/alignments/rfam/RF00002_combined.fa.gz
            url = "https://rna.bgsu.edu/data/" + filename_gz
            url = url.replace("\\","/")
            print("Downloading alignment for %s from %s" % (family,url))

            urlretrieve(url,filename_gz)


        # read the compressed alignment file
        #try:
        if True:
            print("Reading alignment file %s" % filename_gz)
            f = gzip.open(filename_gz,read_mode)
            counter = 0
            while True:
                header = f.readline()  # the line that starts >
                if not header:
                    f.close()
                    return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list

                counter += 1

                header = header[1:].rstrip("\n") # omit > character, remove \n if present
                sequence = f.readline().rstrip("\n")

                fields = header.split()

                if len(fields) >= 4:
                    # sequence from Rfam
                    if False and family == 'RF00005':
                        accession = fields[0]
                        codon = fields[1]
                        NCBI_taxid = fields[2]
                        domain = fields[3]
                        description = " ".join(fields[4:])
                    else:

                        #print(fields)

                        accession = fields[0]
                        codon = ''
                        NCBI_taxid = fields[1]
                        domain = fields[2]
                        description = " ".join(fields[3:])

                    accession_list.append(accession)
                    NCBI_taxid_list.append(NCBI_taxid)
                    domain_list.append(domain)
                    description_list.append(description)
                    sequence_list.append(sequence)
                    codon_list.append(codon)
                elif len(fields) > 0:
                    # sequence from PDB chain
                    chains = fields[0].split(",")       # chain ids separated by ,
                    for chain in chains:
                        accession_list.append(chain)
                        sequence_list.append(sequence)
                        description_list.append(chain)  # current way to tell that this is a PDB sequence

                        NCBI_taxid_list.append('')
                        domain_list.append('')
                        codon_list.append('')
                else:
                    print('Not sure what is up with line %s' % counter)
                    print(fields)

        # except:
        #     print("Unable to read alignment file %s" % filename_gz)
        #     return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list

    return accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list


def identify_alignment_columns(sequence_ids, Rfam_families, RNA_chains, Pfam_families, protein_chains):
    """
    no longer used

    develops residue columns from alignment for output viewer based on specified positions
    Rfam_families is a list of Rfam families with the same length as sequence_ids
    mapped_structures: list of pdb ids, strings
    sequence_ids: list of experimental sequence position ids derived from input bgsu unit ids, strings
    """

    rna_output_column = {}       # a single dictionary that maps all sequence ids or ranges to alignment columns
    protein_output_column = {}

    RNA_accession_list = {}      # one list for each list in sequence_ids
    RNA_NCBI_taxid_list = {}    # one list for each list in sequence_ids
    RNA_description_list = {}  # one list for each list in sequence_ids

    protein_accession_list = {}      # one list for each list in sequence_ids
    protein_NCBI_taxid_list = {}    # one list for each list in sequence_ids
    protein_description_list = {}  # one list for each list in sequence_ids

    # loop over lists of sequence ids
    for j in range(len(sequence_ids)):

        sequence_id_list = sequence_ids[j]
        Rfam_family = Rfam_families[j]
        RNA_chain = RNA_chains[j]
        Pfam_family = Pfam_families[j]
        protein_chain = protein_chains[j]

        if Rfam_family:
            filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,"%s_full_anchored.sto" % family)
            rna_alignment = AlignIO.read(filename, "stockholm")
        else:
            rna_alignment = []

        if Pfam_family:
            filename = os.path.join(PFAM_ALIGNMENT_DIRECTORY,"%s_full_anchored.sto" % family)
            protein_alignment = AlignIO.read(filename, "stockholm")
        else:
            protein_alignment = []

        RNA_accession_list[j] = []           # one entry for each row of the alignment
        RNA_NCBI_taxid_list[j] = []
        RNA_description_list[j] = []

        protein_accession_list[j] = []
        protein_NCBI_taxid_list[j] = []
        protein_description_list[j] = []

        rna_sequences_to_map = []
        protein_sequences_to_map =[]
        rna_structure_to_sequence = {}
        protein_structure_to_sequence = {}

        anchor_rna_sequence = ""

        for record in rna_alignment:
            name = record.name
            if "/" in name:
                name, NCBI_taxid = name.split("/")
            else:
                NCBI_taxid = None

            RNA_accession_list[j].append(name)
            RNA_NCBI_taxid_list[j].append(NCBI_taxid)
            RNA_description_list[j].append(record.description)

            # if this line of the alignment is the RNA chain we need
            if RNA_chain in record.name:
                rna_anchor_sequence = record.seq

        if len(rna_anchor_sequence) == 0:
            print("No RNA anchor sequence for chain %s found in alignment %s_full_anchored.sto" % (RNA_chain,Rfam_family))
        else:
            print("Found RNA anchor sequence for chain %s:" % RNA_chain)
            print(rna_anchor_sequence)

        print("Dictionary mapping structure to sequence")
        print(rna_structure_to_sequence)

        for record in protein_alignment:
            name = record.name
            if "/" in name:
                name, NCBI_taxid = name.split("/")
            else:
                NCBI_taxid = None

            protein_accession_list[j].append(name)
            protein_NCBI_taxid_list[j].append(NCBI_taxid)
            protein_description_list[j].append(record.description)

            # if this line of the alignment is the RNA chain we need
            if protein_chain in record.name:
                protein_anchor_sequence = record.seq

        r_exp_position = []
        r_ids = []

        p_exp_position = []
        p_ids = []

        for i, sequence_id in enumerate(sequence_id_list):
            fields = sequence_id.split(":").split("|")
            if len(fields) == 6:
                struc, model, chain, monomer, exp_position, sub_exp_position = fields
            if len(fields) == 5:
                struc, model, chain, monomer, exp_position = fields

            #if type(exp_position) != int(exp_position):
                #split somehow for letter ids not separated by bar

            if monomer in aa_translation.keys():   # must be an amino acid
                if exp_position not in p_exp_position:
                    #residue_type = "protein"
                    p_exp_position.append(exp_position)
                    p_ids.append(sequence_id)
            else:
                if exp_position not in r_exp_position:
                    #residue_type = "rna"
                    r_exp_position.append(exp_position)
                    r_ids.append(sequence_id)

            """
            for amino_acid_three_letter in aa_translation.keys():
                #print(amino_acid_three_letter)
                if amino_acid_three_letter in monomer:
                    if structure == struc:
                        if exp_position not in p_exp_position:
                            #residue_type = "protein"
                            p_exp_position.append(exp_position)
                            p_ids.append(sequence_id)
                else:
                    if structure == struc:
                        if exp_position not in r_exp_position:
                        #residue_type = "rna"
                            r_exp_position.append(exp_position)
                            r_ids.append(sequence_id)
            """

        print("Sequence ids")
        print(r_ids)
        print(p_ids)
        print("Sequence positions")
        print(r_exp_position)
        print(p_exp_position)

        #import pdb; pdb.set_trace()

        #obtains experimental positions from the mapped unit ids
        #convert these from variables to lists for multiple pdb structure requests
        #rna_sequences_to_map r_exp_position r_column_mapping_struc

        r_column_mapping_struc = []
        p_column_mapping_struc = []

        r_column_counter, r_column_mapping_struc = mapping_columns_with_gaps(rna_anchor_sequence, r_exp_position)
        #p_column_counter, p_column_mapping_struc = mapping_columns_with_gaps(protein_structure_to_sequence, p_exp_position)
        #determines the position in alignment accounting for gaps

        rna_output_column = developing_output_columns(rna_output_column, r_ids, r_column_mapping_struc, rna_alignment)
        #protein_output_column = developing_output_columns(protein_output_column, p_ids, p_column_mapping_struc, protein_alignment)
        #develops monomer output columns for each position

    return rna_output_column, protein_output_column, accession_list, description_list


def display_columns_file(sequence_ids, rna_output_column, protein_output_column, all_accession_list, all_description_list, unit_id_to_sequence_id):
    # sequence_ids is a list of lists
    # rna_output_column is a dictionary mapping all sequence ids to alignment columns of their particular alignment
    # all_accession_list is a dictionary of lists, telling the species in each row of the alignment
    # all_description_list is a dictionary of lists
    # unit_id_to_sequence_id maps unit id to sequence id

    #print("RNA output column: " + str(rna_output_column) + "\n")
    #print("protein output column: " + str(protein_output_column) + "\n")

    #print unit_id_to_sequence_id
    #print("accession_list: " + str(accession_list) + "\n") #for this test, should be 363 total
    #print("description_list: " + str(description_list) + "\n")

    ####residue_output_column[str(structure + "|" + position)] = residue_columns
    temp_var = ""

    # loop over each set of input residues
    for j, record in enumerate(sequence_ids):           # record is one list of unit ids

        accession_list = all_accession_list[j]          # list of species for the current alignment
        description_list = all_description_list[j]  # list of descriptions for the current alignment

        # possibly very long filename based on sequence ids
        possible_name_string = "_"
        final_column_file_name = possible_name_string.join(record).replace("|","-")

        # for now, just number the output files, force the user to re-name them
        final_column_file_name = os.path.join(OUTPUT_DIRECTORY,"output_file_%d.txt" % j)

        # make a header line called position_order
        position_order = ""

        counter = -1
        for i in record:
            counter += 1

        # loop over sequence ids
        for i, sequence_id in enumerate(record):
            if i < counter:
                position_order = position_order + str(sequence_id + " | ")
            elif i == counter:
                position_order = position_order + str(sequence_id)

        header = ("NCBI ID  | %s | Description" % position_order)

        print("Output file header")
        print(header)

        ordering = []

        print("Sequence ids that are mapped to columns of their corresponding alignment")
        print(rna_output_column.keys())

        # write output file
        with open(final_column_file_name, "w") as f:

            f.write(header+"\n")

            ordering = str(position_order).split(" | ")

            for i, species in enumerate(accession_list): #make the following into a gen fnt

                # if only RNA nucleotides are input ...
                # accumulate RNA sequence across each column
                temp = ""
                for sequence_id in record:
                    temp += rna_output_column[sequence_id][i]

                # print sequence to the screen
                # print("> %s" % species)  # fasta format
                # print(temp)   # row i of the output column for this record
                print(temp + " " + species + description_list[i])

                use_old_output_code = False

                if use_old_output_code:
                    f.write(str(species).ljust(18, ' '))

                    # loop over sequence ids
                    for o in ordering:
                        print(o)

                        # map sequence id back to unit id, it seems
                        # new_key will be observed residue number, not what we want anymore
                        for key, value in unit_id_to_sequence_id.items():
                            if value == o:
                                print(o)
                                print(key)
                                print(value)
                                #import pdb; pdb.set_trace()
                                #print "if statement met"
                                struc, model, chain, mon, new_position = key.split("|")
                                new_key = key.replace(o, new_position)

                                #print new_position

                        #
                        new_key = o

                        if str(o) in protein_output_column:
                            print("in protein")
                            column = protein_output_column[o]
                            for j, monomer in enumerate(column):
                                if j == i:
                                    temp_var = monomer

                            f.write(str(temp_var + "         |        "))


                        if unit_id_to_sequence_id[new_key] in rna_output_column:
                            print("in rna")
                            column = rna_output_column[o]
                            #print column
                            for j, monomer in enumerate(column):
                                if j == i:
                                    temp_var = monomer

                            f.write(str(temp_var + "         |       "))

                    f.write(str(description_list[i]) + "\n")



def main():

    print_instructions = False

    if print_instructions:
        print("Please input IDs in one of the following formats:\n     Author assigned chains and positions for a PDB structure: 4y4o|1|1a|U|1240\n     BGSU RNA Internal Loop IDs: IL_6ZMI_184")
        print("All unit IDs should be typed so that a comma separates each ID.")
        print("If you would like multiple output files, please separate each set of IDs with a semicolon as shown: ")
        print("4v4q|1|A|U|1240,4v4q|1|F|MET|115,4v4q|1|A|A|1239;4y4o|1|FA|U|1240,4y4o|1|LA|ARG|32,4y4o|1|LA|ASN|109\n\n\n")
        print("For more information on how to identify author naming styles, please type \"continue\": ")
        more_info = "" #"continue"
        if more_info == "continue":
            print("Author chain information may be identified on a PDB structure as [Auth]. If using the BGSU RNA site, all \
        ids in the above style are Author formatted IDs.")

    # input_ids is a list of lists of strings of unit ids or unit ids connected by :
    input_ids = [["IL_4Y4O_206"]]
    input_ids = [["IL_4WFL_001"],["HL_4WFL_001"]]  # tripled sheared IL from bacterial SRP Alu domain, RF01854
    input_ids = "4v4q|1|A|U|1240,4v4q|1|F|MET|115,4v4q|1|A|A|1239;4y4o|1|FA|U|1240,4y4o|1|LA|ARG|32,4y4o|1|LA|ASN|109"
    input_ids = [["5J7L|1|CA|A|52"]]
    input_ids = [["IL_4WFL_001"],["HL_4WFL_001"]]  # tripled sheared IL from bacterial SRP Alu domain, RF01854
    input_ids = [["IL_4V9F_106"]]


    # for user convenience, turn input string into list of lists
    if not isinstance(input_ids,list):
        new_input_ids = []
        for input_list in input_ids.split(";"):
            new_input_ids.append(input_list.split(","))
        input_ids = new_input_ids

    print("input_ids: %s" % input_ids)

    # these are old and will probably not be needed
    UPDATE = False # alt: True ; true downloads all files again, reforms alignments manually
    OUTPUT_TYPE = "coupled" #alt: 'concatenated' - no longer used
    phased_alignment = []
    PRESENCE_OF_PROTEINS = detect_presence_of_proteins(input_ids)
    # What happens if one list of unit IDs has amino acids but a different one does not?

    # convert HL, IL, J3 ids for RNA loops to one or more unit id strings with ranges
    input_unit_ids, all_unit_ids = convert_loop_ids_to_unit_ids(input_ids)
    print("input_unit_ids: %s" % input_unit_ids)

    # this may not be needed anymore
    all_PDB_ids, all_chain_ids = extract_PDB_and_chains(all_unit_ids)
    print("all_unit_ids: %s" % all_unit_ids)

    # get one PDB id for each list of unit ids
    ANCHOR_PDB_STRUCTURES = identify_anchor_PDB_structures(input_unit_ids)

    if PRESENCE_OF_PROTEINS and sys.version_info[0] < 3:
        # use pdbx, which is only working with Python 2.7 right now
        print("mapping unit ids to sequence ids for these PDB files")
        print(list(set(ANCHOR_PDB_STRUCTURES)))
        unit_id_to_sequence_id, sequence_id_to_unit_id = map_unit_ids_to_sequence_positions(ANCHOR_PDB_STRUCTURES)
    else:
        # use API, which only works for RNA chains as of 2022-03-30
        unit_id_to_sequence_id, sequence_id_to_unit_id = map_unit_ids_to_sequence_positions_chains(all_chain_ids)

    # load seqres file from PDB ... may be able to be omitted
    #seqres_lines = read_pdb_seqres()

    #protein_anchor_sequences, rna_anchor_sequences = process_anchor_PDB_sequences(seqres_lines)

    # get one Rfam family and chain, one Pfam family and chain, for each list of unit ids
    #Rfam_families, RNA_chains, Pfam_families, protein_chains = determine_families(input_unit_ids)

    # map unit ids strings (and ranges) to the Rfam/Pfam family and chain range in that family
    unit_id_string_to_family_and_chain, families, chains = determine_xfam_families(input_unit_ids)
    print("Mapping of unit id strings to xfam family and chain")
    print(unit_id_string_to_family_and_chain)
    print("All Pfam/Rfam families: %s" % families)
    print("All Pfam/Rfam chains: %s" % chains)

    # load all required alignments
    alignment = {}
    for family in families:
        accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list = load_alignment(family)

        alignment[family] = {}
        alignment[family]["accession"] = accession_list
        alignment[family]["NCBI_taxid"] = NCBI_taxid_list
        alignment[family]["description"] = description_list
        alignment[family]["sequence"] = sequence_list

    # map all sequence ids from the required chains to columns of their alignment
    sequence_id_to_column_mapping = {}

    for unit_id_list in input_unit_ids:
        for unit_id_string in unit_id_list:
            # identify Pfam/Rfam family and PDB experimental sequence chain in that alignment
            family,chain = unit_id_string_to_family_and_chain[unit_id_string]

            # get a dictionary that maps sequence ids to column for the given chain and alignment
            if chain and family and not chain in sequence_id_to_column_mapping:
                sequence_id_to_column_mapping[chain] = map_sequence_ids_to_columns(alignment[family],chain)

                print("Some of the sequence id to column mappings for %s:" % chain)
                sequence_id_keys = sequence_id_to_column_mapping[chain].keys()
                sequence_id_keys = sorted(sequence_id_keys, key=lambda x:int(x.split("|")[4]))
                for i in range(0,min(10,len(sequence_id_keys))):
                    print(sequence_id_keys[i],sequence_id_to_column_mapping[chain][sequence_id_keys[i]])

    # convert unit id list of lists to sequence id list of lists
    sequence_ids = look_up_sequence_ids(input_unit_ids, unit_id_to_sequence_id)
    print("Sequence ids list of lists")
    print(sequence_ids)

    # map unit_id_string to the rows of the alignment to use, in order
    # todo: if one unit_id_list has more than one family, intersect on taxonomy id
    unit_id_string_to_alignment_row = {}
    for unit_id_list in input_unit_ids:
        all_families = set([])
        for unit_id_string in unit_id_list:
            family,chain = unit_id_string_to_family_and_chain[unit_id_string]
            all_families.add(family)

        if len(all_families) == 1:
            for unit_id_string in unit_id_list:
                unit_id_string_to_alignment_row[unit_id_string] = range(0,len(alignment[family]["sequence"]))
        else:
            print("Need to write code to intersect on NCBI taxonomy ID")
            print("Returning one row just to have something, but you can't trust it")
            print(unit_id_list)
            for unit_id_string in unit_id_list:
                unit_id_string_to_alignment_row[unit_id_string] = [1]

    # retrieve columns of the alignment for each unit id string
    unit_id_string_to_sequence_extract = {}
    for unit_id_list in input_unit_ids:
        for unit_id_string in unit_id_list:
            family,chain = unit_id_string_to_family_and_chain[unit_id_string]
            sequences = alignment[family]["sequence"]

            unit_ids = unit_id_string.split(":")

            unit_id = unit_ids[0]
            sequence_id = unit_id_to_sequence_id[unit_id]
            column = sequence_id_to_column_mapping[chain][sequence_id]

            if len(unit_ids) == 1:                   # just one unit id
                unit_id_string_to_sequence_extract[unit_id_string] = [sequence[column] for sequence in sequences]
            else:
                unit_id = unit_ids[-1]
                sequence_id = unit_id_to_sequence_id[unit_id]
                end_column = sequence_id_to_column_mapping[chain][sequence_id]
                unit_id_string_to_sequence_extract[unit_id_string] = [sequence[column:(end_column+1)] for sequence in sequences]

                for unit_id in unit_id_string.split(":"):
                    sequence_id = unit_id_to_sequence_id[unit_id]

                    ids = sequence_id.split(":")        # could be a range

    print("Extracted sequences for these unit id strings:")
    print(unit_id_string_to_sequence_extract.keys())

    # produce an output file for each unit_id_list in input_unit_ids
    for i,unit_id_list in enumerate(input_unit_ids):
        # produce four separate header lines
        header1 = ""
        header2 = "\t\t"
        header3 = "\t\t"
        header4 = "\t\t"
        for unit_id_string in unit_id_list:
            header1 += unit_id_string + "\t"
            family,chain = unit_id_string_to_family_and_chain[unit_id_string]
            header2 += "\t" + family
            header3 += "\t" + chain
            columns = []
            unit_ids = unit_id_string.split(":")
            sequence_id = unit_id_to_sequence_id[unit_ids[0]]
            first_column = sequence_id_to_column_mapping[chain][sequence_id]

            for unit_id in unit_ids:
                sequence_id = unit_id_to_sequence_id[unit_id]
                column = sequence_id_to_column_mapping[chain][sequence_id]
                columns.append(str(column-first_column+1))
            header4 += "\t" + ":".join(columns)

        header1 += "NCBI_taxid\tAccession\tDescription"
        output = header1 + "\n" + header2 + "\n" + header3 + "\n" + header4 + "\n"

        # first columns of output are from the first family and chain
        # but if multiple families are involved, then one might want
        # columns for each of the families.  That will be complicated!
        first_string = unit_id_list[0]
        family,chain = unit_id_string_to_family_and_chain[first_string]
        species = alignment[family]["accession"]
        NCBI_taxid = alignment[family]["NCBI_taxid"]
        description = alignment[family]["description"]

        for j in unit_id_string_to_alignment_row[first_string]:
            for unit_id_string in unit_id_list:
                output += unit_id_string_to_sequence_extract[unit_id_string][j] + "\t"

            output += str(NCBI_taxid[j]) + "\t"
            output += str(species[j]) + "\t"
            output += str(description[j])

            output += "\n"

        print("First 1000 characters of the output")
        print(output[0:(min(1000,len(output)))])

        # write to an output file
        output_filename = os.path.join(OUTPUT_DIRECTORY,"%s_%s.txt" % (input_ids[i][0],family))
        output_filename = output_filename.replace("|","-")
        print("Writing output to %s" % output_filename)

        with open(output_filename,"w") as fid:
            fid.write(output)

    # this is all the code we need; code below is old and probably not needed anymore
    return

    """
    # get a dictionary that maps sequence ids to column for the given RNA chain and alignment
    chain = RNA_chains[j]
    family = Rfam_families[j]
    if chain and family and not chain in sequence_to_column_mapping:
        sequence_to_column_mapping[chain] = map_sequence_ids_to_columns(alignment[family],chain)
        print("sequence id to column mapping for %s" % chain)
        print(sequence_to_column_mapping[chain])

    # get a dictionary that maps sequence ids to column for the given protein chain and alignment
    chain = protein_chains[j]
    family = Pfam_families[j]
    if chain and family and not chain in sequence_to_column_mapping:
        sequence_to_column_mapping[chain] = map_sequence_ids_to_columns(alignment[family],chain)
        print("sequence id to column mapping for %s" % chain)
        print(sequence_to_column_mapping[chain])
    """


    # process RNA-only queries
    if PRESENCE_OF_PROTEINS == False:

        #rna_anchor_sequences = obtain_anchor_PDB_sequences(input_unit_ids, seqres_lines)

        """this if update segment will be deleted after full family expansion is completed and
        identification of pdb linked sequences in the family is solved"""
        if UPDATE:
            cmalign_sequences(RF)
            merge_alignments_RNA(RF)

        #pdb_anchors_aligned_rna = append_cmaligned_sequences(rna_anchor_sequences, RF)

    # process queries to create coupled RNA-protein alignments
    if PRESENCE_OF_PROTEINS == True:
        for i in range(len(Pfam_families)):
            RF = Rfam_families[i]
            PF = Pfam_families[i]

            open_family_files(RF, PF)
            protein_anchor_sequences, protein_struc_header, rna_anchor_sequences = obtain_anchor_PDB_sequences(input_unit_ids, seqres_lines)

            obtain_pfam_sequences(PF)
            write_file_protein_sequences_with_anchors(protein_struc_header, PF)
            if UPDATE:
                hmmbuild(PF)
                hmmalign_sequences(PF)
                merge_alignments_protein(PF)
            pdb_anchors_aligned_protein = append_hmmaligned_sequences(protein_anchor_sequences, PF)

            write_file_rnacentral_sequences_with_anchors(rna_anchor_sequences, RF)
            if UPDATE:
                cmalign_sequences(RF)
                merge_alignments_RNA(RF)
            pdb_anchors_aligned_rna = append_cmaligned_sequences(rna_anchor_sequences, RF)

            create_Phased_id()
            create_Phased_sequences()
            if OUTPUT_TYPE == 'concatenated':
                forming_concatenated_stockholm_file(PF, RF)
            elif OUTPUT_TYPE == 'coupled':
                forming_individual_stockholm_files(PF, RF, pdb_anchors_aligned_protein, pdb_anchors_aligned_rna)
            else:
                print("There is an error in the output method assignment")


    rna_output_column, protein_output_column, accession_list, description_list = identify_alignment_columns(sequence_ids, Rfam_families, RNA_chains,Pfam_families,protein_chains)

    print("Output columns determined")
    # print(rna_output_column)

    display_columns_file(sequence_ids, rna_output_column, protein_output_column, accession_list, description_list, unit_id_to_sequence_id)



if __name__ == '__main__':
    main()


### the following code is what we run in pfam to obtain the pfam seed

# select pf.pfamseq_id, pf.ncbi_taxid
# from pfamA_reg_seed p
# Join pfamseq pf
# on p.pfamseq_acc = pf.pfamseq_acc
# where p.pfamA_acc = "PF00177"

### the following code is what we run in rfam to obtain the rfam seed information

# select rfs.rfamseq_acc, rfs.ncbi_id
# from seed_region sr
# Join rfamseq rfs
# on sr.rfamseq_acc = rfs.rfamseq_acc
# where sr.rfam_acc = "RF00177"