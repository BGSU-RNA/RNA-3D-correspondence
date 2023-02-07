"""
This program takes loop id or unit id as input and maps units to units in other PDB files
It uses alignments of PDB chains to make the mapping.

For debugging, this works:
return [{"text": "%s" % chain_to_equivalence_class_and_position}]

"""

from asyncore import loop
import csv
import gzip
import os.path
import os
import re
#from types import new_class
import requests
import sys
import textwrap
from datetime import datetime

from Bio import AlignIO

from RPFAM import *

from RPFAM_configuration import BASE_DIRECTORY

# manage imports according to Python version present
if sys.version_info[0] < 3:
    # for experimental sequence mapping to resolved sequences
    import pdbx    # currently only really working for Python 2.7
    from pdbx.reader.PdbxReader import PdbxReader


def make_chain_sets(pdb_chain_to_triple_sets):
    """
    Make a list of mappings from pdb_chain to chain triple
    """

    chain_sets = []

    pdb_chains = list(pdb_chain_to_triple_sets.keys())
    M = len(pdb_chains)
    # track where we are in each list
    current_position = [0 for pdb_chain in pdb_chains]

    while True:
        chain_set = {}
        for j in range(0,M):
            pdb_chain = pdb_chains[j]
            if current_position[j] < len(pdb_chain_to_triple_sets[pdb_chain]):
                chain_set[pdb_chain] = pdb_chain_to_triple_sets[pdb_chain][current_position[j]]
            else:
                chain_set[pdb_chain] = (None,None,None)

        chain_sets.append(chain_set)

        # for now, just return one set per PDB file
        # later, change the current_position to cover all combinations of chains in this PDB
        break

    return chain_sets

def map_across_species(nt_lists,scope,resolution,depth):
    """
    nt_lists is a list of lists of strings
    scope is a string like EC, Rfam, molecule
    resolution is a number like 1.5, 2.0, 2.5
    depth is an integer telling how far to go into each equivalence class
    """

    print("Time %s." % datetime.now())

    # convert HL, IL, J3 ids for RNA loops to one or more unit id strings with ranges
    input_unit_ids, all_unit_ids = convert_loop_ids_to_unit_ids(nt_lists)

    if not all_unit_ids:
        raise Exception("Not able to look up unit ids for this request; maybe the PDB id is too new or the loop does not exist.")

    print("input_unit_ids: %s" % input_unit_ids)
    #print("all_unit_ids: %s" % all_unit_ids)

    # load mappings from rfam id to pdb chains
    rfam_family_to_chains, rfam_chain_to_family, pdb_chain_to_rfam_family = read_rfam_to_pdb_mappings()

    # get the current list of equivalence classes at specified resolution and parse it for the chains in them
    print('Time %s. Loading EC at specified resolution' % datetime.now())
    ife_to_equivalence_class, equivalence_class_to_ifes, chain_to_equivalence_class_and_position = get_equivalence_classes(resolution)

    # get the current list of equivalence classes with no resolution cutoff in case given units are not at specified resolution
    print('Time %s. Loading EC at all resolution' % datetime.now())
    ife_to_equivalence_class_all, equivalence_class_to_ifes_all, chain_to_equivalence_class_and_position_all = get_equivalence_classes('all')

    # make a place to store alignments
    alignment = {}
    chain_sequence_id_to_column = {}

    # make a place to store outputs over all inputs
    result_list = []

    # make a place to store the following mappings
    chain_sequence_id_to_unit_id = {}
    chain_sequence_id_to_column = {}
    column_to_sequence_id = {}

    # loop over requested unit id sets
    for unit_id_list in input_unit_ids:
        # make a place to store outputs for this one input
        result = {}
        result["mappings"] = []

        print('Time %s. Current unit_id_list' % datetime.now())
        print(unit_id_list)

        # make a list of all unit ids to map
        unit_ids_to_map = []
        for unit_id_string in unit_id_list:
            unit_ids_to_map += unit_id_string.replace(",",":").split(":")  # just make a list of unit ids
        N = len(unit_ids_to_map)

        print("Time %s. Mapping these %d unit ids:" % (datetime.now(),N))
        print(unit_ids_to_map)

        all_PDB_ids, all_pdb_chains = extract_PDB_and_chains(unit_ids_to_map)
        #print("all_PDB_ids: %s" % all_PDB_ids)
        print("Time %s. All PDB chains: %s" % (datetime.now(),all_pdb_chains))

        # use local file or API to get unit ids and sequence ids and the mapping between them
        for pdb_chain in set(all_pdb_chains):
            if not pdb_chain in chain_sequence_id_to_unit_id:
                print("Time %s. Mapping sequence ids to unit ids for %s.  " % (datetime.now(),pdb_chain))
                unit_id_to_sequence_id, sequence_id_to_unit_id = map_unit_ids_to_sequence_positions_chains([pdb_chain])
                # store according to the chain; better to have many small dictionaries
                chain_sequence_id_to_unit_id[pdb_chain] = sequence_id_to_unit_id

        # map unit ids to the Rfam family and chain range in that family
        unit_id_to_rfam_family_and_chain = map_unit_id_to_rfam([unit_ids_to_map], rfam_family_to_chains)

        # collect mapping information for each unit id
        map_out = []
        for unit_id,rfam_tuple in unit_id_to_rfam_family_and_chain.items():
            pdb_chain = "|".join(unit_id.split("|")[0:3])
            rfam_family,rfam_chain = rfam_tuple

            mapping = {}
            mapping["unit_id"] = unit_id
            mapping["pdb_chain"] = pdb_chain
            mapping["rfam_chain"] = rfam_chain
            mapping["rfam_family"] = rfam_family
            mapping['rfam_families'] = [rfam_family]
            mapping['alignment_name'] = rfam_family

            # read a joint alignment if requested, to work with SSU or LSU instead of just the given Rfam family
            # tRNA and 5S rRNA are already cross domain alignments
            if scope == 'molecule':
                for molecule in joint_alignments.keys():
                    if rfam_family in joint_alignments[molecule]:
                        # replace original family name with name of joint alignment family like SSU or LSU
                        mapping['alignment_name'] = molecule
                        mapping['rfam_families'] = joint_alignments[molecule]
                        print('Time %s. Mapping %s across joint alignment of %s involving %s' % (datetime.now(),mapping['unit_id'],molecule,joint_alignments[molecule]))

            map_out.append(mapping)

            print("Time %s. %s maps to %s in family %s" % (datetime.now(),pdb_chain,rfam_chain,rfam_family))

        # load required alignments or joint alignments
        for i in range(0,N):
            alignment_name = map_out[i]['alignment_name']
            if not alignment_name in alignment.keys():
                filename = os.path.join(RFAM_ALIGNMENT_DIRECTORY,"%s_PDB_chains.sto" % alignment_name)
                accession_list, NCBI_taxid_list, domain_list, description_list, sequence_list, codon_list = load_alignment(alignment_name,filename)

                alignment[alignment_name] = {}
                alignment[alignment_name]["accession"] = accession_list
                alignment[alignment_name]["NCBI_taxid"] = NCBI_taxid_list
                alignment[alignment_name]["domain"] = domain_list
                alignment[alignment_name]["description"] = description_list
                alignment[alignment_name]["sequence"] = sequence_list
                alignment[alignment_name]["codon"] = codon_list

        # for the chain of each given unit id, map sequence ids to columns
        for i in range(0,N):
            rfam_chain = map_out[i]['rfam_chain']
            alignment_name = map_out[i]['alignment_name']
            if not rfam_chain in chain_sequence_id_to_column:
                print('Time %s. Mapping sequence ids to columns for chain %s in %s' % (datetime.now(),rfam_chain,alignment_name))
                start_index = int(rfam_chain.split("_")[2])
                chain_sequence_id_to_column[rfam_chain],column_to_sequence_id[rfam_chain] = map_sequence_ids_to_columns(alignment[alignment_name],rfam_chain,start_index)

        # map given unit ids to sequence ids to columns of the relevant alignment
        for i in range(0,N):
            unit_id = map_out[i]['unit_id']
            rfam_chain = map_out[i]['rfam_chain']
            sequence_id = unit_id_to_sequence_id.get(unit_id,'')
            column = chain_sequence_id_to_column[rfam_chain].get(sequence_id,None)
            map_out[i]['column'] = column

            if column == None:
                print('Not able to map given unit id %s to a column, so something will go wrong' % unit_id)

        # for each chain in the unit ids to map, list all pdb chains in the corresponding Rfam family
        pdb_chain_to_alignable_pdb_chains = {}
        pdb_chain_to_rfam_chain = {}
        pdb_chains_ordered = []
        for i in range(0,N):
            pdb_chain = map_out[i]['pdb_chain']

            if not pdb_chain in pdb_chains_ordered:
                pdb_chains_ordered.append(pdb_chain)

            if not pdb_chain in pdb_chain_to_alignable_pdb_chains.keys():
                pdb_chain_to_alignable_pdb_chains[pdb_chain] = set([])
                for rfam_family in map_out[i]['rfam_families']:
                    for rfam_chain in rfam_family_to_chains[rfam_family]:
                        fields = rfam_chain.split("_")
                        alignable_pdb_chain = fields[0] + "|1|" + fields[1]
                        pdb_chain_to_alignable_pdb_chains[pdb_chain].add(alignable_pdb_chain)
                        pdb_chain_to_rfam_chain[alignable_pdb_chain] = rfam_chain

        if len(pdb_chain_to_alignable_pdb_chains) == 0:
            raise Exception("Did not find any PDB chains to map to")

        # sort the list of alignable pdb chains by position in the EC, leaving out those in no EC at this resolution
        # pdb_chain is a chain that given unit ids come from
        # a triple set is an alignable chain, equivalence class, and position of the chain in the EC
        pdb_chain_to_triple_sets = {}
        pdb_to_pdb_chain_to_triple_sets = {}
        for pdb_chain, alignable_pdb_chains in pdb_chain_to_alignable_pdb_chains.items():
            pdb_chain_to_triple_sets[pdb_chain] = []

            for alignable_pdb_chain in alignable_pdb_chains:
                equivalence_class,position = chain_to_equivalence_class_and_position.get(alignable_pdb_chain,(None,None))
                if equivalence_class:
                    if not scope == 'EC'\
                       or equivalence_class.split("_")[2] == chain_to_equivalence_class_and_position_all.get(pdb_chain,'__')[0].split("_")[2]:
                        pdb_chain_to_triple_sets[pdb_chain].append((alignable_pdb_chain,equivalence_class,position))

                        pdb = alignable_pdb_chain.split("|")[0]
                        if not pdb in pdb_to_pdb_chain_to_triple_sets.keys():
                            pdb_to_pdb_chain_to_triple_sets[pdb] = {}

                        if not pdb_chain in pdb_to_pdb_chain_to_triple_sets[pdb].keys():
                            pdb_to_pdb_chain_to_triple_sets[pdb][pdb_chain] = []

                        pdb_to_pdb_chain_to_triple_sets[pdb][pdb_chain].append((alignable_pdb_chain,equivalence_class,position))

            # let's try this
            #return [{"text": "%s" % pdb_chain_to_triple_sets[pdb_chain]}]


            pdb_chain_to_triple_sets[pdb_chain] = sorted(pdb_chain_to_triple_sets[pdb_chain], key = lambda x : (x[2],x[0]))
            #print(pdb_chain_to_triple_sets[pdb_chain])

            # sort each list by position, or recognize that a list is empty
            for pdb in pdb_to_pdb_chain_to_triple_sets.keys():
                if pdb_chain in pdb_to_pdb_chain_to_triple_sets[pdb].keys():
                    pdb_to_pdb_chain_to_triple_sets[pdb][pdb_chain] = sorted(pdb_to_pdb_chain_to_triple_sets[pdb][pdb_chain], key = lambda x : (x[2],x[0]))
                else:
                    pdb_to_pdb_chain_to_triple_sets[pdb][pdb_chain] = []

        if len(pdb_to_pdb_chain_to_triple_sets) == 0:
            raise Exception("Did not find any PDB files to map to")

        # loop over distinct PDB files to which some alignment can be made
        for pdb in pdb_to_pdb_chain_to_triple_sets.keys():
            print('Time %s. Mapping given units to chains in %s' % (datetime.now(),pdb))

            print("pdb, pdb_to_pdb_chain_to_triple_sets",pdb,pdb_to_pdb_chain_to_triple_sets[pdb])

            # assemble sets of alignable chains from the same PDB file
            chain_sets = make_chain_sets(pdb_to_pdb_chain_to_triple_sets[pdb])

            # Map columns back to unit ids in PDB files and in biological assemblies inside PDB files
            # When one PDB file has multiple tRNAs, it may not be clear which one to map to
            # When one PDB file has multiple ribosomes, it may not be clear which SSU goes which which LSU, 5S
            # Keep the original unit ids to make it possible to identify ones that are very far different

            # track depth in the list somehow; this is not good enough
            d = 0

            # loop over sets of alignable chains from this PDB file
            for chain_set in chain_sets:

                if d >= depth:
                    break

                mapped_unit_ids = []
                missing_unit_id = False

                print('Time %s. Working with chain_set %s' % (datetime.now(),chain_set))

                # loop over unit ids to map
                for i in range(0,N):
                    pdb_chain = map_out[i]['pdb_chain']
                    alignment_name = map_out[i]['alignment_name']
                    column = map_out[i]['column']

                    alignable_pdb_chain = chain_set.get(pdb_chain,[''])[0]

                    if alignable_pdb_chain:
                        rfam_chain = pdb_chain_to_rfam_chain[alignable_pdb_chain]

                        # use local file or API to get unit ids and sequence ids and the mapping between them
                        if not alignable_pdb_chain in chain_sequence_id_to_unit_id:
                            #print("Time %s. Mapping sequence ids to unit ids for %s.  " % (datetime.now(),alignable_pdb_chain))
                            unit_id_to_sequence_id, sequence_id_to_unit_id = map_unit_ids_to_sequence_positions_chains([alignable_pdb_chain])
                            chain_sequence_id_to_unit_id[alignable_pdb_chain] = sequence_id_to_unit_id

                        # map sequence ids to columns for the entire current rfam chain
                        if not rfam_chain in chain_sequence_id_to_column:
                            print('Time %s. Mapping sequence id to column for %s in %s' % (datetime.now(),rfam_chain,rfam_family))
                            start_index = int(rfam_chain.split("_")[2])
                            chain_sequence_id_to_column[rfam_chain],column_to_sequence_id[rfam_chain] = map_sequence_ids_to_columns(alignment[alignment_name],rfam_chain,start_index)

                            if False and chain_sequence_id_to_column[rfam_chain]:
                                print("Some of the sequence id to column mappings for %s:" % rfam_chain)
                                sequence_id_keys = chain_sequence_id_to_column[rfam_chain].keys()
                                sequence_id_keys = sorted(sequence_id_keys, key=lambda x:int(x.split("|")[4]))
                                for i in range(0,min(5,len(sequence_id_keys))):
                                    print(sequence_id_keys[i],chain_sequence_id_to_column[rfam_chain][sequence_id_keys[i]])

                        # map column back to sequence ids and to unit ids for target chains
                        sequence_id = column_to_sequence_id[rfam_chain].get(column,'')

                        if sequence_id:
                            unit_id = chain_sequence_id_to_unit_id[alignable_pdb_chain].get(sequence_id,'Not resolved')
                        else:
                            unit_id = "No alignment"

                        if unit_id == 'NULL':
                            unit_id = 'Not resolved'

                        mapped_unit_ids.append(unit_id)
                        if unit_id == '':
                            missing_unit_id = True

                    else:
                        mapped_unit_ids.append('No map')
                        missing_unit_id = True

                missing_unit_id = False  # suppress these warnings


                if missing_unit_id:
                    print('Missing at least one unit id for %s in %s in %s' % (ife,equivalence_class,rfam_chain_to_family[rfam_chain]))
                else:
                    # every line of the result is a dictionary with certain fields
                    one_result = {}

                    one_result["unit_id_list"] = mapped_unit_ids
                    one_result['pdb'] = pdb  # useful in case there are no mapped units at all

                    one_result['rfam_EC_chain'] = []

                    # loop over chains in the given unit ids to get the chains they map to
                    for pdb_chain in pdb_chains_ordered:
                        alignable_pdb_chain = chain_set.get(pdb_chain,[''])[0]
                        if alignable_pdb_chain:
                            one_result['rfam_EC_chain'].append((pdb_chain_to_rfam_family[alignable_pdb_chain],chain_to_equivalence_class_and_position[alignable_pdb_chain][0],alignable_pdb_chain))
                        else:
                            one_result['rfam_EC_chain'].append(('No Rfam', 'No EC', 'No chain'))

                    result["mappings"].append(one_result)

                    d += 1

        # include the original unit ids in the same format as each mapping
        one_result = {}
        one_result["unit_id_list"] = unit_ids_to_map
        one_result['pdb'] = unit_ids_to_map[0].split("|")[0]

        one_result['rfam_EC_chain'] = []

        for pdb_chain in pdb_chains_ordered:
            one_result['rfam_EC_chain'].append((pdb_chain_to_rfam_family[pdb_chain],chain_to_equivalence_class_and_position_all[pdb_chain][0],pdb_chain))

        result["query"] = one_result

        # produce text output
        output_list = []
        c = 0
        for one_result in result["mappings"]:
            if c == 0:
                # make a header once
                header = ''

                j = 1
                for rfam,EC,pdb_chain in one_result['rfam_EC_chain']:
                    header += 'Rfam family %d\t' % j
                    header += 'EC %d\t' % j
                    header += "Mapped chain %d\t" % j
                    j += 1

                j = 1
                for unit_id in one_result['unit_id_list']:
                    header += 'Position %d\t' % j
                    j += 1

                header += 'Viewer\t'
                c = 1

            output = ''

            for rfam,EC,pdb_chain in one_result['rfam_EC_chain']:
                output += '%s\t%s\t%s\t' % (rfam,EC,pdb_chain)

            view_list = []
            for unit_id in one_result['unit_id_list']:
                output += '%s\t' % unit_id
                if "|" in unit_id:
                    view_list.append(unit_id)

            output += 'http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s' % ",".join(view_list)

            output_list.append(output)

        output_list = [header] + output_list

        result["text"] = "\n".join(output_list)

        # append all the results from this set of unit ids
        result_list.append(result)

    return result_list

    
def main():

    nt_lists = [['7K00|1|A|G|1405:7K00|1|A|U|1406:7K00|1|A|5MC|1407:7K00|1|A|A|1408:7K00|1|A|C|1409', '7K00|1|A|G|1491:7K00|1|A|A|1492:7K00|1|A|A|1493:7K00|1|A|G|1494:7K00|1|A|U|1495:7K00|1|A|C|1496']]
    nt_lists = [["IL_5J7L_035"]]  # bacterial kink turn
    nt_lists = [["IL_4V9F_106"]]
    nt_lists = [['5J7L|1|AA|A|253','5J7L|1|AA|U|273']]
    nt_lists = [['5J7L|1|AA|A|253,5J7L|1|AA|U|273']]
    nt_lists = [['5J7L|1|AA|A|253']]
    nt_lists = [["HL_4TNA_002"]]                      # tRNA anti-codon loop, slow!
    nt_lists = [["IL_5TBW_013"]]                      # each strand is in a different chain, fails
    nt_lists = [["IL_4V9F_089,4V9F|1|0|U|2554"]]      # minor groove with intercalated U
    nt_lists = [["IL_4V9F_007"]]
    nt_lists = [["IL_4V88_434"]]                      # eukaryotic kink turn
    nt_lists = [["3J7P|1|S2|G|547,3J7P|1|S2|G|535"]]  # GG cWW basepair
    nt_lists = [["IL_7RQB_004"]]                      # IL missing nucleotide 101, getting 3 strands
    nt_lists = [["IL_4V9F_007"]]  # archaeal LSU
    nt_lists = [["IL_5J7L_060"]]  # bacterial SSU decoding loop
    nt_lists = [["IL_4V9F_007"]]  # archaeal LSU
    nt_lists = [["IL_5J7L_024"]]  # bacterial SSU
    nt_lists = [["IL_6ZMI_005"]]  # eukaryotic LSU including 5.8S and long chain
    nt_lists = [["IL_6ZMI_009"]]  # eukaryotic LSU kink turn, not well aligned to bacteria
    nt_lists = [["IL_5J7L_014"]]  # bacterial SSU kink turn
    nt_lists = [["IL_6ZMI_178"]]  # eukaryotic 5S symmetric platform motif
    nt_lists = [["5J7L|1|DA|A|1858,5J7L|1|DA|G|1884"]]  # AG tHS in E. coli but GA tSH in T. th.
    nt_lists = [["7K00|1|a|G|1992,7K00|1|a|A|1664"]]
    nt_lists = [["IL_6ZMI_024"]]  # eukaryotic LSU kink turn
    nt_lists = [["6ZMI|1|L5|G|16,6ZMI|1|L8|C|141,6ZMI|1|L8|C|26,6ZMI|1|L5|G|348,6ZMI|1|S2|U|120,6ZMI|1|S2|U|344,6ZMI|1|L7|C|95,6ZMI|1|L7|G|81"]] # four chains from human ribosome
    nt_lists = [["HL_7K00_033"]]                      # tRNA anti-codon loop, slow!  Also fails.
    nt_lists = [["IL_6ZMI_019"]]  # eukaryotic LSU + 5.8S big IL

    scope = 'Rfam'     # same Rfam family as each individual unit
    scope = 'EC'       # equivalence class, for each unit
    scope = 'molecule' # same molecule like SSU or LSU where available

    resolution = '2.5A'   # resolution cutoff for equivalence classes
    resolution = '3.0A'   # resolution cutoff for equivalence classes

    if not resolution in ['1.5A','2.0A','2.5A','3.0A','3.5A','4.0A','20.0A','all']:
        resolution = '3.0A'

    depth = 3        # how many IFEs from the equivalence class to map to

    result_list = map_across_species(nt_lists,scope,resolution,depth)

    for i, result in enumerate(result_list):
        #print(result["text"])

        with open('results/map_across_species_%d.txt' % i,write_mode) as f:
            f.write(result["text"])

    return

if __name__ == '__main__':
    main()
