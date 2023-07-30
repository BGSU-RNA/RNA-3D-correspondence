from __future__ import division
from collections import OrderedDict, defaultdict
from discrepancy import matrix_discrepancy, relative_discrepancy
import numpy as np
# import pandas as pd
#from ordering import optimalLeafOrder
import io
import itertools
import re
from ordering_similarity import treePenalizedPathLength
import math
from operator import itemgetter
import requests

EXCLUDE_LIST = ['No alignment', 'Not resolved', 'No chain']

def check_valid_single_chain_query(units):
    all_chains = [unit.split("|")[2] for unit in units]
    if len(set(all_chains)) == 1:
        return True
    else:
        return False

def get_disc(d, first, second):
    return d.get(second, {}).get(first, 0.0)

def get_nt_position_index(data):
    nt_position_index = {}
    for ife, nt_positions in data.iteritems():
        for idx, position in enumerate(nt_positions, 1):
            nt_position_index[position] = idx
    return nt_position_index

def filter_dict_by_list(data, filter_list):
    updated_dict = {}
    for k,v in data.iteritems():
        pdb = k.split("|")[0]
        if pdb in filter_list:
            updated_dict[k] = v
    return updated_dict

def create_all_nt_pairs(data):
    nt_pairs_dict = OrderedDict()
    for ife, nt_positions in data.iteritems():
        nt_pairs_dict[ife] = list(itertools.combinations(nt_positions, 2))
    return nt_pairs_dict

# def format_species_name(name_list):
#     formatted_names = []
#     for species_name in name_list:
#         name_components = species_name.split(" ")
#         if len(name_components) > 2:
#             formatted_species_name = " ".join(species_name.split(" ")[:2])
#             formatted_names.append(str(formatted_species_name))
#         else:
#             formatted_names.append(str(species_name))
#     return formatted_names

# def format_species_name(name_list):
#     return [ " ".join(name.split()[:2]) if len(name.split()) > 2 else name for name in name_list ]

def format_species_name(name):
    if (name is not None) and (len(name.split()) > 2):
        return " ".join(name.split()[:2])
    elif (name is None) or (name == "synthetic construct"):
        return " "
    else:
        return name

def get_name_count(data):
    data = [i for i in data if i.strip()]
    count_dict = {}
    for name in data:
        count_dict[name] = count_dict.get(name, 0) + 1
    count_values = [(k,v) for k,v in count_dict.iteritems() if k is not None]
    return sorted(count_values, key=lambda x: x[1], reverse=True)

def format_query_units(query_units):
    formatted_units = []
    for unit in query_units:
        unit = "|".join(unit.split("|")[3:])
        formatted_units.append(str(unit))
    return ", ".join(formatted_units)

def get_correspondence_dict(correspondence):
    correspondence_dict = {}
    for sublist in correspondence:
        ife = "|".join(sublist[0].split("|")[:3])
        correspondence_dict[ife] = sublist
    return correspondence_dict

def get_pdb_and_chain_from_ife(ife):
    pdb, _, chain = ife.split("|")
    return pdb, chain

def get_correspondence_across_species(param_dict):
    base_url = "http://rna.bgsu.edu/correspondence/map_across_species?id="
    suffix_url = "&format=json"
    # depth = param_dict.get('depth', 1)
    complete_url = base_url + str(param_dict['selection']) +  "&scope=" + str(param_dict['scope']) + "&depth=" + str(param_dict['depth']) + "&resolution=" + str(param_dict['resolution']) + suffix_url
    response = requests.get(complete_url)

    if response.status_code == 200:
        response = response.json()
        query_nts_list = response['query']['unit_id_list']
        query_nts_str = format_query_units(query_nts_list)
        rfam_accession = response['query']['rfam_EC_chain'][0][0]
        query_ife = response['query']['rfam_EC_chain'][0][2]
        pdb, chain = get_pdb_and_chain_from_ife(query_ife)

        query_info = {
            "pdb": str(pdb),
            "chain": str(chain),
            "resolution_limit": str(param_dict['resolution']),
            "rfam_accession": str(rfam_accession),
            "rfam_description": "",
            "query_nts_list": query_nts_list,
            "query_nts_str": query_nts_str,
            "query_ife": query_ife
        }

        correspondence_list = []
        equivalence_class_dict = {}
        for entry in response["mappings"]:
            correspondence_list.append(entry["unit_id_list"])
            equivalence_class_dict[entry["rfam_EC_chain"][0][2]] = entry["rfam_EC_chain"][0][1]

        correspondence_list = [sublist for sublist in correspondence_list if not any(x in sublist for x in EXCLUDE_LIST)]
        return query_info, equivalence_class_dict, correspondence_list
    else:
        return None, None, None

def generate_sequence_logo_data(data):
    sequence_dict = OrderedDict()
    nt_list = ['A', 'C', 'U', 'G']
    count_list = []
    total_seq_count = sum(n for seq, n in data)

    for idx, val in enumerate(data[0][0], 1):
        sequence_dict.setdefault(idx, {})['A'] = 0
        sequence_dict.setdefault(idx, {})['C'] = 0
        sequence_dict.setdefault(idx, {})['U'] = 0
        sequence_dict.setdefault(idx, {})['G'] = 0

    for elem in data:
        seq = elem[0]
        count = elem[1]
        for idx, val in enumerate(seq, 1):
            if val in nt_list:
                sequence_dict[idx][val] += count

    for k, v in sequence_dict.iteritems():
        position_list = []
        for nt in nt_list:
            if v[nt] != 0:
                position_count = int(v[nt])/int(total_seq_count)
                position_count = "%.1f" % position_count
                position_list.append(position_count)
            else:
                position_list.append(0)
        count_list.append(position_list)

    return count_list

def get_sequence_variability(query_units_string):
    api_url = "http://rna.bgsu.edu/correspondence/variability?id={query}&format=unique".format(query=query_units_string)
    response = requests.get(api_url)
    data = response.text
    lines = list(io.StringIO(data, newline=''))
    lines = lines[2:]
    
    count_dict = {}
    for line in lines:
        test = re.split(r'\t+', line)
        count = test[-1].strip()
        sequence = "".join(test[:-1])
        count_dict[sequence] = count

    count_dict = [(str(k),int(v)) for k,v in count_dict.iteritems() if "-" not in k]
    return sorted(count_dict, key=lambda x: x[1], reverse=True)
    
def format_correspondence_display(corr):

	display_str = ""
	for ife, correspondence in corr.iteritems():
		display_str += ",".join(correspondence) + "</br>"

	return display_str


def format_pairwise_correspondence_display(correspondence):

    display_str = ""
    for elem in correspondence:
        display_str += elem[0] + " aligns_to " + elem[1] + "</br>"

    return display_str


def check_missing_correspondence(data, query_length):

    entries_with_missing_data = [k for k, v in data.iteritems() if len(v) != query_length]

    updated_correspondence = {}
    for k, v in data.iteritems():
        if len(v) == int(query_length):
            updated_correspondence[k] = v

    return entries_with_missing_data, updated_correspondence

    # correspondence_lengths = [len(v) for k,v in corr_std.iteritems()]

    # most_common_len = max(correspondence_lengths, key = correspondence_lengths.count)

    # missing_correspondence = [ife for ife, units in corr_complete.iteritems() if len(units) != query_length]
  
    # for ife, correspondence in corr_std.iteritems():
    #     if len(correspondence) != most_common_len:
    #         missing_correspondence.append(ife)

	# if missing_correspondence:
	# 	for ife in missing_correspondence:
	# 		del corr_complete[ife]
	# 		del corr_std[ife]



def check_missing_correspondence_new(corr_complete, corr_std, chain):

    correspondence_lengths = [len(v) for k,v in corr_std.iteritems()]

    most_common_len = max(correspondence_lengths, key = correspondence_lengths.count)

    entries_to_remove = []
  
    for ife, correspondence in corr_std.iteritems():
        if len(correspondence) != most_common_len:
            entries_to_remove.append(ife)

    if entries_to_remove:
        for k in entries_to_remove:
            corr_std.pop(k, None)
            corr_complete.pop(k, None)

    return entries_to_remove, corr_complete, corr_std


def order_data(rot, ctr):

      rot_keys = list(rot.keys())
      ctr_keys = list(rot.keys())

      if set(rot_keys) == set(ctr_keys):
            common_keys = rot_keys
      else:
            common_keys = list(set(rot_keys).intersection(ctr_keys))  

      rotation_ordered = [rot.get(k) for k in common_keys]
      center_ordered = [ctr.get(k) for k in common_keys]
      #rotation_ordered_dict = OrderedDict()

    #   missing_data = check_missing_rotation_or_center(rotation_ordered, center_ordered)

      #return rotation_ordered, center_ordered, common_keys, missing_data

    #   if missing_data:
    #       rotation_ordered, center_ordered, common_keys = remove_instances_with_missing_data(common_keys, rotation_ordered, center_ordered, missing_data)
      #for k in common_keys:
          #rotation_ordered_dict[k] = rot.get(k)
    
      return rotation_ordered, center_ordered, common_keys, None


def check_missing_rotation_or_center(rot, ctr):

    missing_rotation_index = set()
    for sublist in rot:
        if None in sublist:
            missing_rotation_index.add(rot.index(sublist))

    missing_center_index = set()
    for sublist in ctr:
        if None in sublist:
            missing_center_index.add(ctr.index(sublist))

    missing_data = missing_rotation_index.union(missing_center_index)

    return list(missing_data)


def remove_instances_with_missing_data(common_keys, rot_ordered, ctr_ordered, missing_data_idxs):

    for elem in reversed(missing_data_idxs):
        del common_keys[elem]
        del rot_ordered[elem]
        del ctr_ordered[elem]

    return rot_ordered, ctr_ordered, common_keys

def calculate_relative_disc(ife_list, center_data, core_len, query_len):
    distances = defaultdict(lambda: defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(a + 1, len(ife_list)):
            try:
                disc = relative_discrepancy(center_data[a], center_data[b],
                                            core_len, query_len)
                distances[ife_list[a]][ife_list[b]] = disc
            except:
                distances[ife_list[a]][ife_list[b]] = None

    return distances 


def calculate_geometric_disc(ife_list, rotation_data, center_data):
    distances = defaultdict(lambda: defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(a + 1, len(ife_list)):
            try:
                disc = matrix_discrepancy(center_data[a], rotation_data[a], center_data[b],
                                        rotation_data[b])
                distances[ife_list[a]][ife_list[b]] = disc
            except:
                distances[ife_list[a]][ife_list[b]] = None

    return distances 

      
      
def order_similarity(ife_list, distances):
    '''
    Calculates the order by similarity between the instances
    Input:
    ife_list: a list of ifes
    distances: a nested dictionary containing the discrepancy values between pairs of ifes

    Output:
    ifes_ordered: a list of tuples where the first element of the tuple would be the order while
                  the second element would be ife
    '''

    dist = np.zeros((len(ife_list), len(ife_list)))
    for index1, member1 in enumerate(ife_list):
        curr = distances.get(member1, {})
        for index2, member2 in enumerate(ife_list):
            dist[index1, index2] = curr.get(member2, 0)

    dist = (dist + np.swapaxes(dist, 0, 1))

    # ordering, _, _ = orderWithPathLengthFromDistanceMatrix(dist, 10, scanForNan=True)
    #disc_order = optimalLeafOrder(dist)

    disc_order = treePenalizedPathLength(dist, 10, 39873)

    ifes_ordered = [(idx, ife_list[order]) for idx, order in enumerate(disc_order)]

    '''
    ordering = []
    order_idx = []

    for idx, order in enumerate(disc_order):
        ordering.append(ife_list[order])
        order_idx.append(idx)

    ifes_ordered = zip(order_idx, ordering)
    ''' 

    return ifes_ordered

def build_heatmap_data_new(distances, ifes_ordered):
    ife1 = []
    ife2 = []
    row_labels = []

    for member1 in ifes_ordered:
        row_labels.append(member1[1])
        for member2 in ifes_ordered:
            ife1.append(member1[1])
            ife2.append(member2[1]) 

    ife_pairs = zip(ife1, ife2)

    disc_ordered = [get_disc(distances, first, second) or get_disc(distances, second, first) for first, second in ife_pairs]

    disc_formatted = []
    for disc in disc_ordered:
        disc = '%.4f' % disc
        disc_formatted.append(disc)

    a = np.array(disc_formatted)
    a = a.astype(np.float)
    max_disc = max(a)

    # How many elements each 
    # list should have 
    sublist_len = len(ifes_ordered)
    
    # using list comprehension 
    disc_list = [disc_formatted[i:i + sublist_len] for i in range(0, len(disc_formatted), sublist_len)]

    return disc_ordered, disc_list, row_labels

def percentile(data, perc):
    size = len(data)
    return sorted(data)[int(math.ceil((size * perc) / 100)) - 1]

def partition_list(lst, n):
    return [lst[i:i+n] for i in range(0, len(lst), n)]

def build_heatmap_data(distances, ifes_ordered):
    index1 = []
    index2 = []
    ife1 = []
    ife2 = []

    for member1 in ifes_ordered:
        for member2 in ifes_ordered:
            index1.append(member1[0])
            ife1.append(member1[1])
            index2.append(member2[0])
            ife2.append(member2[1])

    ife_pairs = zip(ife1, ife2)

    disc_ordered = [get_disc(distances, first, second) or get_disc(distances, second, first) for first, second in ife_pairs]

    disc_formatted = []
    for disc in disc_ordered:
        # This will lead to error if not a number
        if disc:
            disc = '%.4f' % disc
            disc_formatted.append(disc)
        else:
            disc_formatted.append(disc)

    disc_lists = partition_list(disc_formatted, len(ifes_ordered))
    ife_list = [str(x[1]) for x in ifes_ordered]

    cleaned_disc = [float(x) for x in disc_formatted if x is not None]
    cleaned_disc.sort()

    # a = np.array(disc_formatted)
    # a = a.astype(np.float)

    # disc_filtered = [float(x) for x in disc_formatted if x != 'nan']
    percentile_score = percentile(cleaned_disc, 95)
    maximum_disc = '{0:.2f}'.format(max(cleaned_disc))
    # percentile_score = '%.4f' % percentile_score

    #return sorted_disc, None

    # .percentile/mean/mediab/amax doesn't appear to be working
    #percentile = np.percentile(map(float,a), 95)
    #mean = "{:.3f}".format(np.mean(a))
    #median = "{:.3f}".format(np.median(a))
    #max_disc = "{:.3f}".format(np.amax(a))
    # max_disc = max(a)

    # heatmap_data = [
    #     {"ife1": if1, "ife1_index": if1_index, "ife2": if2, "ife2_index": if2_index, "discrepancy": discrepancy}
    #     for if1, if1_index, if2, if2_index, discrepancy in zip(ife1, index1, ife2, index2, disc_formatted)
    # ]

    heatmap_data = ["#chart", disc_lists, ife_list]

    #disc_pairwise = zip(ife1, ife2, disc_filtered)

    return heatmap_data, percentile_score, maximum_disc

def build_heatmap_data_revised(distances, ifes_ordered):
    ife_indices, ife_values = zip(*ifes_ordered)
    num_ifes = len(ifes_ordered)

    disc_ordered = []
    for i in range(num_ifes):
        for j in range(num_ifes):
            ife_pair = (ife_values[i], ife_values[j])
            disc = get_disc(distances, *ife_pair) or get_disc(distances, *ife_pair[::-1])
            disc_ordered.append(disc)

    disc_formatted = ['%.4f' % d if d is not None else d for d in disc_ordered]

    disc_lists = partition_list(disc_formatted, num_ifes)

    cleaned_disc = [d for d in disc_ordered if d is not None]
    cleaned_disc.sort()

    percentile_score = percentile(cleaned_disc, 95)
    maximum_disc = '%.2f' % max(cleaned_disc)

    heatmap_data = ["#chart", disc_lists, [str(x) for x in ife_values]]

    return heatmap_data, percentile_score, maximum_disc


def build_coord_data(ifes_ordered, corr_data):
  
  coord_ordered = OrderedDict()

  # what if the correspondence is not available? Is it possible?
  table_residue_rows = []

  for ife in ifes_ordered:
    correspondence = corr_data.get(ife[1])
    table_residue_rows.append(correspondence)
    correspondence = ','.join(correspondence)
    coord_ordered[ife[1]] = correspondence
   
  return coord_ordered, table_residue_rows


def build_coord_data_unordered(corr_data):

    coord_unordered = OrderedDict()

    for k, v in corr_data.iteritems():
        correspondence = ','.join(v)
        coord_unordered[k] = correspondence

    return coord_unordered


def get_sorted_units(units):
    unsorted_units = units.split(',')
    # This assumes the last element after the split operation to be an integer
    sorted_units = sorted(unsorted_units, key=lambda x: int(x.split('|')[4]))
    return sorted_units

# here we are assuming there are only 2 methods
def get_exp_method_name(exp_method):

    if exp_method.lower() == 'xray':
        return 'X-RAY DIFFRACTION'
    elif exp_method.lower() == 'em': 
        return 'ELECTRON MICROSCOPY'
    else:
        return 'all'


def format_pairwise_interactions_display(corr_data, nt_pairs_positions, pairwise_data, query_len):

    nt_pairs = list(nt_pairs_positions)
    nt_pairs.sort()

    ife_list = pairwise_data.keys()

    interactions_dict = {}
    for chain in ife_list:
        pairwise_interactions_ordered = OrderedDict()
        for unique_res_pair in nt_pairs:
            pairwise_interactions_ordered[unique_res_pair] = ''

        interactions_dict[chain] = pairwise_interactions_ordered

    for chain in ife_list:
        for res_pair in nt_pairs:
            interactions_dict[chain][res_pair] = pairwise_data.get(chain, {}).get(res_pair)

    header_index = [str(x) for x in range(1, query_len+1)]
    header_index = '    '.join(header_index)
    header_res_pairs = '    '.join(nt_pairs)
    header = header_index + '   ' + header_res_pairs

    display_str = ""
    display_str += header + " </br>"

    for ife, interaction_data in interactions_dict.iteritems():
        part_1 = '  '.join(corr_data[ife])
        part_2 = '  '.join(str(interaction) for interaction in interaction_data.values())
        display_str += part_1 + "   " + part_2 + " </br>"
     

    return display_str, nt_pairs


def format_pairwise_interactions_single_display(pairwise_interactions, nt_pairs, query_units):

    residue_pairs = list(nt_pairs)
    residue_pairs = sorted(residue_pairs, key=lambda x: (x[0], x[1]))
    #residue_pairs.sort()

    pairwise_interactions_sorted = [pairwise_interactions.get(residue_pair, None) for residue_pair in residue_pairs]

    residue_pairs_annotated = ['nt'+str(x[0])+'-nt'+str(x[1]) for x in residue_pairs]

    query_len = len(query_units)
    header_index = [str(x) for x in range(1, query_len+1)]
    header_index = '    '.join(header_index)
    header_res_pairs = '    '.join(residue_pairs_annotated)
    header = header_index + '   ' + header_res_pairs

    row_part_1  =  '    '.join(query_units)
    row_part_2  =  '    '.join(pairwise_interactions_sorted)
    row = row_part_1  +  '    '  + row_part_2

    display_str = ""
    display_str += header + " </br>"
    display_str += row 

    return str(display_str)


def format_pairwise_interactions(data):
    pairwise_data, nt_pairs_positions = data
    nt_pairs = sorted(list(nt_pairs_positions))
    ifes_list = pairwise_data.keys()

    interactions_dict = {}
    for ife in ifes_list:
        pairwise_interactions_ordered = OrderedDict.fromkeys(nt_pairs, "")
        interactions_dict[ife] = pairwise_interactions_ordered

    for ife in ifes_list:
        for nt_pair in nt_pairs:
            interactions_dict[ife][nt_pair] = pairwise_data.get(ife, {}).get(nt_pair, "")

    return (nt_pairs, interactions_dict)            


def get_chain_id(query_units):
    chain_id = "|".join(query_units[0].split("|")[:3])
    return chain_id


def build_loop_name(corr_complete):

    loop_names = []
    for ife, loop_range in corr_complete.iteritems():
        range_start = loop_range[0].split('|')[-1]
        range_end = loop_range[1].split('|')[-1]
        pdb = loop_range[0].split('|')[0]
        model_num = loop_range[0].split('|')[1]
        chain = loop_range[0].split('|')[2]
        name = model_num + '/' + chain + '/' + range_start + ':' + range_end
        loop_names.append((pdb, name))

    return loop_names


def display_loop_correspondences(loop_correspondences):

    display_str = ''

    for loop in loop_correspondences:
        display_str += loop[0] + '  ' + loop[1] + '</br>'

    return display_str


def get_correspondence_positions(correspondence):
    
    positions_ref = {}
    for ife, corr in correspondence.iteritems():
        positions = {}
        for count, value in enumerate(corr, start=1):
            positions[count] = "|".join(value.split("|")[3:])

        positions_ref[ife] = positions

    return positions_ref


def get_positions_header(query_len):

    return [ int(x) for x in range(1, query_len+1)]


def process_query_units(query_units):

    pdb = query_units[0].split('|')[0]
    model = query_units[0].split('|')[1]
    chain = query_units[0].split('|')[2]
    ife = "|".join(query_units[0].split('|')[:3])
    units_str = ['|'.join(unit.split('|')[3:]) for unit in query_units]
    units_str = ', '.join(units_str)
    units_length = len(query_units)

    return {"pdb": pdb, "model": model, "chain": chain, "units_list": query_units, "units_str": units_str, "units_length": units_length, "ife": ife}

def get_resolution_data_ordered(ifes_ordered, resolution_dict):
    resolution_data = OrderedDict()
    for ife in ifes_ordered:
        pdb = ife[1].split("|")[0]
        resolution_data[str(ife[1])] = resolution_dict.get(pdb, "")
    return resolution_data


'''
def percentile(N, percent, key=lambda x:x):
    """
    Find the percentile of a list of values.

    @parameter N - is a list of values. Note N MUST BE already sorted.
    @parameter percent - a float value from 0.0 to 1.0.
    @parameter key - optional key function to compute value from each element of N.

    @return - the percentile of the values
    """
    if not N:
        return None
    k = (len(N)-1) * percent
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return key(N[int(k)])
    d0 = key(N[int(f)]) * (c-k)
    d1 = key(N[int(c)]) * (k-f)
    return d0+d1
'''