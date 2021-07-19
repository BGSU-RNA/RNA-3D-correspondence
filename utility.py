from collections import OrderedDict, defaultdict
from discrepancy import matrix_discrepancy
import numpy as np
from ordering import optimalLeafOrder
from ordering_similarity import treePenalizedPathLength

def get_disc(d, first, second):
    return d.get(second, {}).get(first, 0.0)


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


def check_missing_correspondence(corr_complete, corr_std):

    correspondence_lengths = [len(v) for k,v in corr_std.iteritems()]

    most_common_len = max(correspondence_lengths, key = correspondence_lengths.count)

    missing_correspondence = []
  
    for ife, correspondence in corr_std.iteritems():
        if len(correspondence) != most_common_len:
            missing_correspondence.append(ife)

	if missing_correspondence:
		for ife in missing_correspondence:
			del corr_complete[ife]
			del corr_std[ife]

	return missing_correspondence, corr_complete, corr_std


def order_data(rot, ctr):

      rot_keys = list(rot.keys())
      ctr_keys = list(rot.keys())

      if set(rot_keys) == set(ctr_keys):
            common_keys = rot_keys
      else:
            common_keys = list(set(rot_keys).intersection(ctr_keys))  

      rotation_ordered = [rot.get(k) for k in common_keys]
      center_ordered = [ctr.get(k) for k in common_keys]
    
      return rotation_ordered, center_ordered, common_keys

def calculate_geometric_disc(ife_list, rotation_data, center_data):
    distances = defaultdict(lambda: defaultdict(int))

    for a in range(0, len(ife_list)):
        for b in range(a + 1, len(ife_list)):
            disc = matrix_discrepancy(center_data[a], rotation_data[a], center_data[b],
                                      rotation_data[b])
            distances[ife_list[a]][ife_list[b]] = disc

    return distances 
      

def order_similarity(ife_list, distances):
    '''
    Calculates the order the similarity between the instances
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
        disc = '%.4f' % disc
        disc_formatted.append(disc)

    a = np.array(disc_formatted)
    a = a.astype(np.float)

    # .percentile/mean/mediab/amax doesn't appear to be working
    #percentile = np.percentile(map(float,a), 95)
    #mean = "{:.3f}".format(np.mean(a))
    #median = "{:.3f}".format(np.median(a))
    #max_disc = "{:.3f}".format(np.amax(a))
    max_disc = max(a)

    heatmap_data = [
        {"ife1": if1, "ife1_index": if1_index, "ife2": if2, "ife2_index": if2_index, "discrepancy": discrepancy}
        for if1, if1_index, if2, if2_index, discrepancy in zip(ife1, index1, ife2, index2, disc_formatted)
    ]

    #disc_pairwise = zip(ife1, ife2, disc_formatted)

    return max_disc, heatmap_data


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


def format_pairwise_interactions_display(corr_data, res_pairs_ref, pairwise_data, query_len):

    pairwise_residue_pairs = list(res_pairs_ref)
    pairwise_residue_pairs.sort()

    chains_list = pairwise_data.keys()

    pairwise_interactions_collection = {}
    for chain in chains_list:
        pairwise_interactions_ordered = OrderedDict()
        for unique_res_pair in pairwise_residue_pairs:
            pairwise_interactions_ordered[unique_res_pair] = ''

        pairwise_interactions_collection[chain] = pairwise_interactions_ordered

    for chain in chains_list:
        for res_pair in pairwise_residue_pairs:
            pairwise_interactions_collection[chain][res_pair] = pairwise_data.get(chain, {}).get(res_pair)

    header_index = [str(x) for x in range(1, query_len+1)]
    header_index = '    '.join(header_index)
    header_res_pairs = '    '.join(pairwise_residue_pairs)
    header = header_index + '   ' + header_res_pairs

    display_str = ""
    display_str += header + " </br>"

    for ife, interaction_data in pairwise_interactions_collection.iteritems():
        part_1 = '  '.join(corr_data[ife])
        part_2 = '  '.join(str(interaction) for interaction in interaction_data.values())
        display_str += part_1 + "   " + part_2 + " </br>"
     

    return display_str, pairwise_residue_pairs


def format_pairwise_interactions_table(res_pairs_ref, pairwise_data):

    pairwise_residue_pairs = list(res_pairs_ref)
    pairwise_residue_pairs.sort()

    chains_list = pairwise_data.keys()

    pairwise_interactions_collection = {}
    for chain in chains_list:
        pairwise_interactions_ordered = OrderedDict()
        for unique_res_pair in pairwise_residue_pairs:
            pairwise_interactions_ordered[unique_res_pair] = ''

        pairwise_interactions_collection[chain] = pairwise_interactions_ordered

    for chain in chains_list:
        for res_pair in pairwise_residue_pairs:
            pairwise_interactions_collection[chain][res_pair] = pairwise_data.get(chain, {}).get(res_pair)

    return pairwise_interactions_collection, pairwise_residue_pairs            


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
    chain = query_units[0].split('|')[2]
    residues = ['|'.join(unit.split('|')[3:]) for unit in query_units]
    residues = ','.join(residues)

    return (pdb, chain, residues)