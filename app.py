from flask import Flask
from flask import render_template, request
from flask import make_response              # for motif variability, to return plain text
import process_input as pi
import query_service as qs
import equivalence_class_service as ec
import correspondence_service as cs
import pairwise_service as ps
from rotation import get_rotation
from center import get_center
import logging
import utility as ui
import json
import time
from collections import OrderedDict, defaultdict
from discrepancy import matrix_discrepancy
import numpy as np
import sys
import time
from get_neighboring_chains import test_run


#from flask_cors import CORS   # for circular

app = Flask(__name__, template_folder='templates')
#CORS(app,expose_headers=["x-suggested-filename"])  # for circular


accepted_resolutions = ['1.5', '2', '2.0', '2.5', '3', '3.0', '3.5', '4', '4.0', 'all']
accepted_disc_methods = ['geometric', 'relative']

@app.route('/')
def home():
    return render_template("index.html")

@app.route('/SVS')
def SVS_home():
    # Manisha: in the future, we will have an argument like "input_form=True"
    # When that happens on this route, we will process the arguments that are provided
    # and supply those to the template below.  But we don't know how to do that!
    # http://rna.bgsu.edu/correspondence/SVS?input_form=True&selection=IL_4V9F_007 
    return render_template("variability.html",input_parameters=request.args)


@app.route('/list')
def display_correspondence():

    query_parameters = request.args

    exp_method = query_parameters.get('exp_method')
    resolution = str(query_parameters.get('resolution_threshold'))
    chain_id = query_parameters.get('chain')
    loop_id = query_parameters.get('loop_id')
    unit_id = query_parameters.get('unit_id')
    res_num = query_parameters.get('res_num')

    if resolution not in accepted_resolutions:
        return 'Please enter a valid resolution. The accepted resolution values \
                are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

    input_type = pi.check_input_type(loop_id, unit_id, res_num)

    if input_type == 'res_num' and chain_id is None:

        return "Please enter the chain parameter"

    query_units = qs.get_query_units_new(input_type, loop_id, unit_id, \
                                         res_num, chain_id)

    if input_type == 'unit_id' or input_type == 'loop_id':
        chain_id = ui.get_chain_id(query_units)

    exp_method = ui.get_exp_method_name(exp_method)

    # Get the equivalence class members that the query ife belongs to
    # Note: 2022-01-06 CLZ this has 4 inputs, but the method only takes 3
    # Likely to crash
    members, ec_name, nr_release = ec.get_ec_members(query_units, resolution, exp_method, chain_id)

    # Get the correspondences for all the members in the equivalence class
    correspondence, corr_complete, corr_std = cs.get_correspondence(query_units, members)

    corr_display = ui.format_correspondence_display(corr_complete)

    return corr_display


@app.route('/pairwise_structure')
def pairwise_correspondence():

    query_parameters = request.args

    chain1 = query_parameters.get('chain1')
    chain2 = query_parameters.get('chain2')

    correspondence = cs.pairwise_structure_correspondence(chain1, chain2)

    corr_display = ui.format_pairwise_correspondence_display(correspondence)

    return corr_display

def correspondence_between_species(parameters_dict):

    start = time.time()

    parameters_dict['selection'], query_error_message = qs.get_query_units_modified(parameters_dict)

    if query_error_message:
        return str(query_error_message)

    query_info, equivalence_class_dict, correspondence_list = ui.get_correspondence_across_species(parameters_dict)

    if not query_info:
        return "Your query returned no results. Please try a different query"

    query_info['chain_name'] = ec.get_chain_standardized_name(query_info)

    corr_complete = ui.get_correspondence_dict(correspondence_list)

    pdb_entries = [k.split("|")[0] for k in corr_complete]

    unique_pdb_entries = list(set(pdb_entries))

    exp_method = ui.get_exp_method_name(parameters_dict['exp_method'])

    if (exp_method != "all"):
        filtered_members = ec.get_members_across_species(unique_pdb_entries, exp_method)
        if not filtered_members:
            return "No PDB entries matches the given selection. Please try a different query"
        corr_complete = ui.filter_dict_by_list(corr_complete, filtered_members)

    correspondence = [item for sublist in correspondence_list for item in sublist]

    nt_position_index = ui.get_nt_position_index(corr_complete)

    possible_nt_pairs = ui.create_all_nt_pairs(corr_complete)

    pairwise_data = ps.get_pairwise_interactions_new(possible_nt_pairs, nt_position_index)

    # pairwise_data = ps.get_pairwise_interactions(corr_complete)

    formatted_pairwise_data = ui.format_pairwise_interactions(pairwise_data)

    correspondence_positions = ui.get_correspondence_positions(corr_complete)

    positions_header = ui.get_positions_header(len(query_info['query_nts_list']))

    # Get rotation data
    rotation_data = get_rotation(correspondence, corr_complete)

    # Get center data
    center_data = get_center(correspondence, corr_complete)

    # Order rotation and center data before computing discrepancy
    rotation_ordered, center_ordered, ife_list, missing_data = ui.order_data(rotation_data, center_data)

    # Calculate geometric discrepancy
    disc_data = ui.calculate_geometric_disc(ife_list, rotation_ordered, center_ordered)

    # Get the instances ordered according to similarity
    ifes_ordered = ui.order_similarity(ife_list, disc_data)

    # Get discrepancy statistics and build the heatmap data for display
    heatmap_data, percentile_score, max_discepancy = ui.build_heatmap_data_revised(disc_data, ifes_ordered)

    # Update the query info dict
    query_info['max_discrepancy'] = max_discepancy
    query_info['percentile_score'] = percentile_score

    # Build coord data
    coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

    units_string = [str(v) for k, v in coord_data.iteritems()]

    # Get the neighboring chains
    neighboring_chains = test_run(units_string)

    # Get the ordered chains as a list
    ifes_ordered_keys = list(coord_data.keys())

    # Get all chain-related information for the entries to be displayed
    chain_info = ec.get_chain_info_dict(ifes_ordered_keys, equivalence_class_dict)

    # Zip both the chain and neighboring chains lists into a dict
    neighboring_chains_dict = OrderedDict(zip(ifes_ordered_keys, neighboring_chains))

    neighboring_chains_list = [row[1] for _, v in neighboring_chains_dict.iteritems() if v for row in v]

    neighboring_chains_count = ui.get_name_count(neighboring_chains_list)

    end = time.time()

    time_diff = '{0:.2f}'.format(end-start)

    output_data = {
        "query_info": query_info,
        "data": heatmap_data,
        "coord": coord_data,
        "res_position": correspondence_positions,
        "positions_header": positions_header,
        "pairwise_data": formatted_pairwise_data,
        "chain_info": chain_info,
        "neighboring_chains": neighboring_chains_dict,
        "neighboring_chains_count": neighboring_chains_count,
        "code_time": time_diff
    }

    return output_data

def correspondence_within_species(parameters_dict):

    equivalence_class_dict = {}

    start = time.time()

    complete_query_units, query_error_message = qs.get_query_units_modified(parameters_dict)

    try:
        if query_error_message:
            return query_error_message

        status_text = "Got units<br>"

        if len(complete_query_units) == 0:
            return "Not able to find units in " + str(selection)

        single_chain_query = ui.check_valid_single_chain_query(complete_query_units)

        if not single_chain_query:
            return "Your query contains nucleotide selections from more than a single chain. Please limit your selections to a single chain"

        complete_query_units_str = ",".join(complete_query_units)

        query_info = ui.process_query_units(complete_query_units)

        query_info['chain_name'] = ec.get_chain_standardized_name(query_info)

        status_text += "Got query_info<br>"

        # if input_type == 'unit_id' or input_type == 'loop_id':
        #     chain_id = ui.get_chain_id(complete_query_units)

        # status_text += "Got chain_id<br>"

        exp_method = ui.get_exp_method_name(parameters_dict['exp_method'])

        status_text += "Got exp_method<br>"

        source_organism = ec.get_source_organism(query_info['pdb'], query_info['chain'])

        status_text += "Got source_organism<br>"

        # Get the equivalence class members that the query ife belongs to
        members, ec_name, nr_release, error_msg = ec.get_ec_members(parameters_dict['resolution'], exp_method, query_info['ife'])

        if len(error_msg) > 0:
            return error_msg

        status_text += "Got equivalence class members<br>"

        # Check whether the selection chain has the same exp_method as in the selection
        empty_members, method_equality = ec.check_valid_membership(members, query_info, exp_method)

        status_text += "Checked valid membership<br>"

        if empty_members is True and method_equality is False:
            return "The current selection has no results returned. Please try a different selection"

        # Filter out PDB entries in the exclude parameter if any
        if parameters_dict['exclude'] is not None:
            members = ui.filter_exclude_ids(members, parameters_dict['exclude'])

        # Get the correspondences for all the members in the equivalence class
        correspondence, corr_complete, corr_std = cs.get_correspondence(complete_query_units, members, method_equality)

        status_text += "Got correspondences<br>"

        # Remove ec member/s that have missing correspondence
        entries_with_missing_correspondence, corr_complete = ui.check_missing_correspondence(corr_complete, query_info['units_length'])

        status_text += "Removed missing members<br>"

        nt_position_index = ui.get_nt_position_index(corr_complete)

        possible_nt_pairs = ui.create_all_nt_pairs(corr_complete)

        pairwise_data = ps.get_pairwise_interactions_new(possible_nt_pairs, nt_position_index)

        status_text += "Got pairwise interactions<br>"

        formatted_pairwise_data = ui.format_pairwise_interactions(pairwise_data)

        status_text += "Formatted pairwise interactions<br>"

        correspondence_positions = ui.get_correspondence_positions(corr_complete)

        status_text += "Got correspondence positions<br>"

        positions_header = ui.get_positions_header(query_info['units_length'])

        # Get rotation data
        rotation_data = get_rotation(correspondence, corr_complete)

        # Get center data
        center_data = get_center(correspondence, corr_complete)

        status_text += "Got %s center and %s rotation data<br>" % (len(center_data),len(rotation_data))

        # Order rotation and center data before computing discrepancy
        rotation_ordered, center_ordered, ife_list, missing_data = ui.order_data(rotation_data, center_data)

        status_text += "Ordered center and rotation data<br>"

        disc_data = ui.calculate_geometric_disc(ife_list, rotation_ordered, center_ordered)

    except:
        return status_text + "<br>... and then something went wrong"

    # Get the instances ordered according to similarity
    ifes_ordered = ui.order_similarity(ife_list, disc_data)

    # Get discrepancy statistics and build the heatmap data for display
    heatmap_data, percentile_score, max_disc = ui.build_heatmap_data(disc_data, ifes_ordered)

    # Build coord data
    coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

    units_string = [str(v) for k, v in coord_data.iteritems()]

    # Get the neighboring chains
    neighboring_chains = test_run(units_string)

    # Get the ordered chains as a list
    ifes_ordered_keys = list(coord_data.keys())

    # Get all chain-related information for the entries to be displayed
    chain_info = ec.get_chain_info_dict(ifes_ordered_keys, equivalence_class_dict)

    # Zip both the chain and neighboring chains lists into a dict
    neighboring_chains_dict = OrderedDict(zip(ifes_ordered_keys, neighboring_chains))

    neighboring_chains_list = [row[1] for _, v in neighboring_chains_dict.iteritems() if v for row in v]

    neighboring_chains_count = ui.get_name_count(neighboring_chains_list)

    end = time.time()

    time_diff = '{0:.2f}'.format(end-start)

    output_data = {
        "data": heatmap_data,
        "max_disc": max_disc,
        "coord": coord_data,
        "ec_name": ec_name,
        "nr_release": nr_release,
        "code_time": time_diff,
        "res_position": correspondence_positions,
        "positions_header": positions_header,
        "pairwise_data": formatted_pairwise_data,
        "selection_data": query_info,
        "percentile": percentile_score,
        "organism": source_organism,
        "neighboring_chains": neighboring_chains_dict,
        "chain_info": chain_info,
        "neighboring_chains_count": neighboring_chains_count
    }

    return output_data


@app.route('/comparison')
def geometric_correspondence_new():

    valid_resolutions = ["1.5", "2.0", "2.5", "3.0", "3.5", "4.0"]

    query_parameters = request.args

    selection = query_parameters.get('selection')
    pdb_id = query_parameters.get('pdb', default=None)
    if pdb_id is not None:
        pdb_id = pdb_id.upper()
    chain_id = query_parameters.get('chain', default=None)
    exp_method = query_parameters.get('exp_method', default='all')
    scope = query_parameters.get('scope', default='EC')
    resolution = query_parameters.get('resolution', default='4.0')
    depth = query_parameters.get('depth')
    exclude = query_parameters.get('exclude', default=None)
    input_form = query_parameters.get('input_form', default='false')

    parameters_dict = {'selection': selection, 'pdb': pdb_id, 'chain': chain_id, 'scope': scope, 'resolution': resolution, 'depth': depth, 'exp_method': exp_method, 'exclude': exclude}

    if input_form.lower() == 'true':
        return render_template("index_form.html", input_parameters=parameters_dict)

    if parameters_dict['resolution'] not in valid_resolutions:
        return "Please enter a valid resolution value. Accepted values are " + ", ".join(valid_resolutions)

    if parameters_dict['scope'] == "EC":
        final_output = correspondence_within_species(parameters_dict)
        if isinstance(final_output, str):
            return str(final_output)
        else:
            return render_template("comparison_ec.html", **final_output)
    else:
        final_output = correspondence_between_species(parameters_dict)
        if isinstance(final_output, str):
            return str(final_output)
        else:
            return render_template("comparison_rfam.html", **final_output)


# @app.route('/comparison_across_species')
# def geometric_correspondence_across_species():

#     valid_scope_values = ["Rfam", "EC", "molecule"]
#     valid_resolutions = ["1.5", "2.0", "2.5", "3.0", "3.5", "4.0"]

#     start = time.time()

#     query_parameters = request.args

#     selection = query_parameters.get('selection')
#     exp_method = query_parameters.get('exp_method', default='all')
#     scope = query_parameters.get('scope', default='Rfam')
#     resolution = query_parameters.get('resolution', default='4.0')
#     depth = query_parameters.get('depth', default=1)
#     exclude = query_parameters.get('exclude', default=None)

#     parameters_dict = {'selection': selection, 'scope': scope, 'resolution': resolution, 'depth': depth}

#     if parameters_dict['scope'] not in valid_scope_values:
#         return "Please enter a valid scope value. Accepted values are " + ", ".join(valid_scope_values)

#     if parameters_dict['resolution'] not in valid_resolutions:
#         return "Please enter a valid resolution value. Accepted values are " + ", ".join(valid_resolutions)

#     query_info, equivalence_class_dict, correspondence_list = ui.get_correspondence_across_species(parameters_dict)

#     if not query_info:
#         return "Your query returned no results. Please try a different query"

#     query_info['chain_name'] = ec.get_chain_standardized_name(query_info)

#     corr_complete = ui.get_correspondence_dict(correspondence_list)

#     pdb_entries = [k.split("|")[0] for k in corr_complete]

#     unique_pdb_entries = list(set(pdb_entries))

#     exp_method = ui.get_exp_method_name(exp_method)

#     if (exp_method != "all"):
#         filtered_members = ec.get_members_across_species(unique_pdb_entries, exp_method)
#         if not filtered_members:
#             return "No PDB entries matches the given selection. Please try a different query"
#         corr_complete = ui.filter_dict_by_list(corr_complete, filtered_members)

#     if exclude:
#         corr_complete = ec.exclude_pdb(corr_complete, exclude)

#     correspondence = [item for sublist in correspondence_list for item in sublist]

#     nt_position_index = ui.get_nt_position_index(corr_complete)

#     possible_nt_pairs = ui.create_all_nt_pairs(corr_complete)

#     pairwise_data = ps.get_pairwise_interactions_new(possible_nt_pairs, nt_position_index)

#     # pairwise_data = ps.get_pairwise_interactions(corr_complete)

#     formatted_pairwise_data = ui.format_pairwise_interactions(pairwise_data)

#     correspondence_positions = ui.get_correspondence_positions(corr_complete)

#     positions_header = ui.get_positions_header(len(query_info['query_nts_list']))

#     # Get rotation data
#     rotation_data = get_rotation(correspondence, corr_complete)

#     # Get center data
#     center_data = get_center(correspondence, corr_complete)

#     # Order rotation and center data before computing discrepancy
#     rotation_ordered, center_ordered, ife_list, missing_data = ui.order_data(rotation_data, center_data)

#     # Calculate geometric discrepancy
#     disc_data = ui.calculate_geometric_disc(ife_list, rotation_ordered, center_ordered)

#     # Get the instances ordered according to similarity
#     ifes_ordered = ui.order_similarity(ife_list, disc_data)

#     # Get discrepancy statistics and build the heatmap data for display
#     heatmap_data, percentile_score, max_discepancy = ui.build_heatmap_data_revised(disc_data, ifes_ordered)

#     # Update the query info dict
#     query_info['max_discrepancy'] = max_discepancy
#     query_info['percentile_score'] = percentile_score

#     # Build coord data
#     coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

#     units_string = [str(v) for k, v in coord_data.iteritems()]

#     # Get the neighboring chains
#     neighboring_chains = test_run(units_string)

#     # Get the ordered chains as a list
#     ifes_ordered_keys = list(coord_data.keys())

#     # Get all chain-related information for the entries to be displayed
#     chain_info = ec.get_chain_info_dict(ifes_ordered_keys, equivalence_class_dict)

#     # Zip both the chain and neighboring chains lists into a dict
#     neighboring_chains_dict = OrderedDict(zip(ifes_ordered_keys, neighboring_chains))

#     neighboring_chains_list = [row[1] for _, v in neighboring_chains_dict.iteritems() if v for row in v]

#     neighboring_chains_count = ui.get_name_count(neighboring_chains_list)

#     # species_name_list = [chain_info[k]['source'] for k, _ in chain_info.iteritems()]

#     # species_name_count = ui.get_name_count(species_name_list)

#     end = time.time()

#     time_diff = '{0:.2f}'.format(end-start)

#     return render_template("comparison_rfam.html", query_info=query_info, data=heatmap_data,
#                             coord=coord_data, code_time=time_diff, res_position=correspondence_positions,
#                             positions_header=positions_header, pairwise_data=formatted_pairwise_data,
#                             chain_info=chain_info, neighboring_chains=neighboring_chains_dict,
#                             neighboring_chains_count=neighboring_chains_count)

# @app.route('/comparison')
# def geometric_correspondence():

#     start = time.time()

#     query_parameters = request.args

#     disc_method = query_parameters.get('disc_method', default='geometric')
#     exp_method = query_parameters.get('exp_method')
#     resolution = str(query_parameters.get('resolution'))
#     chain_id = query_parameters.get('chain')
#     loop_id = query_parameters.get('loop_id')
#     unit_id = query_parameters.get('unit_id')
#     res_num = query_parameters.get('res_num')
#     input_type = query_parameters.get('selection_type')
#     selection = query_parameters.get('selection')
#     # disc_method = query_parameters.get('disc_method')
#     core_nts = query_parameters.get('core_res')

#     equivalence_class_dict = {}

#     if resolution not in accepted_resolutions:
#         return 'Please enter a valid resolution threshold. The accepted resolution values \
#                 are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

#     # if disc_method not in accepted_disc_methods:
#     #     return 'Please correct disc_method to be one of %s' % ", ".join(accepted_disc_methods)

#     #input_type = pi.check_input_type(loop_id, unit_id, res_num)

#     if input_type == 'res_num' and chain_id is None:
#         return "Please enter the chain parameter"

#     # return input_type + " " + selection + " " + chain_id

#     # return str(selection)

#     try:
#         if disc_method == 'geometric':
#             complete_query_units, chain_id = qs.get_query_units_modified(selection)
#         elif disc_method == 'relative':
#             query_units = qs.get_query_units_new(input_type, selection, chain_id)
#             core_units = qs.get_query_units_new(input_type, core_nts, chain_id)
#             complete_query_units = core_units + query_units
#     except:
#         return "Not able to find units in %s" % selection

#     try:
#         status_text = "Got units<br>"

#         if len(complete_query_units) == 0:
#             return "Not able to find units in " + str(chain_id) + " " + str(selection)

#         single_chain_query = ui.check_valid_single_chain_query(complete_query_units)

#         if not single_chain_query:
#             return "Your query contains nucleotide selections from more than a single chain. Please limit your selections to a single chain"

#         complete_query_units_str = ",".join(complete_query_units)

#         query_info = ui.process_query_units(complete_query_units)

#         query_info['chain_name'] = ec.get_chain_standardized_name(query_info)

#         status_text += "Got query_info<br>"

#         if input_type == 'unit_id' or input_type == 'loop_id':
#             chain_id = ui.get_chain_id(complete_query_units)

#         status_text += "Got chain_id<br>"

#         exp_method = ui.get_exp_method_name(exp_method)

#         status_text += "Got exp_method<br>"

#         source_organism = ec.get_source_organism(chain_id)

#         status_text += "Got source_organism<br>"

#         # Get the equivalence class members that the query ife belongs to
#         members, ec_name, nr_release, error_msg = ec.get_ec_members(resolution, exp_method, chain_id)

#         if len(error_msg) > 0:
#             return error_msg

#         status_text += "Got equivalence class members<br>"

#         # Check whether the selection chain has the same exp_method as in the selection
#         empty_members, method_equality = ec.check_valid_membership(members, query_info, exp_method)

#         status_text += "Checked valid membership<br>"

#         if empty_members is True and method_equality is False:
#             return "The current selection has no results returned. Please try a different selection"

#         # Get the correspondences for all the members in the equivalence class
#         correspondence, corr_complete, corr_std = cs.get_correspondence(complete_query_units, members, method_equality)

#         #correspondence_data_length = set()
#         #for k, v in corr_std.iteritems():
#             #correspondence_data_length.add(len(v))
#         status_text += "Got correspondences<br>"

#         # Remove ec member/s that have missing correspondence
#         entries_with_missing_correspondence, corr_complete = ui.check_missing_correspondence(corr_complete, query_info['units_length'])

#         status_text += "Removed missing members<br>"

#         nt_position_index = ui.get_nt_position_index(corr_complete)

#         possible_nt_pairs = ui.create_all_nt_pairs(corr_complete)

#         pairwise_data = ps.get_pairwise_interactions_new(possible_nt_pairs, nt_position_index)

#         status_text += "Got pairwise interactions<br>"

#         formatted_pairwise_data = ui.format_pairwise_interactions(pairwise_data)

#         status_text += "Formatted pairwise interactions<br>"

#         correspondence_positions = ui.get_correspondence_positions(corr_complete)

#         status_text += "Got correspondence positions<br>"

#         positions_header = ui.get_positions_header(query_info['units_length'])

#         # Get rotation data
#         rotation_data = get_rotation(correspondence, corr_complete)

#         # Get center data
#         center_data = get_center(correspondence, corr_complete)

#         status_text += "Got %s center and %s rotation data<br>" % (len(center_data),len(rotation_data))

#         # Order rotation and center data before computing discrepancy
#         rotation_ordered, center_ordered, ife_list, missing_data = ui.order_data(rotation_data, center_data)

#         status_text += "Ordered center and rotation data<br>"

#         if disc_method == 'geometric':
#             # Calculate geometric discrepancy
#             disc_data = ui.calculate_geometric_disc(ife_list, rotation_ordered, center_ordered)
#         elif disc_method == 'relative':
#             disc_data = ui.calculate_relative_disc(ife_list, center_ordered, len(core_units), len(query_units))
#     except:
#         return status_text + "<br>... and then something went wrong"

#     # Get the instances ordered according to similarity
#     ifes_ordered = ui.order_similarity(ife_list, disc_data)

#     # Get discrepancy statistics and build the heatmap data for display
#     heatmap_data, percentile_score, max_disc = ui.build_heatmap_data(disc_data, ifes_ordered)

#     # Build coord data
#     coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

#     units_string = [str(v) for k, v in coord_data.iteritems()]

#     # Get the neighboring chains
#     neighboring_chains = test_run(units_string)

#     # Get the ordered chains as a list
#     ifes_ordered_keys = list(coord_data.keys())

#     # Get all chain-related information for the entries to be displayed
#     chain_info = ec.get_chain_info_dict(ifes_ordered_keys, equivalence_class_dict)

#     # Zip both the chain and neighboring chains lists into a dict
#     neighboring_chains_dict = OrderedDict(zip(ifes_ordered_keys, neighboring_chains))

#     neighboring_chains_list = [row[1] for _, v in neighboring_chains_dict.iteritems() if v for row in v]

#     neighboring_chains_count = ui.get_name_count(neighboring_chains_list)

#     end = time.time()

#     time_diff = '{0:.2f}'.format(end-start)

#     return render_template("comparison_ec.html", data=heatmap_data, max_disc=max_disc, coord=coord_data, ec_name=ec_name,
#                             nr_release=nr_release, code_time=time_diff, res_position=correspondence_positions,
#                             positions_header=positions_header, pairwise_data=formatted_pairwise_data,
#                             selection_data=query_info, percentile=percentile_score,
#                             organism=source_organism, neighboring_chains=neighboring_chains_dict, chain_info=chain_info,
#                             neighboring_chains_count=neighboring_chains_count)


@app.route('/pairwise_interactions')
def pairwise_interactions_correspondence():

    start = time.time()

    query_parameters = request.args

    disc_method = query_parameters.get('disc_method')
    exp_method = query_parameters.get('exp_method')
    resolution = str(query_parameters.get('resolution_threshold'))
    chain_id = query_parameters.get('chain')
    loop_id = query_parameters.get('loop_id')
    unit_id = query_parameters.get('unit_id')
    res_num = query_parameters.get('res_num')
    input_type = query_parameters.get('selection_type')
    selection = query_parameters.get('selection')

    if resolution not in accepted_resolutions:
        return 'Please enter a valid resolution. The accepted resolution values \
                are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

    #input_type = pi.check_input_type(loop_id, unit_id, res_num)

    if input_type == 'res_num' and chain_id is None:
        return "Please enter the chain parameter"

    #query_units = qs.get_query_units_new(input_type, loop_id, unit_id, res_num, chain_id)

    query_units = qs.get_query_units_new(input_type, selection, chain_id)

    query_data = ui.process_query_units(query_units)

    if input_type == 'unit_id' or input_type == 'loop_id': chain_id = ui.get_chain_id(query_units)

    exp_method = ui.get_exp_method_name(exp_method)

    # Get the equivalence class members that the query ife belongs to
    members, ec_name, nr_release = ec.get_ec_members(resolution, exp_method, chain_id)

    # Check whether the selection chain has the same exp_method as in the selection
    empty_members, method_equality = ec.check_valid_membership(members, query_data, exp_method)

    if empty_members is True and method_equality is False: return "The current selection has no results returned. Please try a different selection"

    # Get the correspondences for all the members in the equivalence class
    correspondence, corr_complete, corr_std = cs.get_correspondence(query_units, members, method_equality)

    # Remove ec member/s that have missing correspondence
    missing_data, corr_complete, corr_std = ui.check_missing_correspondence(corr_complete, corr_std)

    pairwise_data, pairwise_residue_pairs_reference = ps.get_pairwise_interactions(corr_complete)

    pairwise_interactions_display, res_pairs = ui.format_pairwise_interactions_display(corr_complete, pairwise_residue_pairs_reference, pairwise_data, len(query_units))

    end = time.time()

    time_diff = '{0:.2f}'.format(end-start)

    return str(pairwise_interactions_display)


@app.route('/pairwise_interactions_single')
def pairwise_interactions_single():

    query_parameters = request.args

    chain_id = query_parameters.get('chain')
    input_type = query_parameters.get('selection_type')
    selection = query_parameters.get('selection')

    if input_type == 'res_num' and chain_id is None:
        return "Please enter the chain parameter"

    query_units = qs.get_query_units_new(input_type, selection, chain_id)

    if input_type == 'unit_id' or input_type == 'loop_id': chain_id = ui.get_chain_id(query_units)

    pairwise_interactions, pairwise_residue_pairs = ps.get_pairwise_interactions_single(query_units)

    pairwise_interactions_display = ui.format_pairwise_interactions_single_display(pairwise_interactions, pairwise_residue_pairs, query_units)

    return str(pairwise_interactions_display)


@app.route('/loop')
def loop_correspondence():

    start = time.time()

    query_parameters = request.args

    exp_method = query_parameters.get('exp_method')
    resolution = str(query_parameters.get('resolution_threshold'))
    loop_id = query_parameters.get('loop_id')
    unit_id = query_parameters.get('unit_id')
    res_num = query_parameters.get('res_num')
    chain_id = query_parameters.get('chain')

    if resolution not in accepted_resolutions:
        return 'Please enter a valid resolution. The accepted resolution values \
                are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

    input_type = pi.check_input_type(loop_id, unit_id, res_num)

    if input_type == 'res_num' and chain_id is None:
        return "Please enter the chain parameter"

    query_units = qs.get_query_units_new(input_type, loop_id, unit_id, res_num, chain_id)

    loop_range = [query_units[0], query_units[-1]]

    if input_type == 'unit_id' or input_type == 'loop_id': chain_id = ui.get_chain_id(query_units)

    exp_method = ui.get_exp_method_name(exp_method)

    # Get the equivalence class members that the query ife belongs to
    members, ec_name, nr_release = ec.get_ec_members(query_units, resolution, exp_method, chain_id)

    # Get the correspondences for all the members in the equivalence class
    correspondence, corr_complete, corr_std = cs.get_correspondence(loop_range, members)

    loop_names = ui.build_loop_name(corr_complete)

    loop_correspondences = cs.get_loop_correspondence(loop_names)

    display_str = ui.display_loop_correspondences(loop_correspondences)

    return display_str

@app.route('/variability')
def variability():

    try:
        from motif_variability import get_sequence_variability
    except Exception as inst:
        output = "%s" % type(inst)
        output += "%s" % inst.args
        output += "%s" % inst
        return output

    query_parameters = request.args

    id = query_parameters.get('id')            # better to just allow some kind of id
    #loop_id = query_parameters.get('loop_id')
    #unit_id = query_parameters.get('unit_id')       # if not given, value is False
    extension = query_parameters.get('extension')
    output_format = query_parameters.get('format')
    domain = query_parameters.get('domain')
    codon = query_parameters.get('codon')
    count = query_parameters.get('count')

    #id = id.replace(" ","")                 # in case this helps
    #loop_id = loop_id.replace(" ","")
    #unit_id = unit_id.replace(" ","")

    to_do = []

    if domain:
        domain = domain.split(",")
    else:
        domain = []

    if codon:
        codon = codon.split(",")
    else:
        codon = []

    if extension:
        extension = int(extension)
    else:
        extension = 0

    if not output_format:
        output_format = "full"

    if not count:
        count = ""

    # error checking, prevent code injection
    # if not "top" in count and not "multiplicity" in count and not "_" in count:
    #     count = ""

    output = "Basic usage\nformat=[full,unique,fasta]\ncount=top_90_percent,top_80_percent,top_10_sequences,top_10_sequences,multiplicity_5_or_more,multiplicity_9_or_more\nextension=1 or 2\ndomain=[A,B,C,E,M,V] comma separated\ncodon=comma separated list\n\n"
    problem = False

    # if unit_id:
    #     output += 'New request made for unit_id %s with extension %s and output format %s\n' % (unit_id,extension,output_format)
    #     if not "|" in unit_id:
    #         output += 'Invalid unit id %s\n' % unit_id
    #         problem = True
    #     else:
    #         to_do = [unit_id.split(",")]
    # elif loop_id:
    #     output += 'New request made for loop_id %s with extension %s and output format %s\n' % (loop_id,extension,output_format)
    #     if not len(loop_id) == 11 or not len(loop_id.split("_")) == 3:
    #         output += 'Invalid loop id %s\n' % loop_id
    #         problem = True
    #     else:
    #         to_do = [[loop_id]]
    if id:
        output += 'New request made for id %s with extension %s and output format %s\n' % (id,extension,output_format)
        if not "|" in id and not "_" in id:
            output += 'Invalid id %s\n' % id
            problem = True
        else:
            to_do = [[id]]

    #count = "top_10_sequences"
    source = "Rfam"

    if not problem:
        if output_format in ["full","unique","fasta","top_motif_models"]:

            if output_format == "top_motif_models" and not loop_id:
                output += "top_motif_models only available for loop input\n"

            try:
                output_list = get_sequence_variability(to_do,extension,output_format,count,domain,codon,source)

                output += "output_list length %d\n" % len(output_list)
                output += "output_list[0] length %d\n" % len(output_list[0])

                output = output_list[0][0]
                #family = output_list[0][1]

                if len(output) > 30000000:
                    output = output[:30000000]
                    output += "\nOutput truncated to 30,000,000 characters.  Restrict the domain or count or codon to get more meaningful output."

            except Exception as inst:
                output += "%s\n" % type(inst)
                output += "%s\n" % inst.args
                output += "Something went wrong with this request\n"
                output += "%s\n" % inst

        else:

            output += "Unknown output format"

    response = make_response(output, 200)
    response.mimetype = "text/plain"

    return response


@app.route('/map_across_species')
def map_across_species():

    # http://rna.bgsu.edu/correspondence/map_across_species?id=HL_7K00_033
    # ln /usr/local/alignment/rpfam/alignments/rfam rfam

    try:
        from map_across_species import map_across_species as m_a_s
    except Exception as inst:
        output = "%s" % type(inst)
        output += "%s" % inst.args
        output += "%s" % inst
        return output

    query_parameters = request.args

    id = query_parameters.get('id').replace("'","")            # better to just allow some kind of id

    scope = query_parameters.get('scope','Rfam')
    resolution = query_parameters.get('resolution','3.0A')
    depth = int(query_parameters.get('depth','5'))
    format = query_parameters.get('format','text')
    match = query_parameters.get('match','full')

    output = "Basic usage; defaults indicated in parentheses\nid=loop id or comma-separated list of unit ids\n&scope=[EC,(Rfam),molecule]\n&resolution=[1.5A,2.0A,2.5A,(3.0A),3.5A,4.0A,20.0A,all]\n&depth=integer (5)\n&format=[(text),json]\n&match=[(full),partial]\n\n"
    problem = False

    to_do = []

    if id:
        output += 'New request made for id %s with scope %s resolution %s depth %d format %s match %s\n' % (id,scope,resolution,depth,format,match)
        if not "|" in id and not "_" in id:
            output += 'Invalid id %s\n' % id
            problem = True
        else:
            to_do = [[id]]


    if not problem:
        if scope in ["EC","Rfam","molecule"]:
            if resolution in ['1.5A','2.0A','2.5A','3.0A','3.5A','4.0A','20.0A','all']:

                try:
                    results_list = m_a_s(to_do,scope,resolution,depth,match)
                    if format == 'text':
                        output = results_list[0]["text"]
                    else:
                        output = json.dumps(results_list[0])

                except Exception as e:
                    exception_type, exception_object, exception_traceback = sys.exc_info()
                    line_number = exception_traceback.tb_lineno
                    output += "\nSomething went wrong with this request on line %s with error type %s\n" % (line_number,exception_type)
                    output += "%s\n" % type(e)
                    output += "%s\n" % exception_traceback
                    #output += "%s\n" % inst.args
                    output += "%s\n" % e
                    problem = True

            else:
                output += "Unknown resolution %s" % resolution
                problem = True

        else:

            output += "Unknown scope %s" % scope
            problem = True

    if problem:
        response = make_response(output, 400)
    else:
        response = make_response(output, 200)

    response.mimetype = "text/plain"

    return response


@app.route('/align_chains')
def align_chains():

    # align two or more chains
    # use Rfam / Infernal alignment, until other options become available

    # http://rna.bgsu.edu/correspondence/align_chains?chains=4V9F|1|0,1S72|1|0,1K9M|1|A
    # http://rna.bgsu.edu/correspondence/align_chains?chains=5J7L|1|AA,1J5E|1|A
    # http://rna.bgsu.edu/correspondence/align_chains?chains=5J7L|1|DA,8B0X|1|a,7K00|1|a  # won't work yet!
    # http://rna.bgsu.edu/correspondence/align_chains?chains=5J7L|1|AA,1J5E|1|A,6TH6|1|Aa,4V88|1|A6,6ZMI|1|S2

    output = "Making progress"

    try:
        from align_chains import align_chains as a_c
    except Exception as inst:
        output = "%s" % inst
        return output

    query_parameters = request.args

    chains = query_parameters.get('chains').replace("'","")

    output = "Specify two or more chains separated by commas, for example, 5J7L|1|AA,4Y4O|1|1a\n"
    problem = False

    if chains:
        output += 'New request made for chains %s\n' % (chains)
        if not "|" in chains or not "," in chains:
            output += 'Invalid chains %s\n' % chains
            problem = True

    if not problem:
        try:
            output = a_c(chains.split(","))

        except Exception as e:
            exception_type, exception_object, exception_traceback = sys.exc_info()
            line_number = exception_traceback.tb_lineno
            output += "\nSomething went wrong with this request on line %s with error type %s\n" % (line_number,exception_type)
            output += "%s\n" % type(e)
            output += "%s\n" % exception_traceback
            #output += "%s\n" % inst.args
            output += "%s\n" % e
            problem = True

    if problem:
        response = make_response(output, 400)
    else:
        response = make_response(output, 200)

    response.mimetype = "text/plain"

    return response



@app.route('/circular')
def circular_diagram():

    from flask import send_file
    import os
    import circular_diagram_14

    query_parameters = request.args

    chains_string = str(query_parameters.get('chains',default='',type=str))

    try:
        filename = ""
        filename = circular_diagram_14.main([None,chains_string])
        iii = 3
    except Exception as e:
        return str(e)

    # return 'Trying to create the circular diagram for %s' % (chains_string)


    ps_file  = os.path.join('/var/www/correspondence/circular_diagram_pdf',filename+".ps")
    pdf_file = os.path.join('/var/www/correspondence/circular_diagram_pdf',filename+".pdf")

    os.remove(ps_file)

    try:
        #return send_file(pdf_file, attachment_filename=filename+".pdf")  # wrong download name
        return send_file(pdf_file, attachment_filename=filename+".pdf", as_attachment=True)  # instant download

        #response = send_file(pdf_file, attachment_filename = filename+".pdf")
        #response.headers["x-suggested-filename"] = filename+".pdf"
        #response.headers["Access-Control-Expose-Headers"] = 'x-suggested-filename'
        #return response
    except Exception as e:
        return str(e)


if __name__ == '__main__':
    app.debug = False
    app.run()