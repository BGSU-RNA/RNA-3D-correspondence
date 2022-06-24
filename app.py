from flask import Flask
from flask import render_template, request
import process_input as pi
import query_service as qs
import equivalence_class_service as ec
import correspondence_service as cs
import pairwise_service as ps
from rotation import get_rotation
from center import get_center
import utility as ui
import json
import time
from collections import OrderedDict, defaultdict
from discrepancy import matrix_discrepancy
import numpy as np
import sys
import time

app = Flask(__name__, template_folder='templates')

accepted_resolutions = ['1.5', '2', '2.0', '2.5', '3', '3.0', '3.5', '4', '4.0', 'all']
accepted_disc_methods = ['geometric', 'relative']

@app.route('/')
# Need to work on the homepage
def home():
    
    return render_template("index_new.html")


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


@app.route('/comparison')
def geometric_correspondence():

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
    disc_method = query_parameters.get('disc_method')
    core_nts = query_parameters.get('core_res')

    if resolution not in accepted_resolutions:
        return 'Please enter a valid resolution threshold. The accepted resolution values \
                are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

    if disc_method not in accepted_disc_methods:
        return 'Please correct disc_method to be one of %s' % ", ".join(accepted_disc_methods)

    #input_type = pi.check_input_type(loop_id, unit_id, res_num)

    if input_type == 'res_num' and chain_id is None:
        return "Please enter the chain parameter"

    #return input_type + " " + selection + " " + chain_id

    try:
        if disc_method == 'geometric':
            complete_query_units = qs.get_query_units_new(input_type, selection, chain_id)
        elif disc_method == 'relative':
            query_units = qs.get_query_units_new(input_type, selection, chain_id)
            core_units = qs.get_query_units_new(input_type, core_nts, chain_id)
            complete_query_units = core_units + query_units
    except:
        return "Not able to find units in %s" % selection

    try:
        status_text = "Got units<br>"

        if len(complete_query_units) == 0:
            return "Not able to find units in " + str(chain_id) + " " + str(selection)

        query_data = ui.process_query_units(complete_query_units)

        status_text += "Got query_data<br>"

        if input_type == 'unit_id' or input_type == 'loop_id':
            chain_id = ui.get_chain_id(complete_query_units)

        status_text += "Got chain_id<br>"

        exp_method = ui.get_exp_method_name(exp_method)

        status_text += "Got exp_method<br>"

        source_organism = ec.get_source_organism(chain_id)

        status_text += "Got source_organism<br>"

        # Get the equivalence class members that the query ife belongs to
        members, ec_name, nr_release, error_msg = ec.get_ec_members(resolution, exp_method, chain_id)

        if len(error_msg) > 0:
            return error_msg

        status_text += "Got equivalence class members<br>"

        # Check whether the selection chain has the same exp_method as in the selection
        empty_members, method_equality = ec.check_valid_membership(members, query_data, exp_method)

        status_text += "Checked valid membership<br>"

        if empty_members is True and method_equality is False:
            return "The current selection has no results returned. Please try a different selection"

        # Get the correspondences for all the members in the equivalence class
        correspondence, corr_complete, corr_std = cs.get_correspondence(complete_query_units, members, method_equality)

        #correspondence_data_length = set()
        #for k, v in corr_std.iteritems():
            #correspondence_data_length.add(len(v))
        status_text += "Got correspondences<br>"

        # Remove ec member/s that have missing correspondence
        missing_data, corr_complete, corr_std = ui.check_missing_correspondence(corr_complete, corr_std)

        status_text += "Removed missing members<br>"

        pairwise_data, pairwise_residue_pairs_reference = ps.get_pairwise_interactions(corr_complete)

        status_text += "Got pairwise interactions<br>"

        pairwise_interactions_data, res_pairs = ui.format_pairwise_interactions_table(pairwise_residue_pairs_reference, pairwise_data)

        status_text += "Formatted pairwise interactions<br>"

        correspondence_positions = ui.get_correspondence_positions(corr_complete)

        status_text += "Got correspondence positions<br>"

        query_len = len(complete_query_units)

        positions_header = ui.get_positions_header(query_len)

        # Get rotation data
        rotation_data = get_rotation(correspondence, corr_std)

        # Get center data
        center_data = get_center(correspondence, corr_std)

        status_text += "Got %s center and %s rotation data<br>" % (len(center_data),len(rotation_data))

        # Order rotation and center data before computing discrepancy
        rotation_ordered, center_ordered, ife_list, missing_data = ui.order_data(rotation_data, center_data)

        status_text += "Ordered center and rotation data<br>"

        '''
        rotation_data_length = set()
        for sublist in rotation_ordered:
            rotation_data_length.add(len(sublist))

        center_data_length = set()
        for sublist in center_ordered:
            center_data_length.add(len(sublist))
        '''

        if disc_method == 'geometric':
            # Calculate geometric discrepancy
            disc_data = ui.calculate_geometric_disc(ife_list, rotation_ordered, center_ordered)
        elif disc_method == 'relative':
            disc_data = ui.calculate_relative_disc(ife_list, center_ordered, len(core_units), len(query_units))
    except:
        return status_text + "<br>... and then something went wrong"

    #return "Temporary return string" + "<br>" + str(complete_query_units) + "<br>" + str(disc_data)

    # Get the instances ordered according to similarity
    ifes_ordered = ui.order_similarity(ife_list, disc_data)

    # Get discrepancy statistics and build the heatmap data for display
    max_disc, heatmap_data, percentile_score = ui.build_heatmap_data(disc_data, ifes_ordered)

    # new heatmap method
    #max_disc, heatmap_data, row_labels = ui.build_heatmap_data_new(disc_data, ifes_ordered)

    # Build coord data
    coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

    #table_cols = zip(*table_rows) 

    #return str(correspondence_positions)

    #return str(correspondence_positions['4WOI|1|AA']['Nt1'])

    end = time.time() 

    time_diff = '{0:.2f}'.format(end-start)

    return render_template("comparison_display.html", data=heatmap_data, max_disc=max_disc, coord=coord_data, ec_name=ec_name, 
                            nr_release=nr_release, code_time=time_diff, res_position=correspondence_positions, 
                            positions_header=positions_header, pairwise_interactions=pairwise_interactions_data,
                            interactions_header=res_pairs, selection_data=query_data, percentile=percentile_score,
                            organism=source_organism)


@app.route('/comparison_without_ordering')
def correspondence_without_ordering():

    start = time.time()

    query_parameters = request.args

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
    missing_data, corr_complete, corr_std = ui.check_missing_correspondence_new(corr_complete, corr_std, chain_id)

    pairwise_data, pairwise_residue_pairs_reference = ps.get_pairwise_interactions(corr_complete)

    pairwise_interactions_data, res_pairs = ui.format_pairwise_interactions_table(pairwise_residue_pairs_reference, pairwise_data)

    correspondence_positions = ui.get_correspondence_positions(corr_complete)

    query_len = len(query_units)

    positions_header = ui.get_positions_header(query_len)

    # Build coord data
    coord_data = ui.build_coord_data_unordered(corr_complete)

    #return str(coord_data)

    end = time.time() 

    time_diff = '{0:.2f}'.format(end-start)

    return render_template("comparison_display_without_ordering.html", coord=coord_data, ec_name=ec_name, 
                            nr_release=nr_release, code_time=time_diff, res_position=correspondence_positions, 
                            positions_header=positions_header, pairwise_interactions=pairwise_interactions_data,
                            interactions_header=res_pairs, selection_data=query_data)



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

@app.route('/motif_variability')
def motif_variability():

    #from motif_variability import get_sequence_variability

    query_parameters = request.args

    loop_id = query_parameters.get('loop_id')
    output_format = query_parameters.get('output_format')

    if len(loop_id) > 11:
        output = 'Invalid loop id %s' % loop_id
        return output

    # output = get_sequence_variability(loop_id,output_format)

    output = 'New request made for loop_id %s with output format %s' % (loop_id,output_format)

    return output


if __name__ == '__main__':
    app.debug = False
    app.run()