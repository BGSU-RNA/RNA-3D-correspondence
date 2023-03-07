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

@app.route('/comparison_across_species')
def geometric_correspondence_across_species():

    valid_scope_values = ["Rfam", "EC", "molecule"]
    valid_resolutions = ["1.5A", "2.0A", "2.5A", "3.0A", "3.5A", "4.0A"]

    start = time.time()

    query_parameters = request.args
    
    loop_id = query_parameters.get('loop_id')
    scope = query_parameters.get('scope')
    resolution = query_parameters.get('resolution')

    parameters_dict = {'loop_id': loop_id, 'scope': scope, 'resolution': resolution}

    if None in parameters_dict.values():
        return "Please specify the values for all the required parameters - id, scope and resolution"

    if parameters_dict['scope'] not in valid_scope_values:
        return "Please enter a valid scope value. Accepted values are " + ", ".join(valid_scope_values)

    if parameters_dict['resolution'] not in valid_resolutions:
        return "Please enter a valid resolution value. Accepted values are " + ", ".join(valid_resolutions)
    
    query_info, equivalence_class_dict, correspondence_list = ui.get_correspondence_across_species(parameters_dict)

    corr_complete = ui.get_correspondence_dict(correspondence_list)

    correspondence = [item for sublist in correspondence_list for item in sublist]

    pairwise_data = ps.get_pairwise_interactions(corr_complete)

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
    heatmap_data, percentile_score, max_discepancy = ui.build_heatmap_data(disc_data, ifes_ordered)

    # Update the query info dict
    query_info['max_discrepancy'] = max_discepancy
    query_info['percentile_score'] = percentile_score

    # Build coord data
    coord_data, table_rows = ui.build_coord_data(ifes_ordered, corr_complete)

    # Get the ordered chains as a list
    ifes_ordered_keys = list(coord_data.keys())

    # Get all chain-related information for the entries to be displayed
    chain_info = ec.get_chain_info_dict(ifes_ordered_keys, equivalence_class_dict)

    species_name_list = [chain_info[k]['source'] for k, _ in chain_info.iteritems()]

    species_name_count = ui.get_name_count(species_name_list)

    end = time.time() 

    time_diff = '{0:.2f}'.format(end-start)

    return render_template("comparison_test_new.html", query_info=query_info, data=heatmap_data, 
                            coord=coord_data, code_time=time_diff, res_position=correspondence_positions, 
                            positions_header=positions_header, pairwise_data=formatted_pairwise_data,
                            chain_info=chain_info, name_count=species_name_count)

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

    equivalence_class_dict = {}

    if resolution not in accepted_resolutions:
        return 'Please enter a valid resolution threshold. The accepted resolution values \
                are 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 and all'

    if disc_method not in accepted_disc_methods:
        return 'Please correct disc_method to be one of %s' % ", ".join(accepted_disc_methods)

    #input_type = pi.check_input_type(loop_id, unit_id, res_num)

    if input_type == 'res_num' and chain_id is None:
        return "Please enter the chain parameter"

    # return input_type + " " + selection + " " + chain_id

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

        complete_query_units_str = ",".join(complete_query_units)

        formatted_query_units = ui.format_query_units(complete_query_units)

        # sequence_count_dict = ui.get_sequence_variability(complete_query_units_str)

        # sequence_logo_data = ui.generate_sequence_logo_data(sequence_count_dict)

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
        missing_data, corr_complete, corr_std = ui.check_missing_correspondence(corr_complete, corr_complete)

        status_text += "Removed missing members<br>"

        pairwise_data = ps.get_pairwise_interactions(corr_complete)

        status_text += "Got pairwise interactions<br>"

        formatted_pairwise_data = ui.format_pairwise_interactions(pairwise_data)

        status_text += "Formatted pairwise interactions<br>"

        correspondence_positions = ui.get_correspondence_positions(corr_complete)

        status_text += "Got correspondence positions<br>"

        query_len = len(complete_query_units)

        positions_header = ui.get_positions_header(query_len)

        # Get rotation data
        rotation_data = get_rotation(correspondence, corr_complete)

        # Get center data
        center_data = get_center(correspondence, corr_complete)

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

    # debugging purpose
    # return str(disc_data)

    # Get the instances ordered according to similarity
    ifes_ordered = ui.order_similarity(ife_list, disc_data)

    # Get discrepancy statistics and build the heatmap data for display
    heatmap_data, percentile_score, max_disc = ui.build_heatmap_data(disc_data, ifes_ordered)

    # new heatmap method
    #max_disc, heatmap_data, row_labels = ui.build_heatmap_data_new(disc_data, ifes_ordered)

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

    return render_template("comparison_test.html", data=heatmap_data, max_disc=max_disc, coord=coord_data, ec_name=ec_name, 
                            nr_release=nr_release, code_time=time_diff, res_position=correspondence_positions, 
                            positions_header=positions_header, pairwise_data=formatted_pairwise_data,
                            selection_data=query_data, percentile=percentile_score,
                            organism=source_organism, neighboring_chains=neighboring_chains_dict, chain_info=chain_info,
                            neighboring_chains_count=neighboring_chains_count, query_units=formatted_query_units)


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
    loop_id = query_parameters.get('loop_id')
    unit_id = query_parameters.get('unit_id')       # if not given, value is False
    extension = query_parameters.get('extension')
    output_format = query_parameters.get('format')
    domain = query_parameters.get('domain')
    codon = query_parameters.get('codon')

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

    output = "Basic usage\nformat=[full,unique,top_95_percent,top_10_sequences,multiplicity_5_or_more]\nextension=1 or 2\ndomain=[A,B,C,E,M] comma separated\ncodon=comma separated list\n\n"
    problem = False

    if unit_id:
        output += 'New request made for unit_id %s with extension %s and output format %s\n' % (unit_id,extension,output_format)
        if not "|" in unit_id:
            output += 'Invalid unit id %s\n' % unit_id
            problem = True
        else:
            to_do = [unit_id.split(",")]
    elif loop_id:
        output += 'New request made for loop_id %s with extension %s and output format %s\n' % (loop_id,extension,output_format)
        if not len(loop_id) == 11 or not len(loop_id.split("_")) == 3:
            output += 'Invalid loop id %s\n' % loop_id
            problem = True
        else:
            to_do = [[loop_id]]
    elif id:
        output += 'New request made for id %s with extension %s and output format %s\n' % (id,extension,output_format)
        if not "|" in id and not "_" in id:
            output += 'Invalid id %s\n' % id
            problem = True
        else:
            to_do = [[id]]


    if not problem:
        if output_format in ["full","unique","fasta","top_motif_models"] or "top" in output_format or "multiplicity" in output_format:

            if output_format == "top_motif_models" and not loop_id:
                output += "top_motif_models only available for loop input\n"

            try:
                output_list, families = get_sequence_variability(to_do,extension,output_format,domain,codon)

                if len(output_list) > 30000000:
                    output = output_list[:30000000]
                    output += "\nOutput truncated to 30,000,000 characters.  Restrict the domain or count or codon to get more meaningful output."
                else:
                    # at the moment, output_list is not actually a list
                    output = output_list

            except Exception as inst:
                #output += "%s\n" % type(inst)
                #output += "%s\n" % inst.args
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

    output = "Basic usage; defaults indicated in parentheses\nid=loop id or comma-separated list of unit ids\n&scope=[EC,(Rfam),molecule]\n&resolution=[1.5A,2.0A,2.5A,(3.0A),3.5A,4.0A,20.0A,all]\n&depth=integer (5)\n&format=[(text),json]\n\n"
    problem = False

    to_do = []

    if id:
        output += 'New request made for id %s with scope %s resolution %s depth %d\n' % (id,scope,resolution,depth)
        if not "|" in id and not "_" in id:
            output += 'Invalid id %s\n' % id
            problem = True
        else:
            to_do = [[id]]


    if not problem:
        if scope in ["EC","Rfam","molecule"]:
            if resolution in ['1.5A','2.0A','2.5A','3.0A','3.5A','4.0A','20.0A','all']:

                try:
                    results_list = m_a_s(to_do,scope,resolution,depth)
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