
'''
import equivalence_class_service as ec
from models import LoopInfo, UnitCorrespondence, UnitInfo, UnitRotation, UnitCenter
from database import db_session
from collections import OrderedDict
from sqlalchemy import case, tuple_
from itertools import groupby
import numpy as np
import time
from discrepancy import matrix_discrepancy
#import utility as ui

def get_loop_range_correspondence(reference_loops, members):

    with db_session() as session:
        
        loop_range_correspondences = OrderedDict()
        for loop, units in reference_loops.iteritems():

            ordering = case (
                       #{unit: index for index, unit in enumerate(units)},
                       dict((unit, index) for index, unit in enumerate(units)),
                       value=UnitCorrespondence.unit_id_1
                       )

            query = session.query(UnitCorrespondence) \
                           .filter(UnitCorrespondence.unit_id_1.in_(units)) \
                           .order_by(ordering) \
                           .filter(tuple_(UnitCorrespondence.pdb_id_2, UnitCorrespondence.chain_name_2) \
                           .in_(members))


            loop_range_correspondences[loop] = [row.unit_id_2 for row in query]


    return loop_range_correspondences


def get_sorted_units(units):
    unsorted_units = units.split(',')
    # This assumes the last element after the split operation to be an integer
    sorted_units = sorted(unsorted_units, key=lambda x: int(x.split('|')[4]))
    return sorted_units 


def group_loop_range_correspondence(loops_range_correspondence):

    grouped_loops_range_correspondence = OrderedDict()
    for loop, correspondence in loops_range_correspondence.iteritems():
        key_ife = lambda x: '|'.join(x.split('|')[:3])
        #sort according to ife
        correspondence.sort(key=key_ife)
        grouped_loops_range_correspondence[loop] = dict((grp, list(items)) for grp, items in groupby(correspondence, key_ife))

    return grouped_loops_range_correspondence


def build_loop_names(loops_range_correspondence):
    
    loop_names = OrderedDict()
    for ref_loop, corresponding_loop in loops_range_correspondence.iteritems():
        inner_range = []
        for chain_id, range_info in corresponding_loop.iteritems():
            pdbid = range_info[0].split("|")[0]
            model = range_info[0].split("|")[1]
            chain = range_info[0].split("|")[2]
            start_pos = range_info[0].split("|")[-1]
            end_pos = range_info[1].split("|")[-1]
            loop_name = model + "/" + chain + "/" + str(start_pos) + ":" + str(end_pos)
            inner_range.append((pdbid, loop_name))

        loop_names[ref_loop] = inner_range

    return loop_names


def get_loop_correspondence(loops_range):
    
    with db_session() as session:

        loop_correspondence = OrderedDict()
        for ref_loop, corresponding_ranges in loops_range.iteritems():
            inner_correspondence = []
            
            query = session.query(LoopInfo) \
                           .filter(tuple_(LoopInfo.pdb_id, LoopInfo.loop_name) \
                           .in_(corresponding_ranges))  

            for row in query:
                loop_units = row.unit_ids
                # Here we are sorting based on the residue num. This will fail
                # if the residue_num has insertions etc
                sorted_loop_units = get_sorted_units(loop_units)
                inner_correspondence.append((row.loop_id, sorted_loop_units))           

            loop_correspondence[ref_loop] = inner_correspondence
        
    return loop_correspondence


def get_loop_rotation(loop_correspondence):
    
    with db_session() as session:

        rotation_data_collection = OrderedDict()
        for ref_loop, correspondences in loop_correspondence.iteritems():
            rotation_dict = {}
            for correspondence in correspondences:

                rotation_data = []
                rotation_data2 = []
                
                loop_id = correspondence[0]
                units = correspondence[1]

                ordering = case (
                           #{unit: index for index, unit in enumerate(units)},
                           dict((unit, index) for index, unit in enumerate(units)),
                           value=UnitRotation.unit_id
                           )

                query = session.query(UnitRotation) \
                               .filter(UnitRotation.unit_id.in_(units)) \
                               .order_by(ordering)

                for row in query:
                    rotation_data.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2],
                                                   [row.cell_1_0, row.cell_1_1, row.cell_1_2],
                                                   [row.cell_2_0, row.cell_2_1, row.cell_2_2]]))

                    rotation_data2.append(row.unit_id)
                    rotation_dict[loop_id] = rotation_data  

            rotation_data_collection[ref_loop] = rotation_dict
            
        return rotation_data_collection 


def get_reference_rotation(reference_loops):

    with db_session() as session:

        rotation_dict = OrderedDict()
        for loop_id, units in reference_loops.iteritems():

            rotation_values = []
            rotation_values2 = []

            ordering = case (
                           #{unit: index for index, unit in enumerate(units)},
                           dict((unit, index) for index, unit in enumerate(units)),
                           value=UnitRotation.unit_id
                           )

            query = session.query(UnitRotation) \
                               .filter(UnitRotation.unit_id.in_(units)) \
                               .order_by(ordering)

            for row in query:
                    rotation_values.append(np.array([[row.cell_0_0, row.cell_0_1, row.cell_0_2],
                                                   [row.cell_1_0, row.cell_1_1, row.cell_1_2],
                                                   [row.cell_2_0, row.cell_2_1, row.cell_2_2]]))
                    
                    rotation_values2.append(row.unit_id)
                    rotation_dict[loop_id] = rotation_values2

    return rotation_dict



def get_reference_center(reference_loops):

    with db_session() as session:

        center_dict = OrderedDict()
        for loop_id, units in reference_loops.iteritems():

            center_values = []
            center_values2 = []

            ordering = case (
                           #{unit: index for index, unit in enumerate(units)},
                           dict((unit, index) for index, unit in enumerate(units)),
                           value=UnitCenter.unit_id
                           )

            query = session.query(UnitCenter) \
                           .filter(UnitCenter.name == "base") \
                           .filter(UnitCenter.unit_id.in_(units)) \
                           .order_by(ordering)

            for row in query:
                    center_values.append(np.array([row.x, row.y, row.z]))
                    
                    center_values2.append(row.unit_id)
                    center_dict[loop_id] = center_values2

    return center_dict


def get_loop_center(loop_correspondence):
    
    with db_session() as session:

        center_data_collection = OrderedDict()
        for ref_loop, correspondences in loop_correspondence.iteritems():
            center_dict = {}
            for correspondence in correspondences:

                center_data = []
                center_data2 = []
                
                loop_id = correspondence[0]
                units = correspondence[1]

                ordering = case (
                           #{unit: index for index, unit in enumerate(units)},
                           dict((unit, index) for index, unit in enumerate(units)),
                           value=UnitCenter.unit_id
                           )

                query = session.query(UnitCenter) \
                               .filter(UnitCenter.name == 'base') \
                               .filter(UnitCenter.unit_id.in_(units)) \
                               .order_by(ordering)

                for row in query:
                    center_data.append(np.array([row.x, row.y, row.z]))

                    center_data2.append(row.unit_id)
                    center_dict[loop_id] = center_data2 

            center_data_collection[ref_loop] = center_dict
            
        return center_data_collection       


def get_reference_loops(chain):

    pdb_id, _, chain_id = chain.split("|")

    with db_session() as session:
        query = session.query(LoopInfo) \
                       .filter(LoopInfo.pdb_id == pdb_id) \
                       .filter(LoopInfo.type == 'HL')

        loop_reference = OrderedDict()
        for row in query:
            position = row.loop_name
            units = row.unit_ids
            model, chain_name, loop_range = position.split("/")
            if chain_name == chain_id:
                range_start_pos, range_end_pos = loop_range.split(":")
                loop_reference[row.loop_id] = (range_start_pos, range_end_pos, row.unit_ids)

    return loop_reference


def format_loop_ranges(reference_loops):
    
    formatted_loop_ranges = OrderedDict()
    reference_loop_units = OrderedDict()
    for loop_id, range_info in reference_loops.iteritems():
        start_pos = range_info[0]
        end_pos = range_info[1]
        units = range_info[2]
        sorted_units = get_sorted_units(units)
        start_unit = None
        end_unit = None
        for unit in sorted_units:
            residue_num = unit.split('|')[-1]
            if residue_num == start_pos:
                start_unit = unit
            if residue_num == end_pos:
                end_unit = unit
        formatted_loop_ranges[loop_id] = [start_unit, end_unit]
        reference_loop_units[loop_id] = sorted_units
        
    return formatted_loop_ranges, reference_loop_units


def calculate_geometric_disc(ref_rotation, ref_center, rotation_data, center_data):
    disc_dict = {}
    
    for k1, v1 in ref_rotation.iteritems():
        for k2, v2 in rotation_data.iteritems():
            for k3, v3 in v2.iteritems():
                print k1 + " " + k3

    return None


chain_id = '5J7L|1|AA'

start = time.time()

reference_loops = get_reference_loops(chain_id)

reference_loops, reference_loops_units = format_loop_ranges(reference_loops)

# Get the equivalence class members that the query ife belongs to
members, ec_name, nr_release = ec.get_ec_members('4.0', 'X-RAY DIFFRACTION', chain_id)

loop_range_correspondences = get_loop_range_correspondence(reference_loops, members)

grouped_loops_range_correspondence = group_loop_range_correspondence(loop_range_correspondences)

loop_names = build_loop_names(grouped_loops_range_correspondence)

loop_correspondences = get_loop_correspondence(loop_names)

rotation_data = get_loop_rotation(loop_correspondences)

center_data = get_loop_center(loop_correspondences)

reference_rotation = get_reference_rotation(reference_loops_units)

reference_center = get_reference_center(reference_loops_units)

disc_data = calculate_geometric_disc(reference_rotation, reference_center,
                                     rotation_data, center_data)

end = time.time() 

time_diff = '{0:.2f}'.format(end-start)

print disc_data
'''
