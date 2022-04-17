from collections import OrderedDict
from models import UnitInfo, LoopInfo
from database import db_session
from sqlalchemy import or_
import utility as ui
import itertools
import re

def get_units(incomplete_units):

    with db_session() as session:    
        
        complete_units = []
        for unit in incomplete_units:
            query = session.query(UnitInfo).filter(UnitInfo.unit_id.like(unit))
            for row in query:
                complete_units.append(row.unit_id)
    
        return complete_units

def get_single_range_units(range_positions, pdb_id, chain):

    with db_session() as session:

        complete_units = []
        query = session.query(UnitInfo).filter_by(pdb_id=pdb_id) \
                                       .filter_by(chain=chain) \
                                       .filter(UnitInfo.chain_index.between(range_positions[0], range_positions[1])) \
                                       .order_by(UnitInfo.chain_index)

        for row in query:
            complete_units.append(row.unit_id)

        return complete_units


def get_multiple_range_units(range_positions_list, pdb_id, chain):

    with db_session() as session:
        
        complete_units = []
        for single_range in range_positions_list:
            single_range_units = []
            query = session.query(UnitInfo).filter_by(pdb_id=pdb_id) \
                                       .filter_by(chain=chain) \
                                       .filter(UnitInfo.chain_index.between(single_range[0], single_range[1])) \
                                       .order_by(UnitInfo.chain_index)

            for row in query:
                single_range_units.append(row.unit_id)

            complete_units.append(single_range_units)

        #merged_units = list(itertools.chain(*complete_units))
        # return the flattened list of lists
        return list(itertools.chain(*complete_units))

def get_loop_units(loop_id):
    
    with db_session() as session:

        query = session.query(LoopInfo).filter_by(loop_id=loop_id)

        return query[0].unit_ids

def get_query_units(query_type, query_list, query_ife):

    if query_type == 'units_str':
        incomplete_units = [unit[0] for unit in query_list]
        ife = query_ife + '|%|'
        incomplete_units = [ife + unit for unit in incomplete_units]
        complete_units = get_units(incomplete_units)

    elif query_type == 'single_range':
        pass

    elif query_type == 'multiple_ranges':
        pass
        
    elif query_type == 'loop_id':
        loop_id = query_list[0][0]
        unsorted_units = get_loop_units(loop_id)
        complete_units = ui.get_sorted_units(unsorted_units)

    return complete_units


def get_query_units_new(input_type, selection, chain_id):

    if input_type == 'loop_id':
        loop_id = selection
        unsorted_units = get_loop_units(loop_id)
        complete_units = ui.get_sorted_units(unsorted_units)

    elif input_type == 'unit_id':
        # complete_units = selection.split(",")
        complete_units = re.split(',|\t', selection)

    elif input_type == 'res_num':
        incomplete_units = selection.split(",")
        chain = chain_id + '|%|'
        incomplete_units = [chain + unit for unit in incomplete_units]
        complete_units = get_units(incomplete_units)

    elif input_type == 'single_range':
        pdb_id, model, chain = chain_id.split('|')
        range_positions = selection.split(':')
        complete_units = get_single_range_units(range_positions, pdb_id, chain)

    elif input_type == 'multiple_ranges':
        pdb_id, model, chain = chain_id.split('|')
        ranges_selection = selection.split(",")
        range_positions_list = [single_range.split(':') for single_range in ranges_selection]
        complete_units = get_multiple_range_units(range_positions_list, pdb_id, chain)

    return complete_units


