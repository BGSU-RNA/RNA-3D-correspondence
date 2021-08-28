from collections import OrderedDict
from models import UnitInfo, LoopInfo
from database import db_session
from sqlalchemy import or_
import utility as ui

def get_units(incomplete_units):

    with db_session() as session:    
        
        complete_units = []
        for unit in incomplete_units:
            query = session.query(UnitInfo).filter(UnitInfo.unit_id.like(unit))
            for row in query:
                complete_units.append(row.unit_id)
    
        return complete_units


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
        complete_units = selection.split(",")

    elif input_type == 'res_num':
        incomplete_units = selection.split(",")
        chain = chain_id + '|%|'
        incomplete_units = [chain + unit for unit in incomplete_units]
        complete_units = get_units(incomplete_units)


    return complete_units


