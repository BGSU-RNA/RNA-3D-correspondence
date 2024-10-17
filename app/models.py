import logging
logging.info("Importing models 1")

from sqlalchemy.sql.schema import Column
logging.info("Importing models 2")
from sqlalchemy import String, Table
logging.info("Importing models 3")
from database import Base
logging.info("Importing models 4")

# This is a view. Had to make each col as primary key for the query to work
class UnitCorrespondence(Base):
    __table__ = Table('correspondence_units', Base.metadata,
                Column('unit_id_1', String, primary_key=True),
                Column('chain_name_1', String, primary_key=True),
                Column('pdb_id_1', String, primary_key=True),
                Column('unit_id_2', String, primary_key=True),
                Column('chain_name_2', String, primary_key=True),
                Column('pdb_id_2', String, primary_key=True),
                )


class UnitRotation(Base):
    __table__ = Base.metadata.tables['unit_rotations']


class UnitCenter(Base):
    __table__ = Base.metadata.tables['unit_centers']


class NrChains(Base):
    __table__ = Base.metadata.tables['nr_chains']


class NrClassRank(Base):
    __table__ = Base.metadata.tables['nr_class_rank']


class NrClasses(Base):
    __table__ = Base.metadata.tables['nr_classes']


class NrReleases(Base):
    __table__ = Base.metadata.tables['nr_releases']


class UnitInfo(Base):
    __table__ = Base.metadata.tables['unit_info']


class IfeInfo(Base):
    __table__ = Base.metadata.tables['ife_info']


class PDBInfo(Base):
    __table__ = Base.metadata.tables['pdb_info']


class LoopInfo(Base):
    __table__ = Base.metadata.tables['loop_info']


class LoopPositions(Base):
    __table__ = Base.metadata.tables['loop_positions']


class UnitPairsInteractions(Base):
    __table__ = Base.metadata.tables['unit_pairs_interactions']


class UnitPairsInteractions2024(Base):
    __table__ = Base.metadata.tables['unit_pairs_interactions_2024']


# class UnitPairsDistances(Base):
#     __table__ = Base.metadata.tables['unit_pairs_distances']


class ChainInfo(Base):
    __table__ = Base.metadata.tables['chain_info']


class ChainPropertyValue(Base):
    __table__ = Base.metadata.tables['chain_property_value']


# class ProteinChainProperty(Base):
#     __table__ = Base.metadata.tables['protein_chain_property']


