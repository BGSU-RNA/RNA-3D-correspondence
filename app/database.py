import logging
logging.info("Importing database 1")
import os

# the next two lines are for the new server in 2024
import pymysql
pymysql.install_as_MySQLdb()

logging.info("Importing database 2")
# new for 2024, add MetaData
from sqlalchemy import create_engine, MetaData
logging.info("Importing database 3")
from sqlalchemy.ext.declarative import declarative_base
logging.info("Importing database 4")
from sqlalchemy.orm import sessionmaker, scoped_session
logging.info("Importing database 5")
from contextlib import contextmanager
logging.info("Importing database 6")


with open("/usr/local/pipeline/flask.config","r") as f:
	t = f.readlines()
logging.info("Importing database 7")

database_connection_string = t[0].replace("\n","")

# logging.info("database_connection_string %s" % database_connection_string)
logging.info(str(os.environ))
logging.info("Importing database 8")

engine = create_engine(database_connection_string, convert_unicode=True, echo=False)
Base = declarative_base()

# lines below are new for 2024
metadata = MetaData()
logging.info("Importing database 9a")
try:
    metadata.reflect(bind=engine)
except Exception as e:
    logging.exception("Error reflecting metadata %s" % e)
logging.info("Importing database 9b")
Base.metadata = metadata

logging.info("Importing database 9c")

# old line
# SessionLocal = scoped_session(sessionmaker(bind=engine))

# new line for 2024
SessionLocal = scoped_session(sessionmaker(autocommit=False, autoflush=False, bind=engine))

logging.info("Importing database 10")

@contextmanager
def db_session():
	db = SessionLocal()
	try:
		yield db
		db.commit()
	except:
		db.rollback()
		raise
	finally:
		db.close()
