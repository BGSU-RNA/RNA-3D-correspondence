from sqlalchemy import create_engine
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker, scoped_session
from contextlib import contextmanager

'''
engine = create_engine(SQLALCHEMY_DATABASE_URL)
session = sessionmaker(autocommit=False, autoflush=False, bind=engine)
'''

with open("/usr/local/pipeline/flask.config","r") as f:
	t = f.readlines()

database_connection_string = t[0].replace("\n","")

# pool_pre_ping=True was tried but did not work. it would be used to check the connection before using it; maybe that will reduce Internal Service Errors

engine = create_engine(database_connection_string, convert_unicode=True, echo=False)
Base = declarative_base()
Base.metadata.reflect(engine)

SessionLocal = scoped_session(sessionmaker(bind=engine))

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
