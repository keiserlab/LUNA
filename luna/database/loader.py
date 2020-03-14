# import MySQLdb
from sqlalchemy import create_engine, MetaData, Table
from sqlalchemy.orm import mapper, sessionmaker
from sqlalchemy import inspect

from luna.util.exceptions import IllegalArgumentError

import logging
logger = logging.getLogger()


class DBLoader():

    def __init__(self, user, passwd, host, port, db):
        self._user = user
        self._passwd = passwd
        self._host = host
        self._db = db
        self._port = int(port)

        self._engine = None
        self._metadata = None
        self._session = None

        self._start()

    @property
    def session(self):
        self.new_session()
        return self._session

    @property
    def metadata(self):
        return self._metadata

    def _start(self):
        if (self._engine is None):
            try:
                # Engine: mysql-python
                self._engine = create_engine('mysql+mysqldb://%s:%s@%s:%i/%s' %
                                             (self._user, self._passwd, self._host, self._port, self._db))
                self._metadata = MetaData(self._engine)
            except Exception as e:
                logger.exception(e)
                raise

    def get_table(self, tbl_name, autoload=True):
        return Table(tbl_name, self._metadata, autoload=autoload)

    def new_mapper(self, trgt_class, table, properties=None):
        # Return if already mapped
        if inspect(trgt_class, raiseerr=False):
            return

        if isinstance(table, str):
            table = self.get_table(tbl_name=table, autoload=True)

        if not isinstance(table, Table):
            raise IllegalArgumentError("Parameter 'table' should be a table name or a Table object.")

        # TODO: check if table exists
        if (properties is None):
            mapper(trgt_class, table)
        else:
            mapper(trgt_class, table, properties=properties)

    def new_session(self):
        if (self._session is None):
            Session = sessionmaker(bind=self._engine)
            self._session = Session()

    def remove_session(self):
        self._session = None

    def approve_session(self):
        try:
            self._session.commit()
        except Exception as e:
            logger.exception(e)
        finally:
            self.remove_session()
