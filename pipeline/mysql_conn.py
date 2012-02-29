#! /usr/bin/python2


import os
import MySQLdb

class MySQL_conn:
    def __init__(self, hostname, database, username, password):
        
        
        self.connect = MySQLdb.connect(host   = hostname,
                               db     = database,
                               user   = username,
                               passwd = password)
        self.cursor = self.connect.cursor()
    

    def connection(self):
        return self.connect


