#!/usr/bin/env python
#--coding:utf-8--
"""
utils.py
Utilities for cLoops2
"""

__author__ = "CAO Yaqiang"
__email__ = "caoyaqiang0410@gmail.com"

#sys library
import os
import sys
import time
import gzip
import logging
import argparse

#glob settings
#epilog for argparse
EPILOG = "Any bug is welcome reported to caoyaqiang0410@gmail.com"


def isTool(name):
    """
    Check if a tool is on PATh and marked as executable.

    Parameters
    ---
    name: str

    Returns
    ---
    True or False
    """
    from distutils.spawn import find_executable
    if find_executable(name) is not None:
        return True
    else:
        return False


def timer(func, REPEATS=5):
    """
    TIMER for estimate the running time of a funciton, can be used as decorator.

    Parameters
    ---
    func: funciton
    REPEATS: int, repeat time to run target function. 

    Usage
    ---
    @timer
    def run():
        for i in range(1,10000):
            pass
    """

    def wrapper(*args, **kwargs):
        start = time.time()
        for _ in range(REPEATS):
            v = func(*args, **kwargs)
        t = (time.time() - start) / REPEATS
        print("{} time elapsed: {}".format(func.__name__, t))
        return v

    return wrapper


def getLogger(fn=os.getcwd() + "/" + os.path.basename(__file__) + ".log"):
    """
    Setting up the logger systems.

    Parameters
    ----
    fn: str, file name to store the logging infor, default is genreated with time and script name

    Returns
    ----
    logging.loger 
    """
    #get the current time
    date = time.strftime(' %Y-%m-%d', time.localtime(time.time()))
    #set up logging, both write log info to console and log file
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(name)-6s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        filename=fn,
        filemode='a')
    logger = logging.getLogger(
        "cLoops2"
    )  #here "cLoops2" or anyother name can supress the logging (DEBUG) information from 3rd package
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.NOTSET)
    return logger


def callSys(cmds, logger=None):
    """
    Call systematic commands without return.

    Parameters
    ---
    cmds: list, commands to run in the bash
    logger: logging.loger
    """
    for c in cmds:
        if logger is not None:
            logger.info(c)
        else:
            print(c)
        try:
            os.system(c)
        except:
            if logger is not None:
                logger.error(c)
            else:
                print("ERROR for cmd %s" % c)


def cFlush(r):
    """
    One line flush to show the programming process.

    Parameters
    ---
    r: str
    """
    sys.stdout.write("\r%s" % r)
    sys.stdout.flush()
