#!/usr/bin/env python
import sys,os,shutil
import optparse

class EmptyFormatter(optparse.IndentedHelpFormatter):
    '''
    This class provide a formmater to avoid optparse format the help text like remove line break, enable the user to control the display of text.
    '''
    def __init__(self,*args,**kwargs):
        optparse.IndentedHelpFormatter.__init__(self,*args,**kwargs)

    #Below parts are modified from optparser.py in python 2.6 lib to make it work
    def format_description(self,description):
        if (description):
            return description
        else:
            return ""

    def format_epilog(self,epilog):
        if ( epilog):
            return epilog
        else:
            return ""


