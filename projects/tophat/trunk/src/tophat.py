#!/usr/bin/env python
# encoding: utf-8
"""
tophat.py

Created by Cole Trapnell on 2008-12-25.
Copyright (c) 2008 Cole Trapnell. All rights reserved.
"""

import sys
import getopt


help_message = '''
TopHat maps short sequences from spliced transcripts to whole genomes.
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				output = value
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2


if __name__ == "__main__":
	sys.exit(main())
