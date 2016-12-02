#!/usr/bin/python                                                                                               

from optparse import OptionParser

__version__ = "0.1"
__date__ = "26 Mar 2012"
__maintainer__ = "Ellen Leffler"
__email__ = "emleffler@uchicago.edu"

def parse_args():
  """parses option arguments.  ensures there are enough required arguments.
  exits and reports an error if too vew arguments, or incorrect options.                                                                                            
  """

  usage = "%prog -r reference_file -b bed_file"
  description = "".join(["from a bed file, ",
                         "pulls out the sequence from a reference",
                         "for those positions."])
  version = "%%prog %s (released %s)" % (__version__, __date__)

  parser = OptionParser(usage=usage, description=description, version=version)
  parser.add_option("-r","--reference_file", dest="reference_file",
                    help="Specify file with reference sequence")
  parser.add_option("-b","--bed_file", dest="bed_file",
                   help="Specify files with columns chr, position1 and position2")
  return parser.parse_args()

def load_pseudo_fasta_file(fname):
  if fname[-3:] == '.gz':
    f = gzip.open(fname, "r")
  else: 
    f = open(fname, "r")
  header = f.readline().strip()
  values = "".join([line.strip() for line in f])
  return (header, values)

options, args = parse_args()

#load reference
reference = load_pseudo_fasta_file(options.reference_file)
sites = open(options.bed_file)
#No header, exclude these lines:
#header=sites.readline().split()
#print header
#print "\t".join(header)

for pos in sites:
  pos = pos.split()
  bpos=int(pos[1])-1
  epos=int(pos[2])
  seqint = reference[1][bpos:epos]
  pos.append(seqint.upper())
  print "\t".join(pos)
