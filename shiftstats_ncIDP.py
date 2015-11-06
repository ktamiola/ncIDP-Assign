# ncIDP-Library Based Repositioning Plugin for Sparky NMR Assignment Software
# Authors:
#			Kamil Tamiola (k.tamiola@rug.nl)
#			Frans A.A. Mulder (f.a.a.mulder@rug.nl)

import string
import sys
import os

import sputil
import atomnames
import sequence
import sparky
import pyutil

class atom_stats:
  """
  	****************************
	v 1.0 | 11th March 2010 | kT
	****************************
	Atom Statistics Class adopted from Sparky:
	  gsym            -> amino acid one letter code
	  atom_name       -> type of nucleus
	  group_numerb    -> sequence number
	  average_shift   -> sequence specific random coil shift
	  shift_deviation -> sequence specific chemical shift deviation
  """
  def __init__(self, gsym, aname, gnumber, fave, fdev):

    self.group_symbol = gsym
    self.atom_name = aname
    self.group_number = gnumber
    self.average_shift = fave
    self.shift_deviation = fdev


# -----------------------------------------------------------------------------


def get_sparky_home():
  """
  Get the SPARKY HOME variable
  """
  if os.environ.has_key('HOME'):
    # Check if we have SPARKYHOME
    if os.path.exists(os.environ['HOME']+'/Sparky'):
      return (os.environ['HOME']+'/Sparky')
  elif os.environ.has_key('SPARKYHOME'):
    return os.environ['SPARKYHOME']

def atom_statistics(group_number, group_symbol, atom_name, ast):
  key = (group_number, group_symbol, atom_name)
  if ast.has_key(key):
    return ast[key]
  
  return None

def resonance_statistics(reslist):
  """
    Compute global, sequence-specific ncIDP-based chemical shifts using resonance list
  """
  r = reslist[0]
  seq = sequence.molecule_sequence(r.condition.molecule).one_letter_codes
  path = get_sparky_home()+'/Python'
  prediction = compute_RC(path+'/ncIDP-cs.tab', path+'/ncIDP-dev.tab', seq)
  return prediction

def single_resonance_statistics(r):
  """
    Compute global, sequence-specific ncIDP-based chemical shifts using resonance list
  """
  seq = sequence.molecule_sequence(r.condition.molecule).one_letter_codes
  path = get_sparky_home()+'/Python'
  prediction = compute_RC(path+'/ncIDP-cs.tab', path+'/ncIDP-dev.tab', seq)
  return prediction

def sequence_resonance_statistics(tmp):
  """
    Compute global, sequence-specific ncIDP-based chemical shifts using FASTA sequence
  """
  seq = tmp.one_letter_codes
  prediction = compute_RC('ncIDP-cs.tab', 'ncIDP-dev.tab', seq)
  return prediction

def read_File(filename):
	"""
		****************************
		v 1.0 | 11th March 2010 | kT
		****************************
		Read an ACSII file and return a list o lines without the new-line character
	"""
	star_data=open(filename).readlines()						# Start reading file-by-file
	raw_data=[]													# Declare empty lists, for future use
	separated=[]
	for p in star_data[0:]:										# Delete '\n' new-line signs from list 
		raw_data.append(p[:-1])
	for p in range (len(raw_data)):								# Separate elements
		separated.append(raw_data[p].split())	
	return separated

def read_ncIDP_table(fi):
	"""
		****************************
		v 1.0 | 11th March 2010 | kT
		****************************
		Read ncIDP chemical shift table and return structured dictionary
		Important:
		  - anything after REMARK keyword is ignored
		  - first column contains amino acid type
		  - there are 18 columns organized in triplets
		  	- within every triplet column 1 is the correction for the left-neighbour
			- 2nd is the value for the central residue
			- 3rd is the correction for the right-neighbour
		  - triplet order corresponds to:
		      H HA C CA CB N
		    nuclei order! Keep it that way! Otherwise modify cs(list) within read_ncIDP_table function!
	"""
	cs = ['H','HA','C','CA','CB','N']
	table = {}
	data = read_File(fi)
	for l in data:
		if not l[0]=='REMARK':
			table[l[0]]={}
			for s in range(len(cs)):
				table[l[0]][cs[s]]={}
				table[l[0]][cs[s]]['p']=float(l[1+s*3])
				table[l[0]][cs[s]]['c']=float(l[2+s*3])
				table[l[0]][cs[s]]['n']=float(l[3+s*3])
	return table

atom_stat_table = None
def compute_RC(fi_cs, fi_dev, seq):
	"""
		****************************
		v 1.0 | 11th March 2010 | kT
		****************************
		Compute neighbor-corrected random coil chemical shifts and their deviations for an arbitrary 
		protein sequence using ncIDP library
	"""
	cs = ['H','HA','C','CA','CB','N']
	RC={}
	ncIDP_shifts = read_ncIDP_table(fi_cs)
	ncIDP_devs = read_ncIDP_table(fi_dev)
	for res in range(1,len(seq)+1):
		i = res - 1 # This is the actual python index
		RC[res]={}
		RC[res]['symbol']=seq[i]
		for n in cs: 
			RC[res][n]={}
			if i!=0 and i<len(seq)-1:
				RC[res][n]['cs']=ncIDP_shifts[seq[i-1]][n]['p']+ncIDP_shifts[seq[i]][n]['c']+ncIDP_shifts[seq[i+1]][n]['n']
				RC[res][n]['dev']=ncIDP_devs[seq[i-1]][n]['p']+ncIDP_devs[seq[i]][n]['c']+ncIDP_devs[seq[i+1]][n]['n']
			else:
				if i==0:
					RC[res][n]['cs']=ncIDP_shifts[seq[i]][n]['c']+ncIDP_shifts[seq[i+1]][n]['n']
					RC[res][n]['dev']=ncIDP_devs[seq[i]][n]['c']+ncIDP_devs[seq[i+1]][n]['n']
				else:
					RC[res][n]['cs']=ncIDP_shifts[seq[i-1]][n]['p']+ncIDP_shifts[seq[i]][n]['c']
					RC[res][n]['dev']=ncIDP_devs[seq[i-1]][n]['p']+ncIDP_devs[seq[i]][n]['c']				
	# Now change the structure of the dictionary so we can directly inject it into Sparky
	global atom_stat_table
	if atom_stat_table == None:
		atom_stat_table = {}
	for resid in RC.keys():
		for n in cs:
			if (seq[resid-1]=='G' and n=='CB') or (seq[resid-1]=='P' and n=='H'): continue
                        gsym = RC[resid]['symbol']
			gnumber = resid
			aname = n
			fave = RC[resid][n]['cs']
			fdev = RC[resid][n]['dev']
			atom_stat_table[(resid, RC[resid]['symbol'], n)]=atom_stats(gsym, aname, gnumber, fave, fdev)
	return atom_stat_table
