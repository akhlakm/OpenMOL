#!/usr/bin/env python3

import sys
import json
sys.path.append("..")

try:
	import openmol
except ImportError:
	print("Failed to import openmol. Please add to PYTHONPATH, or update sys.path.")
	raise

import tripos_mol2 as mol2 

def print_usage():
	print("%s <DSV system mol2>" %sys.argv[0])
	exit(1)

def load_json(json_file):
	with open(json_file, 'r') as fp:
		obj = json.load(fp)
	return obj

## Simple Logging
## -----------------------------------------------
def alert(msg):
	print(" -- %s\n" %(msg))

nest_index = 0

def say(msg, alert=False):
	global nest_index
	nest_index += 1
	nest = nest_index * " --"
	print("%s %s" %(nest, msg))

def done(msg=""):
	global nest_index
	nest = nest_index * " --"
	print("%s OK. %s" %(nest, msg))
	if nest_index > 0: nest_index -= 1

## Initialize
## -----------------------------------------------

# Use hash table pointing each other
# for better searching
pc_orig = load_json('pc_orig.json')
orig_pc = load_json('orig_pc.json')

unknown_pc_index = []

## MOL2 processing
## -----------------------------------------------
molfile = sys.argv[1]
unit = mol2.read(molfile)
count = 0

say("Restoring pc_index")
for i, a_name in enumerate(unit['atom_name']):
	pc_index 	= "%.4f" %unit['atom_q'][i]

	if pc_index in pc_orig:
		count += 1
		unit['atom_type'][i] = pc_orig[pc_index][0]
		unit['atom_q'][i] = pc_orig[pc_index][1]
	else:
		unknown_pc_index.append(pc_index)
		alert("Unknown pc_index: %s atom ID %d" %(pc_index, i+1))

	## Do additional processing here if needed
	## -----------------------------------------------
	# unit['atom_resname'][i] = 'LIG'


done("Applied pc_index to %d atoms." %count)

if len(unknown_pc_index) > 0:
	alert("Please update the script to handle the unknown atoms or update manually.")
else:
	done("All pc_index applied successfully.")

unit['title'] = "PC index applied DSV system for tleap"
unit = mol2.build(unit)
mol2.Writer(unit, 'tleap_'+molfile.split("/")[-1]).write()

print("Please use the written tleap_dsv_ MOL2 files for system building in tleap.")
