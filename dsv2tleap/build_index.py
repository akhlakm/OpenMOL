#!/usr/bin/env python3

## Usage:
def print_usage():
	print("%s <unit mol2 1> [<unit mol2 2>] [<unit mol2 3>] ..." %sys.argv[0])
	exit(1)


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

import sys
sys.path.append("..")

# Sanity check
if len(sys.argv) < 2:
	print_usage()

try:
	import openmol
except ImportError:
	print("Failed to import openmol. Please add to PYTHONPATH, or update sys.path.")
	raise

import tripos_mol2 as mol2 

# Use hash table pointing each other
# for better searching
pc_orig = {}
orig_pc = {}

# Initial index value, increment
index_val 	= 5.0000
d_index 	= 0.0001

# Counter
count = 0

## MOL2 processing
## -----------------------------------------------
for molfile in sys.argv[1:]:
	say("Processing %s" %molfile)
	unit = mol2.read(molfile)
	count = 0
	
	say("Building pc_index")
	for i, a_name in enumerate(unit['atom_name']):
		a_type = unit['atom_type'][i]
		a_charge = unit['atom_q'][i]
		a_resname = unit['atom_resname'][i]

		orig_key = "%s;%.4f" %(a_type, a_charge)

		if not orig_key in orig_pc:
			index_val += d_index
			orig_pc[orig_key] = "%.4f" %index_val
			orig_items = [a_type, a_charge, a_resname]
			pc_orig["%.4f" %index_val] = orig_items
			count += 1
			if len(pc_orig) > 9999:
				alert("WARN! Too many unique atom types.")

		new_charge = orig_pc[orig_key]

		# update charge
		unit['atom_q'][i] = float(new_charge)

	done("Added %d new indices" %count)
	done()
	unit['title'] = "PC_index MOL2 with partial charge as reference."
	unit = mol2.build(unit)
	mol2.Writer(unit, 'dsv_'+molfile.split("/")[-1]).write()

for pc in pc_orig:
	print(pc, " <-- ", pc_orig[pc])

say("Saving pc_orig")
openmol.write_json(pc_orig, "pc_orig.json")
done()

say("Saving orig_pc")
openmol.write_json(orig_pc, "orig_pc.json")
done()

print("Please use the written dsv_ MOL2 files for system building.")
print("Then use ./apply_index.py <system.mol2> to apply original types and charges.")
