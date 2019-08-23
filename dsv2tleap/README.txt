
Prepare MOL2 file(s) for system building in DSV which later can be
imported in tleap.

1. Prepare the units:
	Use ./build_index.py on all the MOL2 files to set partial charges
	replaced by some indices (pc_index).
	This will save the MOL2 files with dsv_ prefix.
	The index will also be saved in two json files.

2. Load dsv_*.mol2 files in DSV:
	do editing
	add hydrogens, they should have 0.000 partial charge.
	copy paste units, the pc index should be kept the same.
	save as the final system mol2 file. (single)


3. Restore original types and charges:
	Use ./apply_index.py <system.mol2>
	It will parse the atoms, match the partial charge with the saved
	indices in the pc_orig.json file and restore
	the original types and partial charges.
	If you need additional updates, for example setting charges of the
	newly added H atoms, you can edit apply_index.py for loop.
	The final structure will be saved with a tleap_dsv_ prefix.
	Now you can load the MOL2 in tleap to solvate or make parm7 files.

NOTE: 	pc_index uses partial charges >5.0000 for maintaining index.
	You can add manual partial charges if they are not >5.0,
	which is very unlikely.

