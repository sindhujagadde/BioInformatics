filename = input('Enter input file name: ')
pdb = list(open(filename))
sheets = list()
locs = dict()
atoms = dict()

for line in pdb:
	data = line.split()
	if not data:
		continue
	if data[0] == "SHEET":
		#this is a beta sheet
		start = data[6]
		end = data[9]
		sheets.append((start, end))
	elif data[0] == "ATOM" and data[2] == "CA":
		if str(data[5]) not in locs:
			locs[str(data[5])] = data[6:9]
		if str(data[5]) not in atoms:
			atoms[str(data[5])] = data[3]
for sheet in sheets:

	
	print("NAMES: ", atoms[str(sheet[0])], atoms[str(sheet[1])])
	print("Indices: ", sheet)
	print("Vector: ", locs[str(sheet[0])], locs[str(sheet[1])])
	