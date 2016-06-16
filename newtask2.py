import numpy as np
import sys
from math import sin, cos, pi, acos
class BetaSheet:
	def __init__(self, start, end):
		self.start = start
		self.end = end
	def load(self, atoms):
		self.atoms = atoms
		CA1 = atoms[self.start]
		CA3 = atoms[str(int(self.start) + 2)]
		self.vector = (CA1.loc, CA3.loc)
		self.direction = np.array(CA3.loc) - np.array(CA1.loc)
	def calculate_transforms(self, other):
		translation = np.array(self.vector[0]) + np.array(other.vector[0])
		rotation_axis = np.cross(self.direction, other.direction)
		with np.errstate(invalid='ignore'):
			rotation_axis_normed = rotation_axis / np.linalg.norm(rotation_axis)
			cos_angle = self.direction.dot(other.direction) / np.linalg.norm(self.direction) / np.linalg.norm(other.direction)
			angle = acos(cos_angle)
		return translation, (rotation_axis, angle)
	def move_to(self, other):
		translation, rotation = self.calculate_transforms(other)
		for index in range(int(self.start), int(self.end) + 1):
			atom = self.atoms[str(index)]
			atom.translate(translation)
			atom.rotate(*rotation)
class Atom:
	def __init__(self, line):
		tokens = line.split()
		self.tokens = tokens
		if tokens[0] != "ATOM":
			self.tokens = []
			return None
		self.index = int(tokens[5])
		try:
			self.loc = (float(tokens[6]), float(tokens[7]), float(tokens[8]))
		except ValueError:
			line = line.replace('-', ' -')
			Atom.__init__(self, line)
	def translate(self, transform):
		self.loc = tuple(np.array(self.loc) + np.array(transform))
	def rotate(self, axis, angle):
		x, y, z = axis
		point = self.loc
		I = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
		CP = np.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
		Tensor = np.matrix([[x * x, x * y, x * z], [x * y, y * y, z * y], [x * z, y * z, z * z]])
		rotation = cos(angle) * I + sin(angle) * CP + (1 - cos(angle)) * Tensor
		point_vector = np.matrix([[point[0]], [point[1]], [point[2]]])
		self.loc = [x[0] for x in (rotation * point_vector).tolist()]
	def __repr__(self):
		self.tokens[6] = str(self.loc[0])
		self.tokens[7] = str(self.loc[1])
		self.tokens[8] = str(self.loc[2])
		return " ".join(self.tokens)
def translate(transform, point):
	xt, yt, zt = transform
	xp, yp, zp = point
	xp += xt
	yp += yt
	zp += zt
	return (xp, yp, zp)
def rotate(direction, angle, point):
	x, y, z = direction
	I = np.matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
	CP = np.matrix([[0, -z, y], [z, 0, -x], [-y, x, 0]])
	Tensor = np.matrix([[x * x, x * y, x * z], [x * y, y * y, z * y], [x * z, y * z, z * z]])
	rotation = cos(angle) * I + sin(angle) * CP + (1 - cos(angle)) * Tensor
	point_vector = np.matrix([[point[0]], [point[1]], [point[2]]])
	return (rotation * point_vector).transpose()

filename = input('Enter input file name: ')
pdb = list(open(filename))


sheets = list()
locs = dict()
atoms = dict()
sheet_stuff = list()
sys.stdout=open('final.pdb','w')
for line in pdb:
	data = line.split()
	#print(data)
	if not data:
		continue
	if data[0] == "SHEET":
		#print("Beta sheet")
		sheet_stuff.append(line)
		start = data[6]
		end = data[9]
		sheets.append(BetaSheet(start, end))
	elif data[0] == "ATOM" and data[2] == "CA":
		locs[str(data[5])] = data[6:9]
		atoms[str(data[5])] = Atom(line)
for sheet in sheets:
	sheet.load(atoms)
	align_to = sheets[0]
for sheet in sheets[1:]:
	sheet.move_to(align_to)
for line in sheet_stuff:
	#print(line[:-1])
	for atom in atoms:
		print(atoms[atom])

def transform_vector(start1, end1, start2, end2):
	"""accepts 4 points in 3-space specifying the beginning and end of two vectors. returns 2 points in 3-space that represent the second vector translated/rotated to match the first"""
	v1 = tuple(map(lambda *x: x[1] - x[0], list(zip(start1, end1)))) #first vector in standard form
	v2 = tuple(map(lambda *x: x[1] - x[0], list(zip(start2, end2)))) #second vector in standard form
	l1 = sum(tuple(map(lambda c: c**2, v1)))**0.5 #magnitude of first vector
	l2 = sum(tuple(map(lambda c: c**2, v2)))**0.5 #magnitude of second vector
	beginning = start1
	end = tuple(map(lambda c: c * l2 / l1, v1))
	return (beginning, end)
#f.close()