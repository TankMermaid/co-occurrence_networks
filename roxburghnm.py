

import numpy as np
import csv
import matplotlib.pyplot as plt
import sys

print "The species name is: ", sys.argv[1]
filename = str(sys.argv[1])
ifile = open(filename, 'rb')
reader = csv.reader(ifile, delimiter='\t')

mat1=[]
for row in reader:
	mat1.append(row)
	
mat = np.zeros((len(mat1),len(mat1[0])))

for i in range(len(mat1)):
	for j in range(len(mat1[0])):
		mat[i,j] = mat1[i][j]

ifile.close()

import os
for file in os.listdir("."):
    if file.endswith(".dat"):
        print file


# Edge contacts
def edges(mat,i,j,L):
	res=0	
	if i ==0 and (j>0 and j<L-1):
		tot = mat[i+1,j] + mat[i,j-1] + mat[i,j+1]
		res = tot/3
	elif i==L-1 and (j>0 and j<L-1):
		tot = mat[i-1,j] + mat[i,j-1] + mat[i,j+1]
		res = tot/3
	elif (i>0 and i<L-1) and j==0:
		tot = mat[i-1,j] + mat[i,j+1] + mat[i+1,j]
		res = tot/3
	elif (i>0 and i<L-1) and j==L-1:
		tot = mat[i-1,j] + mat[i,j-1] + mat[i+1,j]
		res = tot/3
	elif (i>0 and i<L-1) and (j>0 and j<L-1):
		tot = mat[i-1,j] + mat[i+1,j] + mat[i,j-1] + mat[i,j+1]
		res = tot/4
	return res
	
# Corner contacts
def corners(mat,i,j,L):
	res=0
	if i ==0 and (j>0 and j<L-1):
		tot = mat[i+1,j-1] + mat[i+1,j+1]
		res = tot/2
	elif i==L-1 and (j>0 and j<L-1):
		tot = mat[i-1,j-1] + mat[i-1,j+1]
		res = tot/2
	elif (i>0 and i<L-1) and j==0:
		tot = mat[i-1,j+1] + mat[i+1,j+1]
		res = tot/2
	elif (i>0 and i<L-1) and j==L-1:
		tot = mat[i-1,j-1] + mat[i+1,j-1]
		res = tot/2
	elif (i>0 and i<L-1) and (j>0 and j<L-1):
		tot = mat[i-1,j-1] + mat[i+1,j+1] + mat[i+1,j-1] + mat[i-1,j+1]
		res = tot/4
	return res
	
# Open areas
def opens(mat,i,j,L):
	res = 0
	tot = 1
	if (i>0 and i<L-1) and (j>0 and j<L-1):
		tot = mat[i-1,j-1] + mat[i+1,j+1] + mat[i+1,j-1] + mat[i-1,j+1]
	if tot==0:
		res=1
	return res

# Solid areas
def solids(mat,i,j,L):
	res=0
	tot=1
	if (i>0 and i<L-1) and (j>0 and j<L-1):
		tot = mat[i-1,j-1] + mat[i+1,j+1] + mat[i+1,j-1] + mat[i-1,j+1]
	if tot==4:
		res=1
	return res
sol=0
e=0
for i in range(51):
	for j in range(51):
		if mat[i,j]==1:
			sol+= corners(mat,i,j,51)
			e+=1

def corners_contacts(mat):
	dim = np.shape(mat);sol=0
	for i in range(dim[0]):
		for j in range(dim[0]):
			if mat[i,j]==1:
				sol+= corners(mat,i,j,dim[0])
	return sol
	
def edges_contacts(mat):
	dim = np.shape(mat);sol=0
	for i in range(dim[0]):
		for j in range(dim[0]):
			if mat[i,j]==1:
				sol+= edges(mat,i,j,dim[0])
	return sol
			
def solids_areas(mat):
	dim = np.shape(mat);sol=0
	for i in range(dim[0]):
		for j in range(dim[0]):
			if mat[i,j]==1:
				sol+= solids(mat,i,j,dim[0])
	return sol			
			
def opens_areas(mat):
	dim = np.shape(mat);sol=0
	for i in range(dim[0]):
		for j in range(dim[0]):
			if mat[i,j]==1:
				sol+= opens(mat,i,j,dim[0])
	return sol		



def phi(matnew, edges_obs, corners_obs, solids_obs, opens_obs):
	edges_r = edges_contacts(matnew)
	corners_r = corners_contacts(matnew)
	solids_r = solids_areas(matnew)
	opens_r = opens_areas(matnew)
	if solids_obs==0:
		phir = abs( round((edges_r/edges_obs),2) - 1) + abs(round((corners_r/corners_obs),2) -1) + abs((opens_r/opens_obs) -1)
	else:
		phir = abs( round((edges_r/edges_obs),2) - 1) + abs(round((corners_r/corners_obs),2) -1) + abs((solids_r/solids_obs)-1) + abs((opens_r/opens_obs) -1)
	#print edges_r,corners_r,solids_r,opens_r
	#print edges_obs,corners_obs,solids_obs,opens_obs
	#print phir
	return phir


# Eliminate last row and column

mat1 = np.delete(mat,50,0)
mat = np.delete(mat1,50,1)


#------------------------- Calculate observed edge, corner, solid and open areas
e_obs = edges_contacts(mat)
c_obs = corners_contacts(mat)
s_obs = solids_areas(mat)
o_obs = opens_areas(mat)
			
# Splitting matrix into 10 x 10 blocks
division = 5
temp=[]
for k in range(5):
	for l in range(5):
		print np.arange(k*10,k*10+10),np.arange(l*10,l*10+10)
		temp.append(mat[k*10:k*10+10,l*10:l*10+10])





	#1) Rotation/reflection
rotemp=[]
for i in range(len(temp)):
	rotemp.append(np.roll(np.rot90(temp[i]),5,axis=1))
#2) shuffle blocks
sel = np.arange(len(temp))
np.random.shuffle(sel)
shuffled_rotemp =[]
for i in range(len(rotemp)):
	shuffled_rotemp.append(rotemp[sel[i]])


simulations=1000
sims=[]
phis =[]
for sim in range(simulations):
	print "Simulation: %d" % sim
	if sim > 0:
		rotemp=[]
		for i in range(len(temp)):
			rotemp.append(np.roll(np.rot90(temp[i]),5,axis=1))
		sel = np.arange(len(temp))
		np.random.shuffle(sel)
		shuffled_rotemp =[]
		for i in range(len(rotemp)):
			shuffled_rotemp.append(rotemp[sel[i]])
	n=0
	alpha=0.01
	while n < 1:
		sel = np.arange(len(shuffled_rotemp))
		np.random.shuffle(sel)
		shuffled_rotemp2 =[]
		for i in range(len(shuffled_rotemp)):
			shuffled_rotemp2.append(np.roll(np.rot90(shuffled_rotemp[sel[i]]),np.random.randint(10),axis=0))
		shuffled_rotemp = shuffled_rotemp2
		for s in range(700):
			if s % 1000==0 and s>0:
				a = np.random.randint(len(sel))
				b = np.random.randint(len(sel))
				#print "Block %d and block %d" % (a,b)
				while a != b:
					b = np.random.randint(len(sel))
				tempshuf = shuffled_rotemp[a]
				shuffled_rotemp[a] = np.rot90(np.roll(shuffled_rotemp[b],np.random.randint(10),axis=1))
				shuffled_rotemp[b] = np.rot90(np.roll(tempshuf,np.random.randint(10),axis=1))
			if s == 0:
				temp1 = temp
				matnew = np.vstack( (np.hstack((shuffled_rotemp[0:5])), np.hstack((shuffled_rotemp[5:10])), np.hstack((shuffled_rotemp[10:15])), np.hstack((shuffled_rotemp[15:20])), np.hstack((shuffled_rotemp[20:25])) ) )
				phiold = round(phi(matnew,e_obs,c_obs, s_obs, o_obs),3)
			prob_out = np.random.lognormal(0.1,1)
			while prob_out > 1:
				prob_out = np.random.lognormal(0.1,1)
			if np.random.random() > prob_out:
				counts=0
				b1 = np.random.randint(len(shuffled_rotemp))
				b2 = np.random.randint(len(shuffled_rotemp))
				pres = np.where(shuffled_rotemp[b1]==1)
				abse = np.where(shuffled_rotemp[b2]==0)
				while pres[0].size == 0:
					b1 = np.random.randint(len(shuffled_rotemp))
					pres = np.where(shuffled_rotemp[b1]==1)
				while abse[0].size == 0:
					b2 = np.random.randint(len(shuffled_rotemp))
					abse = np.where(shuffled_rotemp[b2]==0)
				while counts < 10:
					k = np.random.randint(len(pres[0]))
					l = np.random.randint(len(abse[0]))
					oldcoordx = pres[0][k]
					oldcoordy = pres[1][k]
					newcoordx = abse[0][l]
					newcoordy = abse[1][l]
					shuffled_rotemp[b1][oldcoordx,oldcoordy] = 0
					shuffled_rotemp[b2][newcoordx,newcoordy] = 1
					matnew = np.vstack( (np.hstack((shuffled_rotemp[0:5])), np.hstack((shuffled_rotemp[5:10])), np.hstack((shuffled_rotemp[10:15])), np.hstack((shuffled_rotemp[15:20])), np.hstack((shuffled_rotemp[20:25])) ) )
					phinew = round(phi(matnew,e_obs,c_obs, s_obs, o_obs),3)
					if phinew < phiold:
						#print phinew, phiold
						phiold = phinew
						break
					else:
						#print phinew, phiold
						shuffled_rotemp[b1][oldcoordx,oldcoordy] = 1
						shuffled_rotemp[b2][newcoordx,newcoordy] = 0
					counts +=1
					if counts == 3:
						continue
			else:
				counts=0
				i = np.random.randint(len(shuffled_rotemp))
				pres = np.where(shuffled_rotemp[i]==1)
				while pres[0].size == 0:
					i = np.random.randint(len(shuffled_rotemp))
					pres = np.where(shuffled_rotemp[i]==1)
				abse = np.where(shuffled_rotemp[i]==0)
				#print abse
				while counts < 10:
					kk = np.random.randint(len(pres[0]))
					ll = np.random.randint(len(abse[0]))
					oldcoordx = pres[0][kk]
					oldcoordy = pres[1][kk]
					newcoordx = abse[0][ll]
					newcoordy = abse[1][ll]
					shuffled_rotemp[i][oldcoordx,oldcoordy] = 0
					shuffled_rotemp[i][newcoordx,newcoordy] = 1
					matnew = np.vstack( (np.hstack((shuffled_rotemp[0:5])), np.hstack((shuffled_rotemp[5:10])), np.hstack((shuffled_rotemp[10:15])), np.hstack((shuffled_rotemp[15:20])), np.hstack((shuffled_rotemp[20:25])) ) )
					phinew = round(phi(matnew,e_obs,c_obs, s_obs, o_obs),3)
					if phinew < phiold:
						#print phinew, phiold
						phiold = phinew
						break
					else:
						#print phinew, phiold
						shuffled_rotemp[i][oldcoordx,oldcoordy] = 1
						shuffled_rotemp[i][newcoordx,newcoordy] = 0
					counts +=1
					if counts == 10:
						continue
			if phinew <= alpha:
				break
		if phinew <= alpha:
			break
		n+=1
	sims.append(matnew)
	phis.append(phiold)


name = filename.split('.')
outfile = open(name[0]+"_matrices.csv","wb")
wr = csv.writer(outfile)
for f in range(len(sims)):
	wr.writerows(sims[f])
	wr.writerow("\n")
	
outfile.close()


