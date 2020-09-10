
import numpy as np
import random
from matplotlib import pyplot as plt
import copy as cp
from scipy import stats

########################################## FUNCTIONS ###########################################
# N number of filaments
# l length of filaments
# L length of domain

def dynamics_all(N_0,l,L,timestep,finaltime,delta,D,regions,start,end,K_e,K_n,K_d,p_e,type_pe):

	### length of ecad region
	dist_ecad = 0
	for j in range(len(end)):
		dist_ecad += end[j]-start[j]

	### for every filament (1,N) keep track of:
	###  x - centres of mass, p - orientation
	### stress_individual - total stress on individual filament, y - for plotting purposes
	N = N_0
	x = [[random.uniform(0,L) for i in range(N)]]
	p = [[random.choice([-1,1]) for i in range(N)]]
	y = [[random.uniform(0,2) for i in range(N)]]
	stress_individual = [[0 for i in range(N)]]
	
	### initialise time and its indexing
	time = 0 
	k = 0

	### discretise the domain -L,2L to calculate the total stress in space 
	# regions of length L included so that we can include measurements across the periodic BC's
	xpoints = 300 # << make sure this is a multiple of three
	disc_x = np.linspace(-L, 2*L, num=xpoints)
	stress = np.zeros((round(finaltime/timestep)+1,xpoints))
	# this is the resultant stress in (0,L); we flip the ghost regions to calculate the stress in the region
	stress_final = np.zeros((round(finaltime/timestep)+1,int(xpoints/3)))
	# fraction of actin filaments bound
	fraction_bound = np.zeros(round(finaltime/timestep)+1)
	

	while time<finaltime:

		x_i = x[k]
		p_i = p[k]
		
		stress_individual_i = np.zeros((N))
		### determines movement for each filament
		Del = np.zeros((N))

		# Check whether filaments are bound or unbound
		s = [0 for i in range(N)]
		# only initialise for constant stress
		for i in range(0,N):
			if type_pe == 0:
				if ecad(x_i[i],regions,start,end) == 1 and random.random()<p_e:
					s[i] = 1
		# note: if type_pe = 1; 
		# initially no filaments are bound > in the next time step their state will depend on individual stress  

		for i in range(0,N):
			for j in range(i+1,N):
				
				##### check for filaments close enough, but further than delta and one or both unbound!!!
				if dist(x_i[i],x_i[j],L) < l and dist(x_i[i],x_i[j],L) > delta and (s[i]+s[j])<2:
					
					# find whethere xi or xj is smaller >>> PERIODIC BCs <<<<<
					if abs(x_i[i]-x_i[j])<l:
						

						if p_i[i]*p_i[j]==-1:
							if p_i[i]==-1:
								del_i = delta
								del_j = -delta
							else:
								del_i = -delta
								del_j = delta
							xm = random.uniform(max(x_i[i],x_i[j])-l/2,min(x_i[i],x_i[j])+l/2)

						else:
							if x_i[i] < x_i[j]:
								del_i = delta
								del_j = -delta
								if p_i[i]==1:
									xm = x_i[i] + l/2
								else:
									xm = x_i[j] - l/2
							else:
								del_i = -delta
								del_j = delta
								if p_i[j]==1:
									xm = x_i[j] + l/2
								else:
									xm = x_i[i] - l/2

					elif abs(x_i[i]+L-x_i[j])<l:

						x_i[i] += L

						if p_i[i]*p_i[j]==-1:
							if p_i[i]==-1:
								del_i = delta
								del_j = -delta
							else:
								del_i = -delta
								del_j = delta
							xm = random.uniform(max(x_i[i],x_i[j])-l/2,min(x_i[i],x_i[j])+l/2)

						else:
							del_i = -delta
							del_j = delta
							if p_i[j]==1:
								xm = min(x_i[i],x_i[j]) + l/2
							else:
								xm = max(x_i[i],x_i[j]) - l/2

						
					elif abs(x_i[i]-L-x_i[j])<l:

						x_i[i] -=  L

						if p_i[i]*p_i[j]==-1:
							if p_i[i]==-1:
								del_i = delta
								del_j = -delta
							else:
								del_i = -delta
								del_j = delta
							xm = random.uniform(max(x_i[i],x_i[j])-l/2,min(x_i[i],x_i[j])+l/2)

						else:
							del_i = delta
							del_j = -delta
							if p_i[j]==1:
								xm = min(x_i[i],x_i[j]) + l/2
							else:
								xm = max(x_i[i],x_i[j]) - l/2

					########### if the filaments are unbound then filament moves and stress is calculated
					if s[i]==0 and del_i is not 0:
						Del[i] += del_i
						[stress[k,:],stress_individual_i[i]] =  stress_fun(x_i[i],del_i,l,xm,disc_x,stress[k,:],stress_individual_i[i])
					if s[j]==0 and del_j is not 0:
						Del[j] += del_j
						[stress[k,:],stress_individual_i[j]] =	stress_fun(x_i[j],del_j,l,xm,disc_x,stress[k,:],stress_individual_i[j])

		# ## plots for checking whether stresses are calculated correctly
		# plt.clf()
		# fig,axs = plt.subplots(2)
		# for j in range(0,N):
		# 	if p_i[j] == -1:
		# 		axs[0].plot([x_i[j]-l/2,x_i[j]+l/2],[(j)*0.2+0.2,(j)*0.2+0.2],'b-')
		# 		plt.ylim(0,1)
				
		# 	else:
		# 		axs[0].plot([x_i[j]-l/2,x_i[j]+l/2],[(j)*0.2+0.2,(j)*0.2+0.2],'r-')

		# axs[0].plot(disc_x,stress[k,:])
		# axs[0].plot(xm,0.3,'g.')
		# if k > 0:
		# 	axs[1].plot(disc_x,np.mean(stress[0:k,:],axis=0))
		# axs[0].set_ylim(-1,1)
		# axs[1].set_ylim(-1,1)
		# plt.savefig('trial/trial_n_'+str(k)+'.png')
				
		### get rid of ghost areas in stress calculation
		stress_final[k,:] = np.sum(np.transpose(np.reshape(stress[k,:], (3,int(xpoints/3) ))),axis=1)
		fraction_bound[k] = np.sum(s)/N

		### move all the centres of mass 
		x.append(list(Del + np.array(x[k])))
		stress_individual.append(list(stress_individual_i))

		p.append(cp.deepcopy(p[k]))
		y.append(cp.deepcopy(y[k]))

		### diffusion
		for j in range(0,N):
			x[-1][j] = x[-1][j] + random.gauss(0,D*timestep)
			if x[-1][j] < 0:
				x[-1][j] += L
			if x[-1][j] > L:
				x[-1][j] -= L


		### check how many filaments are nucleated in normal & ecad region and how many are removed
		N_n = stats.poisson.rvs(K_n)
		N_e = stats.poisson.rvs(K_e)
		N_d = stats.poisson.rvs(K_d*N)
		### change to total number of filaments
		N += N_n + N_e - N_d

		for i in range(N_n):
			r = random.uniform(0,L-dist_ecad)
	
			for j in range(len(start)):
				if r > start[j]:
					r += end[j]-start[j]
			x[-1].extend([r,])
			p[-1].extend([random.choice([-1,1]),])
			y[-1].extend([random.uniform(0,10),])
			stress_individual[-1].extend([0,])

		for i in range(N_e):
			r = random.uniform(0,dist_ecad) + start[0]

			for j in range(len(start)):
				if r > end[j]:
					r += -end[j]+start[j+1]
			x[-1].extend([r,])
			p[-1].extend([random.choice([-1,1]),])
			y[-1].extend([random.uniform(0,10),])
			stress_individual[-1].extend([0,])

		for i in range(N_d):
			kill = random.randint(0,N)
			del x[-1][kill]
			del y[-1][kill]
			del p[-1][kill]
			del stress_individual[-1][kill]


		### check whether for the next time step filaments will be bound or unbound
		for i in range(0,N):
			# constant stress
			if type_pe == 0:
				if ecad(x[-1][i],regions,start,end) == 1 and random.random()<p_e:
					s[i] = 1
			else:
				if ecad(x[-1][i],regions,start,end) == 1:
					stress_middle = 50;
					sd = 20;
					prob = p_e*np.exp(-(abs(stress_individual[-1][i])-stress_middle)**2/(2*sd**2))
					if random.random()<prob:
						s[i] = 1

		### update timestep and time index
		time += timestep
		k += 1


	return x,p,y,disc_x,stress_final,fraction_bound
	# note: x - centres of mass, p - orientations (+1,-1) and y - rand number between (0,10) for plotting purposes!
	# (disc_x,stress_final) - plotting stress, fraction_bound - franction of bound filaments 


# periodic distance between two filaments
def dist(x1,x2,L):
	return min(abs(x1-x2),abs(x1-x2+L),abs(x1-x2-L))

# establish ecad region
def ecad(x,regions,start,end):
	for i in range(0,regions):
		if x > start[i] and x < end[i]:
			return 1
	return 0

def stress_fun(x,Del,l,xm,dx,stress,stress_individual):
	f = np.zeros(len(dx))
	if Del>0:
		for kk in range(len(dx)):
			if x-l/2 <= dx[kk] and xm >= dx[kk]:
				f[kk] += 1/l*(dx[kk]-x+l/2)
			elif xm < dx[kk] and dx[kk]<=x+l/2:
				f[kk] += 1/l*(dx[kk]-x+l/2) - 1 
	else:
		for kk in range(len(dx)):
			if x-l/2 <= dx[kk] and xm >= dx[kk]:
				f[kk] += -1/l*(dx[kk]-x+l/2)
			elif xm < dx[kk] and dx[kk]<=x+l/2:
				f[kk] += -1/l*(dx[kk]-x+l/2) +1

	stress_individual += np.sum(f)
	stress += f

	return [stress,stress_individual]

########################################## RUN PROGRAM ###########################################
N_0 = 500 # original number of filaments
l = 1 # actin length
L = 100 # region length
delta = 0.1 # move during one timestep (velocity of motor*timestep)
D = 0 # diffusion coeff
regions = 1    # ecad localisation and number of islands
start = [L/2 - 10,] # regions have to be specified in ascending order!
end = [L/2 + 10,]
timestep = 0.1
finaltime = 31
K_e = 0 #1*timestep # rate of filament nucleation in ecad region 
K_n = 0 #100*timestep # rate of filament nucleation in normal region 
K_d = 0 #0.001*timestep # rate of filament decay
pe_list = np.array([0.25,0.5,0.75,1.]) # probability of filament being bound in ecad region
type_pe = 1


for ii in range(len(pe_list)):
	p_e = pe_list[ii]
	x,p,y,disc_x,stress,fraction_bound = dynamics_all(N_0,l,L,timestep,finaltime,delta,D,regions,start,end,K_e,K_n,K_d,p_e,type_pe)
	np.savetxt('/Users/aplochocka/Desktop/numerics_1d_tissue/bound_'+'kn_'+str(K_n)+'_pe2_'+str(int(p_e*100))+'.txt',fraction_bound)
	
	max_stress = max(max(k) for k in stress)
	min_stress = min(min(k) for k in stress)

	for i in range(0,len(x)-1,100):

		y_p = []
		y_n = []


		plt.clf()
		plt.figure(1)
		for j in range(len(x[i])):
			if p[i][j] == -1:
				plt.plot([x[i][j]-l/2,x[i][j]+l/2],[y[i][j],y[i][j]],'b-')
				y_n.extend([x[i][j],])

			else:
				plt.plot([x[i][j]-l/2,x[i][j]+l/2],[y[i][j],y[i][j]],'r-')
				y_p.extend([x[i][j],])

		plt.xlim(0-l/2,L+l/2)
		plt.ylim(0,2)
		plt.savefig('act/actin_ke_'+str(K_e)+'_kn_'+str(K_n)+'_pe2_'+str(int(p_e*100))+'n_'+str(i)+'.png')

		plt.clf()
		fig,axs = plt.subplots(2)
		counts, bins = np.histogram(y_p,bins = np.arange(0,L,2))
		axs[0].hist(bins[:-1],bins,weights=counts/max(counts),color = 'red')
		axs[0].plot([start,start],[0,1],'k-')
		axs[0].plot([end,end],[0,1],'k-')
		axs[0].set_xlim(0,L)
		axs[0].set_ylim(0,1.1)
		counts2, bins2 = np.histogram(y_n,bins = np.arange(0,L,2))
		axs[1].hist(bins2[:-1],bins2,weights=counts2/max(counts2),color = 'blue') 
		axs[1].plot([start,start],[0,1],'k-')
		axs[1].plot([end,end],[0,1],'k-')
		axs[1].set_xlim(0,L)
		axs[1].set_ylim(0,1.1)
		plt.savefig('dist/dist_ke_'+str(K_e)+'_kn_'+str(K_n)+'_pe2_'+str(int(p_e*100))+'n_'+str(i)+'.png')




		plt.clf()
		fig,axs = plt.subplots(2)
		xpoints = 300
		disc_x = np.linspace(0, L, num=xpoints/3)
		axs[0].plot(disc_x,stress[i,:])
		axs[0].plot([start,start],[0,1],'k-')
		axs[0].plot([end,end],[0,1],'k-')
		axs[0].set_xlim(0,L)
		axs[0].set_ylim(min_stress,max_stress)
		if i > 0:
			axs[1].plot(disc_x,np.mean(stress[0:i,:],axis=0))
			axs[1].plot([start,start],[0,1],'k-')
			axs[1].plot([end,end],[0,1],'k-')
			axs[1].set_xlim(0,L)
			axs[1].set_ylim(min_stress,max_stress)
		plt.savefig('stress/actin_ke_'+str(K_e)+'_kn_'+str(K_n)+'_pe2_'+str(int(p_e*100))+'n_'+str(i)+'.png')



