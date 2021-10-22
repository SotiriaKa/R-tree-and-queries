#Sotiria Kastana, 2995
#py assignment2.py data_rectangles.txt

import sys
import csv
import math
from operator import itemgetter

# MEROS I___________________________________________________________________________________________________________________________________________________________________________________________
block_size = 1024														#block size
record_size = 36 														#record size
rects_sorted = list()
rtree = open('rtree.txt','w')

with open(sys.argv[1], 'r') as f:
	rects = list(f)
	rects.sort(key = lambda x: float(x.split('\t')[1]))				    #taksinomisi me vasi to x-low
	for r in rects:
		rects_sorted.append(r)
f.close()
rects_objs = len(rects_sorted)
capacity = block_size/record_size										#plithos eggrafwn pou xwraei kathe komvos																		
leaves = int(rects_objs/capacity) + 1                       			#fylla -> anwfli
current_line = 0
sqr_l = int(math.sqrt(leaves)) + 1										#riza fyllwn -> anwfli
low_pos = sqr_l

while low_pos < rects_objs:												#ana riza(fylla) taksinomisi me vasi to y-low
	sort_y = list(rects_sorted[current_line:low_pos])
	sort_y.sort(key = lambda y: float(y.split('\t')[3]))
	rects_sorted[current_line:low_pos] = sort_y
	current_line = low_pos
	low_pos += sqr_l

if current_line < (rects_objs - 1) :									#gia osa perissevoun
	sort_y = list(rects_sorted[current_line:rects_objs])
	sort_y.sort(key = lambda y: float(y.split('\t')[3]))
	rects_sorted[current_line:rects_objs] = sort_y
#_____________________________________________________________________________________________________________________________________________________________________________________
#layer 1 ->  dhmiourgia fyllwn
layer = 1   
areas = [] 																#mia lista pou se kathe thesi tis apothikevei to emvado kathe eggrafis
numofnodes_meanarea = [] 												#mia lista pou krataei se kathe thesi tis mia lista me 3 theseis (dld 1h thesi gia layer, 2h gia to plithos komvwn se kathe epipedo kai 3h to meso emvado se kathe epipedo)
nodes = []  															#afto ine to array pou prostheto tous komvous
i = 0   																#metritis gia fylla 
j = 0	  																#metritis twn eggrafwn se kathe fyllo
def area(xmin,xmax,ymin,ymax):											#synartisi pou upologizei emvado enos orthogoniou
	x = xmax-xmin
	y = ymax - ymin
	ar = x*y
	return ar

def leaf_record(ln, ARRA, areas):										#synartisi pou epistrefei ton pinaka ARRA prostithontas se afton kainouria eggrafi kai ton pinaka areas prostithontas se afton to emvado tis neas eggrafis
	arra = [0]*5 
	arra[0] = int(rects_sorted[ln].split('\t')[0])
	arra[1] = float(rects_sorted[ln].split('\t')[1])
	arra[2] = float(rects_sorted[ln].split('\t')[2])
	arra[3] = float(rects_sorted[ln].split('\t')[3])
	arra[4] = float(rects_sorted[ln].split('\t')[4])
	ARRA[a]=arra    
	areas.append(area(arra[1],arra[2],arra[3],arra[4])) 						
	return ARRA, areas

#ya ola ta leaves-1	
for i in range(leaves-1): 												#ektos apo to teleftaio fyllo, gt mporei na min gemisei olo
	node_id = i
	ARRA = [0]*32														#enas pinakas 32 thesewn gia kathe komvo
	ARRA[0] = node_id													#stin 1h thesi tou apothikevete to node_id tou komvou(iso me ti thesi tou ston pinaka nodes)
	ARRA[1] = capacity													#stin 2h thesi tou apothikevete to plithos twn eggrafwn tou komvou
	ARRA[2] = 0 														# =0 -> dld den einai o telefteos komvos toy epipedou
	ARRA[3] = layer														#stin 4h thesi apothikevete to epipedo sto opoio anikei o komvos
	a = 4 																#apo tin 5h thesi kai meta se kathe thesi tou yparxei mia eggrafi pou deixnei se kathe antikeimeno
	for j in range(i*capacity,capacity*(i+1)):  						#ya kathe fyllo/komvo, ana capacity eggrafes kalei tin synartisi 
		ARRA, areas = leaf_record(j, ARRA, areas)
		a += 1
	nodes.append(ARRA)													#afou gemisei o komvos prostithetai ston pinaka nodes

#ya to telftaio fyllo poy exei <= capacity 
node_id += 1
ARRA = [0]*32															
ARRA[0] = node_id
ARRA[1] = rects_objs - (leaves-1)*capacity
ARRA[2] = 1																# = 1 -> dld einai o teleftaios komvos toy epipedou
ARRA[3] = layer
a = 4        															
for k in range(((leaves-1)*capacity),rects_objs):
	ARRA, areas = leaf_record(k, ARRA, areas)
	a += 1
nodes.append(ARRA)
sum_areas = 0

for ar in areas:														#ypologismos mesou emvadou olwn twn orthogoniwn kathe epipedou
	sum_areas += ar
mean = sum_areas/len(areas)
numofnodes_meanarea.append([layer,leaves,mean])
layer += 1
areas = []  															 #midenizete ya kathe neo epipedo
#_____________________________________________________________________________________________________________________________________________________________________________________
#layer += 1 ->  dhmiourgia ypoloipwn komvwn
root = 0 
node_id += 1															# an to root einai 0 simenei oti dn exei ftasei akoma stin riza

def node_record(n, a, ARRA, areas):										#synartisi antoistixh me thn leaf_record
	arra = [0]*5 
	m1= 4 + nodes[n][1]													#to m1 vriskete, gt mporei an anaferete ston pio deksi komvo na min einai gemismenes kai oi 32 theseis
	m2 = nodes[n][1] - 1 												
	mlx = sorted(nodes[n][4:m1], key=itemgetter(1))						#taksinomisi me vasi to x-low kai pairnw tin min timi (dld 1i thesi)
	mhx = sorted(nodes[n][4:m1], key=itemgetter(2)) 					#taksinomisi me vasi to x-high kai pairnw tin max timi (dld teleftaia thesi)
	mly = sorted(nodes[n][4:m1], key=itemgetter(3))						#antistoixa gia y
	mhy = sorted(nodes[n][4:m1], key=itemgetter(4))
	max_low_X = mlx[0][1]
	max_high_X = mhx[m2][2]
	max_low_Y = mly[0][3]
	max_high_Y =mhy[m2][4]
	arra[0] = nodes[n][0]
	arra[1] = max_low_X
	arra[2] = max_high_X
	arra[3] = max_low_Y
	arra[4] = max_high_Y
	ARRA[a] = arra
	areas.append(area(arra[1],arra[2],arra[3],arra[4])) 
	return ARRA, areas
 
while (root == 0):														#ya oso den exo dhmiourgisei tin riza
	previous_layer = 0
	n = 0 																#deixnei posous komvous exo prosperasei 
	while previous_layer == 0:											#ya oso vriskome se diaforetiko epipedo apo to proigoumeno moy epipedo (dld mexri na ftaso sto epipedo twn paidiwn mou)
		if nodes[n][3] == layer-1:										#otan vro to prwto komvo pou einai sto proigoumeno epipedo apo mena, dld to prwto paidi mou
			previous_layer = 1
			n -= 1 														# giati tha to ayksisei stin epomeni grammi tou kodika mou, ara gia na meino ousiastika se afti tin thesi pou vrika oti ksekinaei to previous layer moy
		n +=1
	my_layer_first_node = node_id										#kratao to id apo ton 1o komvo tou proigoumenou mou layer gia na to xrisimopoihso parakatw kai na vlepo posoi komvoi mou menoun apo to proigomeno epipedo kai na ksero an eine o teleftaios komvos
	i = 0
	while my_layer_first_node - n > capacity:   						#dld gia kathe komvo tou pinaka nodes, diavazo tous komvous mexris otou na ftaso se kapoion opoy exei < capaicty, ara simenei oti tous ypoloipous mporo na tous prostheso mazi ston telfteo komvo
		root += 1
		ARRA = [0]*32													
		ARRA[0] = node_id
		ARRA[1] = capacity
		ARRA[2] = 0									
		ARRA[3] = layer
		a = 4 
		for l in range(0,capacity): 			
			ARRA, areas = node_record(n, a, ARRA, areas)
			a += 1
			n += 1
		nodes.append(ARRA)
		node_id += 1
		i += 1
	ARRA = [0]*32														#ya ton deksiotero/teleftaio komvo kathe epipedoy									
	ARRA[0] = node_id
	ARRA[1] = my_layer_first_node - n
	ARRA[2] = 1									
	ARRA[3] = layer
	a = 4 
	for k in range(n,my_layer_first_node): 
		ARRA, areas = node_record(n, a, ARRA, areas)
		a += 1
		n += 1
	nodes.append(ARRA)
	sum_areas = 0
	for ar in areas:													
		sum_areas += ar
	mean = sum_areas/len(areas)
	number_of_nodes = node_id + 1 - my_layer_first_node
	numofnodes_meanarea.append([layer,number_of_nodes,mean])
	node_id += 1
	layer += 1 
	areas = []
	if root == 0 :														#an to root meinei iso me 0 simenei oti dn mpike stin while, ara o protos komvos pou dimiourgise einai kai teleftaios tou epipedou, dld exei eggrafes < capacty -> RIZA 
		root = 1
		layer -= 1 														#afairo ena apo to layer giati ine o telfteos komvos(riza) ara dn prepei na prostheso layer ya to epomeno
		node_id -= 1 													#to idio kai ya to node_id
		print("h arithmisi twn nodes ksekinaei apo 0!!")
		print("height of R-tree is: "+str(layer)+" layers")
		for nm in range(0,len(numofnodes_meanarea)):
			print("# nodes in layer: "+str(numofnodes_meanarea[nm][0])+" are: "+str(numofnodes_meanarea[nm][1])+" with mean_area: "+str(numofnodes_meanarea[nm][2]))
	else :
		root = 0														#allios an einai != 0 simenei oti mpike stin proti while, ara oxi riza, ara to ksanakanoume 0 gia na ksanampei stin while
rtree.write(str(node_id)+"\n")
rtree.write(str(layer)+"\n")
nn = 0 
for nn in range(nn,len(nodes)):
	rtree.write(str(nodes[nn][0])+", "+str(nodes[nn][1])+", "+str(nodes[nn][4:32])+"\n")

# MEROS II_____________________________________________________________________________________________________________________________________________________________________________________________
#ksero oti to layer twn fyllwn = 1

with open('query_rectangles.txt', 'r') as q:
	q_rects = list(q)													#h lista me ta query rectangles
q.close()
	
def intersect(q, nod):
	qxl = float(q.split('\t')[1])
	qxh = float(q.split('\t')[2])
	qyl = float(q.split('\t')[3])
	qyh = float(q.split('\t')[4])
	nxl = float(nod[1])
	nxh = float(nod[2])
	nyl = float(nod[3])
	nyh = float(nod[4])
	if qxl > nxh or qxh < nxl:
		return False
	if qyl > nyh or qyh < nyl:		
		return False
	return True

def inside(q, nod):
	qxl = float(q.split('\t')[1])
	qxh = float(q.split('\t')[2])
	qyl = float(q.split('\t')[3])
	qyh = float(q.split('\t')[4])
	nxl = float(nod[1])
	nxh = float(nod[2])
	nyl = float(nod[3])
	nyh = float(nod[4])
	if qxl <= nxl :
		if qxh >= nxh:
			if qyl <= nyl:
				if qyh >= nyh:
					return True
	return False
	
def contains(q, nod):
	qxl = float(q.split('\t')[1])
	qxh = float(q.split('\t')[2])
	qyl = float(q.split('\t')[3])
	qyh = float(q.split('\t')[4])
	nxl = float(nod[1])
	nxh = float(nod[2])
	nyl = float(nod[3])
	nyh = float(nod[4])
	if qxl >= nxl:
		if qxh <= nxh:
			if qyl >= nyl:
				if qyh <= nyh:
					return True
	return False

	
def Range_Intersect_Query(q, node, rects_intersects, N):
	if node[3]== 1: 																					#if node == leaf
		for n in range(4,node[1]+4):																	#for each index entry e in n 
			if intersect(q, node[n]):																	#such that inersect
				rects_intersects+=1
	else:
		for n in range(4,node[1]+4):
			if intersect(q, node[n]):
				N+=1
				rects_intersects, N = Range_Intersect_Query(q, nodes[node[n][0]], rects_intersects, N) 	#dld ksanakalo me ton komvo pou exei node id tou to ptr tis eggrafis
	return rects_intersects, N
	
def Range_Inside_Query(q, node, rects_inside, N):
	if node[3]== 1: 																					#if node == leaf
		for n in range(4,node[1]+4):																	#for each index entry e in n 
			if inside(q, node[n]):																		#such that inside
				rects_inside+=1
	else:
		for n in range(4,node[1]+4):
			if inside(q, node[n]):
				N+=1
				rects_inside, N = Range_Inside_Query(q, nodes[node[n][0]], rects_inside, N) 	
			elif intersect(q, node[n]):
				N+=1
				rects_inside, N = Range_Inside_Query(q, nodes[node[n][0]], rects_inside, N) 	
	return rects_inside, N

def Containment_Query(q, node, rects_contain, N):
	if node[3]== 1: 																					#if node == leaf
		for n in range(4,node[1]+4):																	#for each index entry e in n 
			if contains(q, node[n]):																	#such that contains
				rects_contain+=1
	else:
		for n in range(4,node[1]+4):
			if contains(q, node[n]):
				N+=1
				rects_contain, N = Containment_Query(q, nodes[node[n][0]], rects_contain, N) 			
	return rects_contain, N

for qr in q_rects:
	res_intersection = 2*[0]	
	res_inside = 2*[0]	
	res_containment = 2*[0]	
	N = 1											   					 							   #metritis ya tous komvous pou prospeulanw -> gt ksekinaw sigoura apo 1 komvo!!
	rects_intersects = 0								   											   #metrtis ya ta orthogonia pou periexonte sto q_rect
	res_intersection = Range_Intersect_Query(qr, nodes[len(nodes)-1], rects_intersects, N)
	N = 1											   					
	rects_inside = 0
	res_inside = Range_Inside_Query(qr, nodes[len(nodes)-1], rects_inside, N)
	N = 1											   					
	rects_contain = 0
	res_containment = Containment_Query(qr, nodes[len(nodes)-1], rects_contain, N)
	print("query: "+str(qr.split('\t')[0]))
	print("intersects with "+str(res_intersection[0])+" rectangles  -> "+str(res_intersection[1])+" nodes")
	print("inside to "+str(res_inside[0])+" rectangles -> "+str(res_inside[1])+" nodes")
	print("contains to "+str(res_containment[0])+" rectangles -> "+str(res_containment[1])+" nodes")
	print("___________________________________________________________")