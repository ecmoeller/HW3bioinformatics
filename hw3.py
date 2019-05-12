import math
from operator import itemgetter 

import plotly
plotly.tools.set_credentials_file(username='moell162', api_key='XWeDQVjHzfvEN01WsE7p')

import plotly.plotly as py
import plotly.graph_objs as go
import plotly.tools as tls
import matplotlib.pyplot as plt

# Create random data with numpy
import numpy as np

class Pair:
    def __init__(self, proteinA, proteinB):
        self.proteinA = proteinA
        self.proteinB = proteinB  


def q2(pairs):
    # Measure the degree of each protein with at least 1 interaction in the network (exclude selfinteractions)
    degs = {} #protein name to deg num
    neighbors = {} #protein name to list of protein names who are neighbors

    for p in pairs:

        #exclude self interactions
        if(p.proteinA != p.proteinB):
            #Increment proteinA count in degs if found, otherwise add it
            if (p.proteinA not in degs):
                degs[p.proteinA] = 1
                tempList = [p.proteinB]
                neighbors[p.proteinA] = tempList
            else:
                degs[p.proteinA] = degs[p.proteinA] + 1
                #add proteinB to the list that's already there
                tempList = neighbors[p.proteinA]
                tempList.append(p.proteinB)
                neighbors[p.proteinA] = tempList

            #Increment proteinB count in degs if found, otherwise add it
            if (p.proteinB not in degs):
                degs[p.proteinB] = 1
                tempList = [p.proteinA]
                neighbors[p.proteinB] = tempList
            else:
                degs[p.proteinB] = degs[p.proteinB] + 1  
                #add proteinA to the list that's already there
                tempList = neighbors[p.proteinB]
                tempList.append(p.proteinA)
                neighbors[p.proteinB] = tempList      


    # Plot the degree distribution of the protein-protein interaction network (a histogram is fine)
    fig = plt.figure()
    
    n = list(degs.values())
    bins = 75
    plt.hist(n, bins)
    plt.title("Degree distribution of protein-protein interaction network")
    plt.xlabel("Degree Count")
    plt.ylabel("Frequency")

    plotly_fig = tls.mpl_to_plotly( fig )
    py.plot(plotly_fig, filename='all_probe_data_before_log')


    #finding highest degree protein
    max = 0
    maxP = ""
    for k, v in degs.items():
        if v > max:
            max = v
            maxP = k

    print("Max degree value ", max)
    print("Corresponding max protein ", maxP)

    return degs, neighbors

# degs is a dictionary mapping from the protein name to the degree of interaction of that protein
def q3b(degs, neighbors): 

    #Compute clustering coefficient for every protein in the graph
    print("calculating clustering coefficient")

    #Calculate CC for each node
        #Kv: its degree -> stored in degs
        #Nv: number of links between neighbors of V
        #CC(V) = 2 Nv / Kv (Kv-1)

    coeffs = {} #dict that maps protein name to CC
    #Loop through all nodes
    degCount = 0
    doNothing = True
    for key, value in degs.items():
        #Calculate CC for each node
        kv = value
        
        #Need to caluculate Nv, which is links between neighbors
        nList = neighbors[key] #nList is [protein1, protein2, protein3]
        linkedNeighCount = 0
        linkFound = [] #contains tuples of links (proteinA, proteinB)
        for i, n1 in enumerate(nList):

            j = i
            while (j < len(nList)):
            #for n2 in nList:

                n2 = nList[j]
                #Nodes are the same
                if(n1 == n2):
                    #skip
                    doNothing = True
                    
                #Neighbor nodes don't match
                else:
                    n1Neighbors = neighbors[n1]
                    n2Neighbors = neighbors[n2]
                    #Are n1 and n2 in each others neighbors list? If so, we found a link
                    if(n1 in n2Neighbors and n2 in n1Neighbors):

                        #Check if you already have counted this link
                        found = False
                        for link in linkFound:
                            if((link[0] == n1 and link[1] == n2) or (link[0] == n2 and link[1] == n1)):
                                #Do not add
                                found = True
                                break
                            
                        #Did not find link in list so it's a new link!
                        if(found == False):
                            linkFound.append((n1,n2))
                            linkedNeighCount += 1

                    #Not in each other's list so move on to next 
                    else:
                        doNothing = True
                j += 1
        
        #CC(V) = 2 Nv / Kv (Kv-1)
        if(kv != 0 and kv != 1):
            cc = (2 * linkedNeighCount) / (kv * (kv - 1))
        else:
            cc = 0
        coeffs[key] = cc
        degCount += 1

    #Print all coefficients 
    for key, value in coeffs.items():
        print("Key and val", key, value)

    return coeffs
            
def q3graph(degs, coeffs, name):

    degsList = list(degs.values())
    coeffsList = list(coeffs.values())


    # Create a trace
    trace = go.Scatter(
        x = degsList,
        y = coeffsList,
        mode = 'markers'
    )

    layout = go.Layout(
        title=go.layout.Title(
            text=name
            
        ),
        xaxis=go.layout.XAxis(
            title=go.layout.xaxis.Title(
                text='Degree of Interaction'
            )
        ),
        yaxis=go.layout.YAxis(
            title=go.layout.yaxis.Title(
                text='Clustering Coefficient'
            )
        )
    )

    data = [trace]

    fig = go.Figure(data=data, layout=layout)
    py.plot(fig, filename=name)
   
def q4a(litDegs, degs):
    # What is the Pearson correlation between interaction degrees in the systematically mapped 
    # network and in the literature curated network?  (use the proteins in common between the 
    # two networks after you exclude self-interactions in the systematically mapped network) 
    # Discuss your interpretation of this

    print("Pearson correlation time")

    ldegrees = []
    degrees = []

    for key, value in litDegs.items():
        #make sure lit degree value is an int

        #need to first find the proteins in common and add them in order
        protein = degs.get(key)
        #In common
        if(protein != None):
            ldegrees.append(int(value))
            degrees.append(degs[key])

    print("Length of ldegrees", len(ldegrees))
    print("Length of degrees", len(degrees))

    pear = pearson(ldegrees, degrees) 
    return pear


def pearson(litDegs, degs):
    #Computing Pearson similarity coefficient
    # The vectors are in the correct format

    input1 = litDegs
    input2 = degs

    samplesize = len(input2)

    xsum = 0
    ysum = 0
    xysum = 0
    xsquaresum = 0
    ysquaresum = 0
    for i in range(len(input1)):
        x = float(input1[i])
        y = float(input2[i])
        xy = x * y
        xsquare = x * x
        ysquare = y * y
        xsum += x
        ysum += y
        xysum += xy
        xsquaresum += xsquare
        ysquaresum += ysquare

    numerator = (samplesize * xysum) - (xsum * ysum)
    denominator = math.sqrt((samplesize * xsquaresum - xsum * xsum) * (samplesize * ysquaresum - ysum * ysum))
    return (numerator / denominator)

def q4b(degs, coeffs, litDegs):
    #Find an example of a protein with more than 10 interactions in the Rolland et al. network, a clustering coefficient
    #of greater than 0.2, and no interactions in the literature curated network
    proteins = []

    for key, value in degs.items():
        #more than 10 interactions in Rolland et al. network
        if(value > 10):
            cc = coeffs[key]
            #a clustering coefficient of greater than 0.2
            if(cc > 0.2):
                ld = litDegs[key]
                #no interactions in the literature curated network
                if(ld == 0):
                    proteins.append(key)


def main():
    print("In main")
    #Reading in data from txt file
    with open('Human_PPI.csv','r') as f:
        data = f.readlines()

    with open("Lit_degrees.csv", "r") as f2:
        data2 = f2.readlines()

    count = 0
    pairs = []
    for line in data:
        line = line.strip()
        if (count == 0):
            print("did nothing")
        if(count >= 1 ):
            pairterms = line.split('\t')
            pair = Pair(pairterms[0], pairterms[1])
            pairs.append(pair)

        count += 1

    count = 0
    degrees = {}
    for line in data2:
        line = line.strip()
        if (count == 0):
            print("did nothing")
        if(count >= 1 ):
            degreeLine = line.split('\t')
            degrees[degreeLine[0]] = int(degreeLine[1])

        count += 1

    #print("First line ", degrees["A1CF"])

    degs, neighbors = q2(pairs)
    print("Length of degs", len(degs))

    coeffs = q3b(degs, neighbors)

    #Graph clustering coefficients and degrees
    q3graph(degs, coeffs, "Degree of Interaction and Clustering Coefficient for all Proteins")

    #Pearson Coefficient 
    pearsonCo = q4a(degrees, degs)
    print("This is the pearson coefficient ", pearsonCo)

    q4b(degs, coeffs, degrees)

    #print("Neighbors of KRTAP9-2")
    #print(neighbors["KRTAP9-2"])

    #Question q3c
    # PSMC3. How many interaction partners does it have and what are the functions of these proteins? 
    # What is its clustering coefficient? 
    # Are these consistent with what is known about the function of PSMC3? 
    print("Number of interaction partners", degs["PSMC3"])
    print("The interaction partners", neighbors["PSMC3"])
    print("Clustering Coefficient", coeffs["PSMC3"])


if __name__ == '__main__':
    main()
