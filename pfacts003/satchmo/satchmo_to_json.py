import json, string
from Bio.SubsMat import MatrixInfo
from Bio.Nexus import Trees, Nodes

blosum62 = MatrixInfo.blosum62

def makeJsonTree(phy_output):

    #Retrieve tree from tree file and make a JSON object
    _writer = json.JsonWriter()

    #Making JSON for the SATCHMO tree in .phy file
    treeFile = phy_output
    f = open(treeFile)
    tree = ''.join([line.rstrip() for line in f.readlines()])
    f.close()

    ret = {}
    tree = tree.replace(',',':1.0,')
    tree = tree.replace(')',':1.0)')
    tree = tree.replace(';',':1.0;')
    ret['__tree__'] = tree
    ret['__seed_protein_sequence_id__'] = 1002
    
    ans = _writer.write(ret)
    return ans


def readTree(phy_output):
    #Retrieve tree from tree file and make a JSON object
    _writer = json.JsonWriter()

    #Making JSON for the SATCHMO tree in .phy file
    treeFile = phy_output
    f = open(treeFile)
    tree = ''.join([line.rstrip() for line in f.readlines()])
    f.close()
    tree = tree.replace(',',':1.0,')
    tree = tree.replace(')',':1.0)')
    #Changed this - removed ;
    tree = tree.replace(';',':1.0')
    return tree

def mptt(tree, leftId, nodeId, nodeIdOfLeftId):
    rightId = leftId + 1
    for s in tree.node(nodeId).get_succ():
        rightId = mptt(tree, rightId, s, nodeIdOfLeftId)
    nodeIdOfLeftId[leftId] = nodeId
    return (rightId + 1)
    
    
def makeJsonAlign(smo_output, phy_output):
    _writer = json.JsonWriter()
    #Making JSON for root alignment in .smo file
    treeStr = readTree(phy_output)
    #Find left and right IDs using MPTT
    tree = Trees.Tree(rooted = True)
    (clades, vals) = tree._parse(treeStr)
    tree._add_subtree(0, clades)
    rootNode = tree.node(tree.root)
    nodeIdOfLeftId = {}
    leftIdOfSatchmoId = {}
    mptt(tree, 1, 0, nodeIdOfLeftId)
    #tree.display()
    #print(nodeIdOfLeftId)

    smoFile = smo_output
    #Read from .smo file
    f_smo = open(smoFile)
    smoString = ''.join([line for line in f_smo.readlines()])
    f_smo.close()
    #Find the number of nodes
    index1 = string.find(smoString, ' ')
    index2 = string.find(smoString, '\n')
    no_of_nodes = int(smoString[index1 + 1:index2])
    #Read sequence alignment at each node
    curr_index = 0
    i = 0
    ret = {}
    keys = {}
    no_of_leaves = {}
    while i < no_of_nodes:
        curr_index = string.find(smoString, '[', curr_index)
        leaf_index = string.find(smoString, "leaf", curr_index)
        internal_index = string.find(smoString, "internal", curr_index)

        if((leaf_index < internal_index) or (internal_index == -1)) and (leaf_index != -1):
            nodeType = "leaf"
            #print("Its a leaf!")
        elif((leaf_index > internal_index) or (leaf_index == -1))and (internal_index != -1):
            nodeType = "internal"
            #print("Its an internal node!")
        else:
            print("Invalid node type!")

        if(nodeType == "leaf"):
            curr_index = string.find(smoString, "index", curr_index)
            newLine_index = string.find(smoString, '\n', curr_index)
            satchmo_id = int(smoString[curr_index + len("index") + 1:newLine_index])
            curr_index = string.find(smoString, '>', curr_index)
            newLine_index = string.find(smoString, '\n', curr_index)
            seq_name = smoString[curr_index+1:newLine_index]
            for id in nodeIdOfLeftId.keys():
                if(tree.is_terminal(nodeIdOfLeftId[id])) and (seq_name in tree.get_taxa(nodeIdOfLeftId[id])):
                    leftIdOfSatchmoId[satchmo_id] = id
                    populate_ret(ret, id, curr_index, keys, no_of_leaves, smoString)

        elif(nodeType == "internal"):
            curr_index = string.find(smoString, "index", curr_index)
            newLine_index = string.find(smoString, '\n', curr_index)
            satchmo_id = int(smoString[curr_index + len("index") + 1:newLine_index])
            target_index = string.find(smoString, "target", curr_index)
            newLine_index = string.find(smoString, '\n', target_index)
            target = int(smoString[target_index + len("target") + 1:newLine_index])
            target = nodeIdOfLeftId[leftIdOfSatchmoId[target]]
            #print("Target = ." + str(target) + ".");
            template_index = string.find(smoString, "template", curr_index)
            newLine_index = string.find(smoString, '\n', template_index)
            template = int(smoString[template_index + len("template") + 1:newLine_index])
            template = nodeIdOfLeftId[leftIdOfSatchmoId[template]]
            #print("Template = ." + str(template) + ".")
            #print(tree.node(nodeIdOfLeftId[id]).get_succ())
            for id in nodeIdOfLeftId.keys():
                #print(tree.node(nodeIdOfLeftId[id]).get_succ())
                if(not(tree.is_terminal(nodeIdOfLeftId[id]))) \
                and (template in tree.node(nodeIdOfLeftId[id]).get_succ()) \
                and (target in tree.node(nodeIdOfLeftId[id]).get_succ()) \
                and (len(tree.node(nodeIdOfLeftId[id]).get_succ()) == 2):
                    #print("Template = ." + str(template) + ".")
                    #print(tree.node(nodeIdOfLeftId[id]).get_succ())
                    #print("Adding to id = " + str(id))
                    leftIdOfSatchmoId[satchmo_id] = id
                    curr_index = string.find(smoString, '>', curr_index)
                    populate_ret(ret, id, curr_index, keys, no_of_leaves, smoString)
                
        i += 1
        #print(leftIdOfSatchmoId)
        #print(ret)
    for id in nodeIdOfLeftId.keys():
        addColorString(ret, keys, no_of_leaves, id)
        #print(ret)
    ans = _writer.write(ret)
    return ans
    #Find where the root node alignment starts
#    index3 = string.find(smoString, "index " + str(no_of_nodes - 1), 0)
#    index3 = string.find(smoString, "{", index3)
    #Find name and alignment (at root) of all sequences
#    index4 = string.find(smoString, '}', index3)
  
#    index5 = string.find(smoString, '>', index3)

#    ret = {}
def populate_ret(ret, id, curr_index, keys, no_of_leaves, smoString):
    #print("Populating ret...")
    ret[id] = {}
    keys[id] = {}
    no_of_leaves[id] = 0
    index5 = curr_index
    while index5 != -1:
        #print("String near index5 = " + smoString[index5:index5+15])
        index6 = string.find(smoString, '\n', index5)
        seq_name = smoString[index5+1:index6]
        last = string.find(smoString, '//', index6)
        index5 = string.find(smoString, '>', index6, last)
        if index5 == -1:
            index5 = string.find(smoString, '//', index6)
        seq_alignment = smoString[index6 + 1:index5]
        seq_alignment = string.replace(seq_alignment, '\n', '')
        index5 = string.find(smoString, '>', index6, last)
        ret[id][seq_name] = {}
        ret[id][seq_name]['aligned_sequence'] = seq_alignment
        #print(seq_alignment)
        keys[id][no_of_leaves[id]] = seq_name
        no_of_leaves[id] += 1
        #print("Index5 = " + str(index5))
        #print("Seq name = " + keys[no_of_leaves - 1])
    #print("Leaves = " + str(no_of_leaves))
    #addColorString(ret, keys, no_of_leaves)
    #ans = _writer.write(ret)
    #return ans

#@param
#ret: A collection of sequences
#keys: An array that stores names of all sequences
#strIndex: Current index/column under consideration in the sequences
#Output: Column Score calculated from Blosum62 matrix and the character
#with maximum frequency
def findAvgColScore(ret, keys, no_of_leaves, strIndex, maxCharList, id):
    start = 0
    maxScore = 0
    retObject = []

    for char in maxCharList:
        charCount = 0
        blosumScoreForCol = 0
        i = 0
        #Calculate the sum of blosum62 of each residue againt the
        #most highly conserved residue
        while i < no_of_leaves[id]:
            charCount += 1
            if(str.islower(ret[id][keys[id][i]]['aligned_sequence'][strIndex])) or (ret[id][keys[id][i]]['aligned_sequence'][strIndex] == "."):
                return 0
            blosumScoreForCol += findBlosumScore(ret[id][keys[id][i]]['aligned_sequence'][strIndex],char)

            if((ret[id][keys[id][i]]['aligned_sequence'][strIndex] == "-")):
                charCount -= 1
            i += 1
        blosumScoreForCol -= findBlosumScore(char, char)

        charCount -= 1
        #Find the character with the best blosum62 score
        if(blosumScoreForCol > maxScore) or (start == 0):
            maxScore = blosumScoreForCol
            maxChar = str(char)
            start = 1

    if (charCount <= 0):
        avgColScore = 0
    else:
        avgColScore = (maxScore / charCount)

    retObject.append(avgColScore)
    retObject.append(maxChar)
    return retObject


#@param
#avgColScore: The average column score calculated from Blosum62 matrix
#Output: The color code in accordance with Belvu coloring
def belvuColor(avgColScore):
    if avgColScore >= 3:
        color = "C"
    elif avgColScore >= 1.5:
        color = "B"
    elif avgColScore >= 0.5:
        color = "G"
    else:
        color = "N"
    return color

#@param
#ret: A collection of sequences
#keys: An array that stores names of all sequences
#strIndex: Current index/column under consideration in the sequences
#Output: A list of all residues that are most conserved in a column
def findMostConsResidue(ret, keys, no_of_leaves, strIndex, id):
    charCount = {}
    i = 0
    while i < no_of_leaves[id]:
        charCount[ret[id][keys[id][i]]['aligned_sequence'][strIndex]] = 0
        i += 1

    #Calculte frequency of each residue
    i = 0
    while i < no_of_leaves[id]:
        if(ret[id][keys[id][i]]['aligned_sequence'][strIndex] != "-"):
            charCount[ret[id][keys[id][i]]['aligned_sequence'][strIndex]] += 1
        i += 1
            
    maxCount = 0
    maxCharList = []
    i = 0
    #Make a list of most highly occuring residues
    while i < no_of_leaves[id]:
        if (charCount[ret[id][keys[id][i]]['aligned_sequence'][strIndex]] > maxCount):
            maxCount = charCount[ret[id][keys[id][i]]['aligned_sequence'][strIndex]]
            maxCharList = []
            maxCharList.append(ret[id][keys[id][i]]['aligned_sequence'][strIndex])
        elif (charCount[ret[id][keys[id][i]]['aligned_sequence'][strIndex]] == maxCount) and (maxCharList.count(ret[id][keys[id][i]]['aligned_sequence'][strIndex]) == 0):
            maxCharList.append(ret[id][keys[id][i]]['aligned_sequence'][strIndex])        
        i += 1
    return maxCharList


#@param
#index1 and 2: Indices of blosum62 matrix
#Output: Corresponding value in bloum62 matrix
def findBlosumScore(index1, index2):
     try:
         score = float(blosum62[(index1, index2)])
     except KeyError, e:
         try:
             score = float(blosum62[(index2, index1)])
         except KeyError, e:
             score = 0
             if((index1 != "-") and (index2 != "-")): 
                 print("Invalid indices for Blosum62 matrix!");
     return score

             
def addColorString(ret, keys, no_of_leaves, id):

    i = 0
    while i < no_of_leaves[id]:
        ret[id][keys[id][i]]['colors'] = ""
        i += 1
            
    j = 0
    #print("Leaves = " + str(no_of_leaves));
    while j < len(ret[id][keys[id][0]]['aligned_sequence']):
        smallFlag = 0
        #Columns with small case letters will not be colored
        if(str.islower(ret[id][keys[id][0]]['aligned_sequence'][j])) or (ret[id][keys[id][0]]['aligned_sequence'][j] == "."):
            i = 0
            while i < no_of_leaves[id]:
                ret[id][keys[id][i]]['colors'] += "N"
                i += 1
            smallFlag = 1
        if(smallFlag == 0):
            #Find most highly conserved residue type
            maxCharList = findMostConsResidue(ret, keys, no_of_leaves, j, id)
            #print(maxCharList)
            #Find blosum62 score for each column
            listObject = findAvgColScore(ret, keys, no_of_leaves, j, maxCharList, id)
            avgColScore = listObject[0]
            maxChar = listObject[1]
            #print("Avg Score = " + str(avgColScore) + " and maxChar = " + maxChar)
            #Decide the color of each column using Belvu style coloring
            color = belvuColor(avgColScore)
            #print("Color = " + color)
            #Make color string for each sequence
            i = 0
            while i < no_of_leaves[id]:
#                if(ret[keys[i]]['aligned_sequence'][j] == "I"):
#                    print("Color now = " + color)
                if(ret[id][keys[id][i]]['aligned_sequence'][j] == "-"):
                    ret[id][keys[id][i]]['colors'] += "N"
                elif (findBlosumScore(ret[id][keys[id][i]]['aligned_sequence'][j], maxChar) >= avgColScore):
                    ret[id][keys[id][i]]['colors'] += color
                else:
                    ret[id][keys[id][i]]['colors'] += "N"
                i += 1
        j += 1
    #print("colors : " + ret[keys[0]]['colors'])
