# -*- coding: utf-8 -*-
"""
A library for parsing and handling of information form nucmer's delta outfile. For use in RepeatM contig matching.

Created on Thu May 30 13:42:15 2019

@author: Daniel Rawlinson, Australian Centre for Ecogenomics (ACE)
@email: daniel.rawlinson@uqconnect.edu.au
"""
'''
TODOS - 
    -function for detecting and dealing with overlap situations
    -partial match situaitons
    -Don't just pass significant Nucmer_Matches into a list. Make a subclass with inherited functions, add threshold level attribute. 

'''


'''
imports
    '''
import os




'''------------------------begin class definition---------------------------'''

#TODO - match.draw() to actually illustrate the alignment?
class Nucmer_Match(object):
    '''a format in which to store each separate nucmer sequence alignment'''
    seqs = None
    lengths = None
    hitstarts_1 = None
    hitstops_1 = None
    hitstarts_2 = None
    hitstops_2 = None
    simerrors = None
    
    
    ''' self is the nucmer_match object, match is the seq names and seq lengths,
    match deets is a dictionary containing separate arrays of hit starts (x2), hit stops (x2),
    and similarityerrors'''
    def __init__(self, match, matchdeets):
        
        #TODO - need float here or integer? integer only works for python 3
        self.seqs =        [match[0],match[1]]
        self.lengths =     list(map(int, [match[2], match[3]]))
        self.hitstarts_1 = list(map(int, [entry[0] for entry in matchdeets]))
        self.hitstops_1 =  list(map(int, [entry[1] for entry in matchdeets]))
        self.hitstarts_2 = list(map(int, [entry[2] for entry in matchdeets]))
        self.hitstops_2 =  list(map(int, [entry[3] for entry in matchdeets]))
        self.simerrors =   list(map(int, [entry[5] for entry in matchdeets]))
        
    def __len__(self):
        '''what happens when you call length of the class object'''
        return len(self.simerrors)
        
    def display(self):
        print('> '+(' '.join(self.seqs))+'\n')
        print('\t\tSeq1 region\tSeq2 region\tSimilarityErrors')
        for i in range(len(self)):
            print('Match {}\t\t{}:{}\t\t{}:{}\t\t{}'.format(i+1, self.hitstarts_1[i], self.hitstops_1[i], self.hitstarts_2[i], self.hitstops_2[i], self.simerrors[i]))

    def get_ani(self):
        identity_list = []
        for i in range(len(self)):
            alignlen = float(self.hitstops_1[i]) - self.hitstarts_1[i]
            identity_list.append((alignlen - self.simerrors[i]) / alignlen)
        
        ANI = sum(identity_list) / len(self)           
            
        return ANI
    
    
    #TODO - fix abs to put in correct place
    def get_align1(self):
        matchlength_1 = sum([float(self.hitstops_1[i]) - self.hitstarts_1[i] for i in range(len(self))])
        alignment_ratio_1 = abs(matchlength_1/self.lengths[0])
        
        return alignment_ratio_1
    
    def get_align2(self):
    
        matchlength_2 = sum([float(self.hitstops_2[i]) - self.hitstarts_2[i] for i in range(len(self))])
        alignment_ratio_2 = abs(matchlength_2/self.lengths[1])
        
        return alignment_ratio_2
    
    def get_lengthratio(self):
        length_ratio = float(self.lengths[0]) / self.lengths[1]
        
        return length_ratio
        
        

    def gen_statistics(self):
        '''produce length ratio and bi-direcitonal alignment statistics to measure closeness of sequences''' 
        #Fractions should all be >1 now. Any deviations from this may constitute interesting alignment relationships
        #One variable in each converted to float so division works properly in Python 2

        length_ratio = float(min(self.lengths)) / max(self.lengths)
        
        alignment_ratio_1 = self.get_align1()
        
        alignment_ratio_2 = self.get_align2()

        ANI = self.get_ani()         
           
        
        return length_ratio, alignment_ratio_1, alignment_ratio_2, ANI
        pass
    
    def apply_threshold(self, threshold = 0.97):
        passthresh = False
        stats = self.gen_statistics()
        if all(stat >= threshold for stat in stats):
            passthresh = True
            
        return passthresh
#TODO - promote to sig_match class
#    def promote(self, threshold):
#        return Sig_Match(self, threshold)
#
##subclass (significant-matches) when passed threshold tests
#class Sig_Match(Nucmer_Match):
#    super().__init__():
#    pass
        
'''--------------------------end class definition---------------------------'''


#'''Temp file for testing'''
file = 'C:\\Users\\dan_r\\Documents\\Honours_Data\\nucmertest\\testsplitvserrors\\splitvserrors_0.2_errors.delta'

''' Init arrays for alignment data of each match.
    Placed together, starts, stops, and errors will form a matrix with number of alignments equal to width
    '''

def deltaread(file): #if reading in a whole bunch of .delta files, record these in a list to preserve function of later funcitons
    '''parse delta files'''
    with open(file, 'r') as f:
        self_match_count = 0
        inputline = f.readline().strip().split(" ")
        deltaname = '---'.join([os.path.basename(elem) for elem in inputline])
        recording = False #"recording" (boolean variable) is a switcher that will determine whether the line is being stored or not (ie. to ignore positions of indels)
        matchdeets = []
        seen_matches = []
        name = False #do we have a match ready to write yet?
        delta = [] #name of dictionary for storage of Nucmer_Match objects
        for line in f.readlines(): #parsing of line with match info
            if line.startswith('>'):
                if name: #skip flushing to dictionary if this is the first match record
                    #FIRST - flush previous hit to a nucmer_match object
                    if [match[1], match[0]] not in seen_matches: #if reverse has not already been read
                        delta.append(Nucmer_Match(match, matchdeets))  
                        seen_matches.append([match[0], match[1]])
                    #remove alignment details of previous match
                    del matchdeets[:]
                #THEN - read in match details
                match = line.replace('>', '').split()
                if match[0] == match[1]: #skip if it's a match to itself
                    self_match_count +=1 #for use in assert. But that's harder to do than originally planned
                    continue
                else:
                    recording = True
                    name = True
                
                continue
            elif line == '0\n':
                recording = True
                continue
            elif recording == True:
                matchdeets.append(line.split()[:-1])
                
                recording = False
    
    #flush last match to  nucmer_match object
    delta.append(Nucmer_Match(match, matchdeets))

    deltadict = {}
    deltadict[deltaname] = delta
    return deltadict

def dict_threshold(deltadict, threshold = 0.90, collate = False, outfile = None):
    '''input = dictionary output from deltaparse function
    apply blanket ratio threshold level and select to process further'''
    thresh_matches = []
    for value in deltadict.values():
        for match in value:
            if match.apply_threshold(threshold) == True:
                thresh_matches.append(match)
    match_dict = {}            
    match_dict[next(iter(deltadict.keys())) +"---"+str(threshold)] = thresh_matches #this won't work if I've read in heaps of delta files into the one dictionary
    
    #optionally output only significant Nucmer_Match objects, for collation into list containing all matches from multiple deltafiles
    if collate == True:
        return thresh_matches
    #optionally write to file   
    elif outfile:
        write_thresh_matches(match_dict, outfile)
    #just return dictionary for interactive use
    else:
        return match_dict

    

def write_thresh_matches(match_dict, filename):
    #TODO - write stats to file
    '''write simple txt file containing sequence matches that pass threshold'''
    header = next(iter(match_dict.keys())).split("---")
    matches = next(iter(match_dict.values()))
    with open(filename, 'w') as f:
        f.write(header[0]+" "+header[1]+" at threshold of "+ header[2]+"\n")
        for match in matches:
            f.write(match.seqs[0] + ' '+ match.seqs[1] + '\n')


#TODO - function to export in JSON format?