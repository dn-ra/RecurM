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
    -Change alignment % from total length to an interval data type
'''


'''
imports
    '''
import os
import intervals



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
        #Earlier method
        #matchlength_1 = sum([float(self.hitstops_1[i]) - self.hitstarts_1[i] for i in range(len(self))])
        #alignment_ratio_1 = abs(matchlength_1/self.lengths[0])
        
        for i in range(len(self)):
            span = intervals.closed(getattr(self, 'hitstarts_{}'.format(1))[i], getattr(self, 'hitstops_{}'.format(1))[i]) #minus one for 0-based indexing
            if i ==0:
                c_span = span
            else:
                c_span = c_span | span
        total_length = sum([ivl.upper - ivl.lower for ivl in c_span]) #must convert list output of upper and lower to int
        alignment_ratio_1 = total_length / self.lengths[0]
        
        
        return alignment_ratio_1
    
    def get_align2(self):
    
        #Earlier method
        #matchlength_2 = sum([float(self.hitstops_2[i]) - self.hitstarts_2[i] for i in range(len(self))])
        #alignment_ratio_2 = abs(matchlength_2/self.lengths[1])
        
        for i in range(len(self)):
            #need to handle reverse complement sequence coordinates
            coords = [self.hitstarts_2[i], self.hitstops_2[i]]
            
            span = intervals.closed(min(coords), max(coords))
            if i ==0:
                c_span = span
            else:
                c_span = c_span | span
        total_length = sum([ivl.upper - ivl.lower for ivl in c_span])
        alignment_ratio_2 = total_length / self.lengths[1]

        return alignment_ratio_2
    
    def get_lengthratio(self):
        length_ratio = float((min(self.lengths)) / max(self.lengths))
        
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
    
    def apply_threshold(self, threshold = 0.97, stats = False): #option to precalculate stats 
        passthresh = False
        if stats == False:
            stats = self.gen_statistics()
            
        if all(stat >= threshold for stat in stats):
            passthresh = True
            
        return passthresh
    
    
    def is_fragment(self, upperthreshold = 0.9, lowerthreshold = 0.9, stats = False): #option to precalculate stats
        if stats == False:
            stats = self.gen_statistics()
            
        min_align = min(stats[1:3])
        max_align = max(stats[1:3])
        ani = stats[3]  
        if max_align > upperthreshold and min_align < lowerthreshold and ani > upperthreshold:
            return True
    
        else:
            return False
    
    def label(self):
        #TODO - tighten definitions
        
        '''
        options = 
            -linear/circular
            -simple/complex (or 3 regions)
            -perfect/imperfect
            -bring a fourth or somehow indicate when aligned terminally (perfect) but there are gaps in the middle
            
        also: flags for outcomes that haven't been programmed for yet
            '''
        
        '''match labels'''
        orientation = None
        complexity = None
        alignment = None
        
        '''arrangement labels'''
        arrange_a = None
        arrange_b = None
        
        
        '''iterate through alignments'''
        for align in range(len(self)): #for each aligned region
            seq_a = [self.hitstarts_1[align], self.hitstops_1[align]]
            seq_b = [self.hitstarts_2[align], self.hitstops_2[align]]
            
            orig_seq_b = seq_b #there has to be a better way than this. keeping orig_seq_b for testing linear perfect alignment
            if seq_b[0] > seq_b[1]: #if reverse complement
                seq_b = [a*-1 for a in seq_b] #switch to negative. This way logic can be treated the same
                b_direction = 'reverse'
            else:
                b_direction = 'forward'
                
            #first alignment block 
            if align == 0:
                orig_b_direction = b_direction
                initial_coords = [seq_a, seq_b]
                if seq_a == [1, self.lengths[0]] and [min(orig_seq_b), max(orig_seq_b)] == [1, self.lengths[1]]: #using orig_seq_b to avoid problems with negative values
                    alignment = 'perfect'
                else:
                    alignment = 'imperfect'
                orientation = 'linear'
                complexity = 'simple'
                
            #second alignment block    
            elif align ==1:
                #if 2 blocks are in different directions
                if orig_b_direction != b_direction:
                    orientation = 'different directions'
                    #flag and break for further programming
                    break
                
                
                #place alignment perfect or imperfect at the end? how do I treat it if the alignment goes all the way to the ends but there is some sequence in the middle that isn't aligned anywhere?
                
                '''seq_a arrangement'''
                
                if min(seq_a) > max(initial_coords[0]) or max(seq_a) < min(initial_coords[0]): #if second alignment is outside the range of the first
                    arrange_a = 'outside'
                elif min(seq_a) > min(initial_coords[0]) and max(seq_a) < max(initial_coords[0]): #if second alignment is completely within the first
                    arrange_a = 'within'
                elif seq_a[0] > initial_coords[0][0] and seq_a[0] < initial_coords[0][1] and seq_a[1] > initial_coords[0][1]: ##second alignment straddles boundaries
                    arrange_a = 'straddle'
                else:
                    arrange_a = '' #if not caught
                

                '''seq_b arrangement'''
                if min(seq_b) > max(initial_coords[1]): # if second alignment of b is also out of range of 1st
                    arrange_b = 'outside right'
                elif  max(seq_b) < min(initial_coords[1]):
                    arrange_b = 'outside left'
                elif max(seq_b) < max(initial_coords[1]) and min(seq_b) > min(initial_coords[1]): #if second alignment is completely within the first
                    arrange_b = 'within'
                elif seq_b[0] > initial_coords[1][0] and seq_b[0] < initial_coords[1][1] and seq_b[1] > max(initial_coords[1]):#straddles the boundary
                    arrange_b = 'straddle_right'
                elif seq_b[1] > initial_coords[1][0] and seq_b[1] < initial_coords[1][1] and seq_b[0] < min(initial_coords[1]):#straddles the boundary
                    arrange_b = 'straddle_left'
                else:
                    arrange_b = '' #if not caught
                
                
                
                '''interpret arrangements'''
                
                if arrange_a == 'within':
                    if arrange_b == 'within':
                        #orientation and alignment should already be set
                        complexity = 'complex'
                    elif arrange_b == 'outside left' or arrange_b == 'outside right':
                        complexity = 'complex'
                        pass
                    elif arrange_b == 'straddle_left' or arrange_b == 'straddle_right':
                        complexity = 'complex'
                    else:
                        complexity = 'not caught'
                
                elif arrange_a =='outside':
                    if arrange_b == 'within':
                        complexity = 'complex'
                    elif arrange_b == 'outside right':
                        complexity = 'complex'
                        orientation = 'linear'
                    elif arrange_b == 'outside left':
                        orientation = 'circular'
                        complexity = 'complex'
                        #this is the complex circular arrangement. no overlap but clearly in opposite directions. Is it actually circular?
                    elif arrange_b == 'straddle_left':
                        complexity = 'complex'
                        orientation = 'circular'
                    elif arrange_b == 'straddle_right':
                        complexity = 'complex'
                        orientation= 'linear'
                    else:
                        complexity = 'not caught'
                
                elif arrange_a =='straddle':
                    if arrange_b == 'within':
                        complexity = 'complex'
                        orientation = 'linear'
                    elif arrange_b == 'outside right':
                        complexity = 'complex'
                        orientation = 'linear' 
                    elif arrange_b == 'outside left':
                        complexity = 'complex'
                        orientation = 'circular'
                        #is this a circle that doesn't fully enclose?
                    elif arrange_b == 'straddle_left':
                        orientation = 'circular'
                        complexity = 'simple'
                    elif arrange_b == 'straddle_right':
                        orientation = 'linear'
                        complexity = 'complex'
                    else:
                        complexity = 'not caught'
                else:
                    complexity = 'not caught'
                    
            #third alignment blocks 
            elif align > 1:
                #always set as complex when there are more than 2 aligned regions?
                complexity = '3 regions'
                #flag and break for further programming
                break
            
        #check if alignment as all the way to the end
        if [min(self.hitstarts_1), max(self.hitstops_1)] == [1, self.lengths[0]] and [min(self.hitstarts_2 + self.hitstops_2), max(self.hitstarts_2 + self.hitstops_2)] == [1,self.lengths[1]]:
            alignment = 'perfect'
        else:
            alignment = 'imperfect'
            
        return [orientation, complexity, alignment]
            
            
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
        match = None
        delta = [] #name of dictionary for storage of Nucmer_Match objects
        for line in f.readlines(): #parsing of line with match info
            if line.startswith('>'):
                if name: #skip flushing to dictionary if this is the first match record
                    #FIRST - flush previous hit to a nucmer_match object
                    '''if [match[1], match[0]] not in seen_matches: #if reverse has not already been read <<< ignore this step. takes too long'''
                    delta.append(Nucmer_Match(match, matchdeets))  
                    seen_matches.append([match[0], match[1]])
                    #remove alignment details of previous match
                    del matchdeets[:]
                #THEN - read in match details
                match = line.replace('>', '').split()
                if match[0] == match[1]: #skip if it's a match to itself << this didn't work! still showing matches matched to itself. Found this because after cluster step I had clusters of size = 1!
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
    if name:
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