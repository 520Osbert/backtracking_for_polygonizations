import numpy as np
import matplotlib.pyplot as plt
import random

'''
fig = plt.figure()
pts_x = np.array([0, 1, 2, 3] * 4)
pts_y = np.array([0] * 4 + [1] * 4 + [2] * 4 + [3] * 4)
fig.add_subplot(1,1,1).scatter(pts_x, pts_y)
plt.show()
'''

# Backtracking pseudo-code:

'''
def backtracking_search(csp) -> solution or failure
    return backtrack({}, csp)

def backtrack(assignment, csp) -> solution or failure
    if assignment is complete:
        return assignment
    var = select_unassigned_variables(csp)
    for each value in order_domain_values(var, assignment, csp)
        if value is consistent with assignment:
            add {var=value} to assignment
            result = backtrack(assignment, csp)
            if result != failure:
                return result
            remove {var=value} from assignment
    return failure
'''

class BT_polygonalization:
    def __init__(self, row, col, pts):
        '''
        Initialization of BT_polygonalization class, taking in a column
        number and contructing a dictionary that represents a n*n graph
        
        Parameters
        ----------
        row : int
            number of rows
        
        col : int
            number of columns
        
        pts : list of str 
            names of the points in the graph; its length must be row*col
            
        '''
        
        self.pts      = pts
        self.row      = row
        self.col      = col
        self.temp     = [pt for pt in pts] # temporary list for point names, eventaully empty
        self.ptDict   = dict()
        self.edgeDict = dict()
        self.polygons = int()
        
        for row in range(self.row):
            for col in range(self.col):
                self.ptDict[self.temp[0]] = (row, col)
                self.edgeDict[self.temp.pop(0)] = ['', '']
                
    def check_intersection(self, pts1, pts2):
        '''
        Check if two input line segments intersect with each other
        
        Parameters
        ----------
        pts1, pts2 : (str, str)
            e.g. ('0', '1')
            tuple of names of the two ending points of a line segment
        
        '''
        def onSeg(p, q, r):
            if q[0] <= max(p[0], r[0]) and q[0] >= min(p[0], r[0]) and q[1] <= min(p[1], r[1]) and q[1] >= min(p[1], r[1]):
                return True
            return False
        
        def direction(p, q, r):
            val = (q[1] - p[1]) * (r[0] - q[0]) - (q[0] - p[0]) * (r[1] - q[1])
            if val == 0:
                return 0
            return 1 if val > 0 else -1 # 1 for clockwise, -1 for counter-clockwise
        
        p1 = self.ptDict[pts1[0]]
        q1 = self.ptDict[pts1[1]]
        p2 = self.ptDict[pts2[0]]
        q2 = self.ptDict[pts2[1]]
        
        o1 = direction(p1, q1, p2)
        o2 = direction(p1, q1, q2)
        o3 = direction(p2, q2, p1)
        o4 = direction(p2, q2, q1)

        if o1 != o2 and o3 != o4:
            return True
        
        if o1 == 0 and onSeg(p1, p2, q1):
            return True
        if o2 == 0 and onSeg(p1, q2, q1):
            return True
        if o3 == 0 and onSeg(p2, p1, q2):
            return True
        if o4 == 0 and onSeg(p2, q1, q2):
            return True
            
        return False
    
    def check_consistency(self, assignment):
        '''
        Check if the current assignment is consistent

        Parameters
        ----------
        assignment : [str]
            max(len(assignment)) = row * col + 1
            e.g. ['0', '1', '2', ...]
            list of vertex names in order of polygonalization
                    
        '''
        if len(assignment) < self.row * self.col + 1:
            if len(assignment) - len(set(assignment)) > 0:
                return False
            if len(assignment) >= 3:
                for i in range(len(assignment)-3):
                    if self.check_intersection([assignment[-2], assignment[-1]], [assignment[i], assignment[i+1]]):
                        return False
        
        elif len(assignment) == self.row * self.col + 1:
            if (len(assignment) - len(set(assignment)) != 1) or (assignment[0] != assignment[-1]):
                return False
            for i in range(1, len(assignment)-3):
                if self.check_intersection([assignment[-2], assignment[-1]], [assignment[i], assignment[i+1]]):
                    return False
        else:
            return False
        return True
    
    def check_completeness(self, assignment):
        if (len(assignment) == self.row * self.col + 1) and self.check_consistency(assignment):
            return True
        return False
    
    def find_next_pts(self, assignment):
        '''
        Find a list of points which are legitimate to form the next edge

        assignment : [str]
            max(len(assignment)) = row * col
            e.g. ['0', '1', '2', ...]
            list of vertex names in order of polygonalization
            
        '''
        result = []
        current_pt = assignment[-1]
        current_pt_x = self.ptDict[current_pt][0]
        current_pt_y = self.ptDict[current_pt][1]
        
        for pt in self.pts:
            pt_x = self.ptDict[pt][0]
            pt_y = self.ptDict[pt][1]
            
            if pt == current_pt:
                continue
            elif (np.abs(current_pt_x-pt_x) >= 3) or (np.abs(current_pt_y-pt_y) >= self.col-1):
                continue
            elif (current_pt_x == pt_x) and (np.abs(current_pt_y-pt_y) >= 2):
                continue
            elif (current_pt_y == pt_y) and (np.abs(current_pt_x-pt_x) >= 2):
                continue
            elif (np.abs(current_pt_x-pt_x) % 2 == 0) and (np.abs(current_pt_y-pt_y) % 2 == 0):
                continue
            elif ((current_pt_x == 0) or (current_pt_x == self.row-1)) and ((pt_y == 0) or (pt_y == self.col-1)) and (current_pt_x != pt_x) and (current_pt_y != pt_y):
                continue
            elif ((current_pt_y == 0) or (current_pt_y == self.col-1)) and ((pt_x == 0) or (pt_x == self.row-1)) and (current_pt_x != pt_x) and (current_pt_y != pt_y):
                continue
            else:
                result.append(pt)
        return result
    
    def BT_search(self):
        '''
        The complete backtracking search algorithm for finding number of polygonalization in row by col grids

        assignment : [str]
            max(len(assignment)) = row * col + 1
            e.g. ['0', '1', '2', ...]
            list of vertex names in order of polygonalization
            
        self.polygons : int
            number of polygons that (row * col) vertices can form

        '''
        
        # initialize count polygons to 0
        self.polygons = 0
        
        def backtrack(assignment):
            if self.check_completeness(assignment):
                self.polygons += 1
                
                fig = plt.figure()
                ax = fig.add_subplot(1,1,1)
                xs = [self.ptDict[i][0] for i in assignment]
                ys = [self.ptDict[i][1] for i in assignment]
                ax.plot(xs, ys)
                plt.show()
                
                pass
            
            next_pts = self.find_next_pts(assignment)
            current_pt = assignment[-1]
            
            for next_pt in next_pts:
                if self.check_consistency(assignment):
                    self.edgeDict[current_pt][1] = next_pt
                    self.edgeDict[next_pt][0]    = current_pt
                    assignment.append(next_pt)
                    
                    backtrack(assignment)
                    
                    self.edgeDict[current_pt][1] = ''
                    self.edgeDict[next_pt][0]    = ''
                    assignment.pop()
        
        first_pt = random.choice(list(self.edgeDict.keys()))
        backtrack([first_pt])
        self.polygons = int(self.polygons / 2)

		
'''
bt = BT_polygonalization(3,3, [str(i) for i in range(9)])
bt.BT_search() 							## Run backtracking search
print (bt.polygons) 						## Print out the # of polygonizations
'''
