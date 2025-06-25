#!/usr/bin/env python
# coding: utf-8

# In[51]:


from itertools import permutations
import numpy as np
import random


# ### Misc Functions

# In[52]:


def print_square(square):
    '''
    Print squre row by row for better readability.
    '''
    for row in square:
        print(row)
    print()


# In[53]:


def switch_first_row(square):
    '''
    Swaps the first row with each row. 
    '''
    n = len(square)
    results = []

    for i in range(1,n):
        # copy orginal square. Deepcopy is safer for list of list
        new_square = [list(row) for row in square]
        
        # Swap row 0 with row i
        new_square[0], new_square[i] = new_square[i], new_square[0]
        results.append(new_square)
    
    return results


# ### 1. Generating Latin Squares

# In[54]:


def generate_all_squares(n: int) -> list:
    '''
    Generates all Latin squares of order n.
    '''
    all_rows = list(permutations(range(1,n+1)))
    squares = []

    def backtrack(square, depth):
        if depth == n:
            squares.append([list(row) for row in square])
            return
        
        for row in all_rows:
            if all(row[col] not in used_in_col[col] for col in range(n)):
                square.append(row)
                for col in range(n):
                    used_in_col[col].add(row[col])

                backtrack(square, depth + 1)

                for col in range(n):
                    used_in_col[col].remove(row[col])
                square.pop()
    used_in_col = [set() for _ in range(n)]

    backtrack([],0)
    return squares

def generate_all_squares_standardized(n: int) -> list:
    '''
    Generates all standardized Latin squares of order n,
    fixing the first row to [1, 2, ..., n].
    '''
    all_rows = list(permutations(range(1, n + 1)))
    squares = []

    def backtrack(square, depth):
        if depth == n:
            squares.append(square[:])
            return

        for row in all_rows:
            if all(row[col] not in used_in_col[col] for col in range(n)):
                square.append(row)
                for col in range(n):
                    used_in_col[col].add(row[col])

                backtrack(square, depth + 1)

                for col in range(n):
                    used_in_col[col].remove(row[col])
                square.pop()

    first_row = tuple(range(1, n + 1))
    used_in_col = [set() for _ in range(n)]
    for col in range(n):
        used_in_col[col].add(first_row[col])  # track first row values
        
    backtrack([first_row], 1)
    return squares
    


# In[55]:


def generate_all_squares_reduced(n:int) -> list:
    '''
    Generates all reduced Latin squares of order n, fixing the first row and first col to [1, 2, ..., n].
    Uses numpy to optimize row/col checks, removes python list overhead(dynamic allocation and resizing),
    and has better memory layout for copying and modifying. 
    Returns a list of numpy arrays, type np.ndarray.
    '''

    results = []
    square = np.zeros((n, n), dtype=int)
    
    # Fix first row and column to [1, 2, ..., n]
    square[0, :] = np.arange(1, n + 1)
    square[:, 0] = np.arange(1, n + 1)
    
    def backtrack(row, col):
        if row == n:
            # copy is a numpy function, more efficient than a deepcopy or manual loops
            results.append(square.copy())
            return
        if col == n:
            backtrack(row + 1, 1)  # move to next row, skip col 0 (fixed)
            return

        for val in range(1, n + 1):
            # Check if val is already in this row or column
            
            # square[row, :col] creates a 1D view, which makes "not in" more efficient
            if val not in square[row, :col] and val not in square[:row, col]:
                square[row, col] = val
                backtrack(row, col + 1)
                square[row, col] = 0  # undo move (backtrack)

    backtrack(1, 1)  # start from cell (1,1)
    return results


# ### 2. Standardize & Reduce Squares

# In[56]:


def standardize(square: list) -> list:
    '''
    Standardize a Latin square so that the first row is [1, 2, ..., n].
    '''
    first_row = square[0]

    # Creates a mapping dict for standardization of square. ex {3:1, 2:3, 1:1} (3's change to 1.. )
    row_mapping = {val: i + 1 for i, val in enumerate(first_row)}

    # Apply the mapping to the entire square
    standardized_square = [[row_mapping[val] for val in row] for row in square]
    '''
    new_square = []
    for row in square:
        new_row = []
        for val in row:
            new_row.append(row_mapping[val])
        new_square.append(new_row)
    '''

    return standardized_square


# In[57]:


def reduce_square(square: list) -> list:
    '''
    Reduces latin square so first row = 1,2,..,n and first col = 1,2,..,n
    '''
    # Step 1: Standardize first row
    square = standardize(square)
    
    # Step 2: Permute rows so first col is 1,2,..,n
    square_reduced = sorted(square, key=lambda row: row[0])     
    # Sorts the rows of the squares based on the first value of each row same as: 
    '''
    def get_first_element(row):
        return row[0]
    '''

    return square_reduced


# ### 3. Permutation Functions (NOT USING THESE)

# In[58]:


def row_col_permutations(square: list) -> list:
    '''
    Generate all Latin squares by permuting rows and columns of the given square.
    Returns a list of all Latin squares generated for row/column permutations of passed square.
    '''
    n = len(square)
    perms = set()

    for row_perm in permutations(range(n)):
        # Apply row permutation
        permuted_rows = [square[i] for i in row_perm]

        for col_perm in permutations(range(n)):
            # Apply column permutation to each row
            permuted_square = [
                [row[j] for j in col_perm]
                for row in permuted_rows
            ]
            # Convert to tuple of tuples for hashing
            perms.add(tuple(tuple(row) for row in permuted_square))

    # Return list of unique permutations
    return [ [list(row) for row in square] for square in perms ]


# In[59]:


def symbol_permutations(square:list , mapping: dict) -> list:
    '''
    Applies a symbol permutation to a Latin square.

    parameters:
        symbol_perm (dict): A mapping from original symbol to new symbol.
    '''
    return [[mapping[val] for val in row] for row in square]


# In[60]:


def all_symbol_permutations(square:list) -> list:
    '''
    Generates all unique Latin squares via symbol permutations.
    '''
    n = len(square)
    perms = permutations(range(1, n + 1))
    seen = set()
    unique_squares = []

    for p in perms:
        mapping = {old: new for old, new in zip(range(1, n + 1), p)}
        permuted = symbol_permutations(square, mapping)
        flat = tuple(tuple(row) for row in permuted)
        if flat not in seen:
            seen.add(flat)
            unique_squares.append(permuted)

    return unique_squares


# ### 4. Find Classes

# In[61]:


def isotopy(n:int):
    '''
    Finds isotopy classes (row/col/symbol permutations).
    '''
    # Step 1. Generate latin squares of order n
    squares = generate_all_squares(n)
    #squares = generate_all_squares_reduced(n)

    
    # covert to tuples to store in set 
    remaining = set(tuple(map(tuple, square)) for square in squares)  # hashable form (immutable)
    classes = []

    row_perms = list(permutations(range(n))) # 0 - 3
    col_perms = list(permutations(range(n))) # 0 - 3
    sym_perms = list(permutations(range(1, n + 1))) # 1 - 4


    while remaining:
        # Step 2: Pick an arbitrary square Sj from the remaining ones
        Sj = random.choice(list(remaining))
        Sj_matrix = [list(row) for row in Sj]

        isotopy_class = set() # set enforces unique squares 

        # Step 3. Generate its  class by row/column/symbol permutations 
        for r_perm in row_perms:
            for c_perm in col_perms:
                for s_perm in sym_perms:
                    # Apply row permutation
                    # ex. r_perm = (2,0,3,1) -> row 0 becomes 2, row 1 becomes 0...
                    row_permuted = [Sj_matrix[r_perm[i]] for i in range(n)]

                    # Apply column permutation
                    # for each row in row_permuted, for each new column position c_perm, take value at original position 
                    col_permuted = [[row[c_perm[j]] for j in range(n)] for row in row_permuted]

                    # Create symbol mapping (dict) from original symbol to permuted symbol
                    sym_map = dict(zip(range(1, n + 1), s_perm))
                    # Apply symbol permutation
                    final = [[sym_map[val] for val in row] for row in col_permuted]

                    # convert back to hashable form and add to set
                    isotopy_class.add(tuple(tuple(row) for row in final))

        # Step 4. Record the class
        # map turns square back to list of list, then create a list of squares and add to classes
        classes.append([list(map(list, square)) for square in isotopy_class])

        # Remove entire isotopy class from remaining
        remaining -= isotopy_class
    
    
    # return sizes also 
    class_sizes = [len(c) for c in classes]
    return classes, class_sizes


# In[62]:


def equivalence(n:int):
    '''
    Finds row/col equivalence classes.
    '''
    # Step 1: Generate latin squares of order n
    squares = generate_all_squares(n)
    #squares = generate_all_squares_reduced(n)
    
    # Standardize all squares. Removes symbol permutations 
    for i in range(len(squares)):
        squares[i] = standardize(squares[i])
        
    # covert to tuples to store in set 
    remaining = set(tuple(map(tuple, square)) for square in squares)  # hashable form
    print(f'{len(remaining)} unique squares after standardizing')

    classes = []
    
    row_perms = list(permutations(range(n))) # 0 - 3
    col_perms = list(permutations(range(n))) # 0 - 3

    while remaining:
        # Step 2: Pick an arbitrary square Sj from the remaining ones
        Sj = random.choice(list(remaining))
        Sj_matrix = [list(row) for row in Sj]  # convert back to list of lists
        
        eq_class = set() # enforces unique squares 
                
        # Step 3. Generate its equivalence class by row/column permutations
        for r_perm in row_perms:
            for c_perm in col_perms:
                # Apply row permutation
                # ex. r_perm = (2,0,3,1) -> row 0 becomes 2, row 1 becomes 0...
                row_permuted = [Sj_matrix[r_perm[i]] for i in range(n)]

                # Apply column permutation
                # for each row in row_permuted, for each new column position c_perm, take value at original position 
                col_permuted = [[row[c_perm[j]] for j in range(n)] for row in row_permuted]
                
                # Step 4. Standardize permuted square
                # removes symbol symmetry from consideration
                final = standardize(col_permuted) 
                
                eq_class.add(tuple(tuple(row) for row in final))

        # Step 5. Record the class
        # map turns square back to list of list, then create a list of squares and add to classes
        classes.append([list(map(list, square)) for square in eq_class])
        remaining -= eq_class # Remove this entire class from the remaining pool        
        
    # return sizes also
    class_sizes = [len(c) for c in classes]
    return classes, class_sizes


# In[63]:


def find_reduced_in_isotopy_class(classes:list):
    
    #Find the amount of reduced squares in each isotopy class passed into the function.    
    
    # find size of squares in isotopy class
    n = len(classes[0][0])
    reduced_squares_in_class = []

    # Find number of reduced squares in each isotopy class
    for iso_class in classes:
        # Reduce all squares.
        for i in range(len(iso_class)):
            iso_class[i] = reduce_square(iso_class[i])
        
        Sj = random.choice(iso_class)
    
        # Gets n squares where first row is interchanged by another
        row_interchange = switch_first_row(Sj)
    
        reduced_squares = set()
        
        # Apply all col permutations to each of the n squares in row_interchange
        col_perms = list(permutations(range(n))) # 0 - 3
        for square in row_interchange:
            for perm in col_perms:
                # Apply column permutation to square
                # for each row in row_interchange, for each new column position c_perm, take value at original position 
                col_permuted = [[row[perm[j]] for j in range(n)] for row in square]
                
                # Call reducing function and add to class
                reduced = reduce_square(col_permuted)
                
                reduced_squares.add(tuple(tuple(row) for row in reduced))
                
        reduced_squares_in_class.append(list(map(list, reduced_squares)))
    
     # return number of reduced squares also
    num_reduced_sizes = [len(c) for c in reduced_squares_in_class]
    return reduced_squares_in_class, num_reduced_sizes 

# ORIGINAL
'''
def find_reduced_in_isotopy_class(n):
    # Step 1: Generate latin squares of order n
    squares = generate_all_squares(n)   

    # Reduce all squares.
    for i in range(len(squares)):
        squares[i] = reduce_square(squares[i])
        
    Sj = random.choice(squares)
    
    # Gets n squares where first row is interchanged by another
    row_interchange = switch_first_row(Sj)

    reduced_squares = set()
    
    # Apply all col permutations to each of the n squares in row_interchange
    col_perms = list(permutations(range(n))) # 0 - 3
    for square in row_interchange:
        for perm in col_perms:
            # Apply column permutation to square
            # for each row in row_interchange, for each new column position c_perm, take value at original position 
            col_permuted = [[row[perm[j]] for j in range(n)] for row in square]
            
            # Call reducing function and add to class
            reduced = reduce_square(col_permuted)
            
            reduced_squares.add(tuple(tuple(row) for row in reduced))
            
    return list(map(list, reduced_squares)), len(reduced_squares)
'''


# #### 4. Test Cases

# Find classes

# In[71]:


#eq_classes_3x3, eq_class_sizes_3x3 = equivalence(3)
#print("Number of equivalence classes:", len(eq_classes_3x3))
#print("Sizes of classes:", eq_class_sizes_3x3)
print("3x3: ")
isotopy_classes_3x3, isotopy_class_sizes_3x3 = isotopy(3)
print("Number of isotopy classes:", len(isotopy_classes_3x3))
print("Sizes of classes:", isotopy_class_sizes_3x3)

reduced_squares_3x3, num_reduced_squares_3x3 = find_reduced_in_isotopy_class(isotopy_classes_3x3)
print()
print("Number of reduced squares in each isotopy class:", num_reduced_squares_3x3)


# In[ ]:


#eq_classes_4x4, eq_class_sizes_4x4 = equivalence(4)
#print("Number of equivalence classes:", len(eq_classes_4x4))
#print("Sizes of classes:", eq_class_sizes_4x4)
print("4x4:")
isotopy_classes_4x4, isotopy_class_sizes_4x4 = isotopy(4)
print("Number of isotopy classes:", len(isotopy_classes_4x4))
print("Sizes of classes:", isotopy_class_sizes_4x4)

reduced_squares_4x4, num_reduced_squares_4x4 = find_reduced_in_isotopy_class(isotopy_classes_4x4)
print()
print("Number of reduced squares in each isotopy class:", num_reduced_squares_4x4)


# In[73]:


#eq_classes_5x5, eq_class_sizes_5x5 = equivalence(5)
#print("Number of equivalence classes:", len(eq_classes_5x5))
#print("Sizes of classes:", eq_class_sizes_5x5)
#print()
'''for eq_class in eq_classes_5x5:
    for _ in range(3):
        print(eq_class[_])
    print('New Class -------')'''
    

print("5x5: ")
isotopy_classes_5x5, isotopy_class_sizes_5x5 = isotopy(5)
print("Number of isotopy classes:", len(isotopy_classes_5x5))
print("Sizes of classes:", isotopy_class_sizes_5x5)

reduced_squares_5x5, num_reduced_squares_5x5 = find_reduced_in_isotopy_class(isotopy_classes_5x5)
print()
print("Number of reduced squares in each isotopy class:", num_reduced_squares_5x5)


# ### 5. "Improved" Jacobson-Matthews Shuffler

# In[67]:


import RandomSquareGenerator # Imported from different file 

latin = RandomSquareGenerator.rlatin(7)
latin.shuffle()
ls = latin.store
print_square(ls) 


# ### 6. Find Classes of Reduced Squares

# In[90]:


def build_reduced_isotopy_classes(reduced_squares):
    n = len(reduced_squares[0])  # Find Order of the Latin square
    reduced_set = {
        tuple(tuple(row) for row in sq): i
        for i, sq in enumerate(reduced_squares)
    } # Create set of all reduced squares with index

    # To keep track of which square belongs to which class
    square_to_class = {} # dict that stores final classification 
    class_id = 0
    seen = set()

    # Reduced squares so row one and col one will be equal 
    row_perms = list(permutations(range(1, n)))  # skip row 0
    col_perms = list(permutations(range(1, n)))  # skip col 0


    for i, square in enumerate(reduced_squares):
        # convert square to hashable key
        key = tuple(tuple(row) for row in square)
        if key in seen:
            continue  # Already assigned to a class, keep going

        # New class
        iso_class_keys = set()
        for r_perm in row_perms:
            # dont permute first col
            row_ix = [0] + list(r_perm) 
            row_perm = [square[row_ix[i]] for i in range(n)]

            for c_perm in col_perms:
                # dont permute first col
                col_ix = [0] + list(c_perm) # ex [0,3,1,2]
                col_perm = [[row[col_ix[j]] for j in range(n)] for row in row_perm]

                # Reduce square, turn to hashable key
                reduced = reduce_square(col_perm)
                reduced_key = tuple(tuple(row) for row in reduced)

                # If new square is in the list of all order n reduced squares 
                if reduced_key in reduced_set: 
                    iso_class_keys.add(reduced_key)

        # Assign class ID to all reduced squares in the found class
        # k is a squre key 
        for k in iso_class_keys:
            # find the index associated to the original square and add it to the dict
            ind = reduced_set[k]
            square_to_class[ind] = class_id
            seen.add(k)

        class_id += 1

    class_list = {}
    for ind, cid in square_to_class.items():
        # if cid is not in the class_list create new dic entry with value [] (#cid: [])
        # add square to the dictionary for final mapping
        class_list.setdefault(cid, []).append(reduced_squares[ind])

    return class_list, square_to_class


# In[91]:


reduced_squares = generate_all_squares_reduced(5)
class_list, square_to_class = build_reduced_isotopy_classes(reduced_squares)

print(f"Found {len(class_list)} reduced isotopy classes.")
for cid, squares in class_list.items():
    print(f"Class {cid}: {len(squares)} squares")
    #for square in squares:
    #    print_square(square)


# In[70]:


# Terminal command to convert notebook to python file
# jupyter nbconvert --to python LatinSquares.ipynb

