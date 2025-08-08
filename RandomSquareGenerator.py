# 5. "Improved" Jacobson-Matthews Shuffler
# https://www.ridgerat-tech.us/latinsq/jmmethod/rlatin5JMi.py
# https://www.ridgerat-tech.us/latinsq/jmmethod/JMMethod.html#FootE
import numpy as np
from numpy.random import Generator, PCG64
#from numpy.random import Generator, MT19937
rg = Generator( PCG64() )
#rg = Generator(MT19937())
"""
This module implements the Jacobson-Mathhews "shuffler" method for
constructing pseudo-random Latin Squares, but in a new manner that does
not directly use the three-dimensional structure that is the theoretical
underpinning.

An initial Latin Square layout is required. "Move" operations are then 
applied to the corners of rectangular sub-structures, a very large number
of times, to obtain a shuffled Latin Square that is (hopefully)
"sufficiently different" from the initial pattern and "sufficiently random"
in all other respects. 

Usage:
  # Construct generator using default configuration options
  import  rlatin5JM  
  rg = rlatin5JM.rlatin( 64 )
  # Build the Latin Square
  rg.shuffle()
  # Access the result
  ls = rg.store

This code is provided Sept 10, 2023, by its author Larry Trammell under 
the CC-BY license; that means you can do pretty much anything you want with 
this as long as you cite appropriate credit to the author for providing 
it.  For details see:
        https://creativecommons.org/licenses/by/4.0/ 
"""
             
### =================================================================== 

# Local utility functions

# Make a "flip a coin" binary choice, returning True or False
def choice( ) :
  rv = rg.choice(2, shuffle=False)
  if rv==0 : return False
  return True


# Pick an integer randomly and uniformly in the range of 0 to Nrc-1 
def randix( NN ) :
  rv = rg.choice(NN, shuffle=False)
  return rv


# Apply a specified permutations to 2D storage contents.  This can 
# be a row permutation, a column permutation, or both.  Which permutations
# are applied depends on which permutation index lists are provided. 
# The rowperm and colperm are index lists of length Nrc containing index 
# values from 0 to N-1 inclusive.  
#
def  permute(storage, rowperm=None, colperm=None) :

   # Permutation not specified
   if rowperm==None and colperm==None :
      # No operation
      return
    
   # Row permutation
   if rowperm != None :
      for row in range(len(rowperm)) :
         newrow = rowperm[row]
         for col in range(len(colperm)) :
            tmp = storage[row,col]
            storage[row,col] = storage[newrow,col]
            storage[newrow,col] = tmp

   # Column permutation
   if colperm != None :
      for col in range(len(colperm)) :
         newcol = colperm[row]
         for row in range(len(rowperm)) :
            tmp = storage[row,col]
            storage[row,col] = storage[row,newcol]
            storage[row,newcol] = tmp

         
### =================================================================== 


# This class includes a 2-D storage area for the constructed Latin 
# Square object.
# 
class rlatin :
   'Defines a class that constructs and shuffles contents of a Latin square'

   # The constructor method reserves storage and size information.
   def  __init__(self, Nrc, Minmoves=0, source=None) :
      # Establish the configuration parameters
      self.Nrc = Nrc
      self.Minmoves = Minmoves
      
      # Minimum number of moves, default heuristic
      if Minmoves==0 :
         self.Minmoves = Nrc * 3

      # Various "flag" constants for content and state
      self.empty = np.short(-1)
      self.filled = np.short(1)
      self.unfilled = np.short(0)

      # Initialize statistics
      self.Moves = 0

      # Define internal state storage for J-M iteration 
      self.clearstate( )

      # Allocate and initialize storage for Latin Square
      self.store = np.zeros((Nrc,Nrc), dtype=np.short)
      self.fill(source)


   # Reset the state variables to return to the initial state
   #
   def  clearstate( self ) :
      # "Algorithm state"
      self.proper = True
      # Operating rectangle location 
      self.rowstart=self.empty
      self.colstart=self.empty
      self.rowdiag=self.empty
      self.coldiag=self.empty
      # Symbols carried from previous move
      self.symstart=self.empty      # Previous start symbol
      self.symcorn=self.empty       # Previous corner symbol
      self.symdiag=self.empty       # Previous diagonal symbol



   # Initialize the Latin Square storage using a predefined Latin 
   # pattern.  It can be an array in memory or a file.  If no pattern is 
   # specified, an artificial pattern is supplied.
   #
   def  fill(self, pattern=None) :
      # Copy from array in memory
      if isinstance(pattern,list):
         # If initializer pattern is specified, copy it
         for row in range(self.Nrc) :
            for col in range(self.Nrc) :
               self.store[row,col] = pattern[row][col]
         return

      # Copy from contents of file         
      if  type(pattern)=="file" :
         # Not yet supported
         raise NotImplementedError
        
      # Default: fill with an artificial band pattern
      for row in range(self.Nrc) :
         for col in range(self.Nrc) :
            valx = (row+col) % self.Nrc
            self.store[row,col] = valx

      # Optional permutation to break location correlations
      self.rcpermute( )


   # Apply a randomized row-column permutations to contents of
   # the Latin Square storage area.  
   #
   def  rcpermute(self) :
      rowix = list(range(self.Nrc))
      rg.shuffle(rowix)
      colix = list(range(self.Nrc))
      rg.shuffle(colix)
      permute( self.store, rowperm=rowix, colperm=colix )
      
   def symbol_permute(self):
      symbol_perm = np.arange(self.Nrc)
      rg.shuffle(symbol_perm)
      symbol_map = {i: symbol_perm[i] for i in range(self.Nrc)}
      self.store = np.vectorize(lambda x: symbol_map[x])(self.store)
      
   def isotopy_shuffle(self):
      self.rcpermute()
      self.symbol_permute()

   # Perform a search in the specified row to find the specified symbol... but
   # but SKIP the excluded column location if given.  Return the column of the   
   # found symbol. If none (a pathological case!) return the special value "empty."
   # If exclusion is empty, any column can provide match.
   #
   def rowmatch( self, row, sym, exclude=-1 ) :
      found = self.empty

      # Search columns of specified row
      for col in range( self.Nrc ):

         # Omit the excluded column location if any location 
         if col == exclude :
            continue

         # Check other columns for a match
         if  self.store[row,col] == sym :
            found = col
            break
 
      # Return the result of search
      return found     


   # Perform a search in the specified column to find the specified symbol... 
   # but SKIP the excluded row location if given.  Return the row of the   
   # found symbol. If none (a pathological case!) return the special value "empty."
   # "empty."  If exclusion is empty, any row can provide the match.
   #
   def colmatch( self, col, sym, exclude=-1 ) :
      found = self.empty

      # Search rows of specified column
      for row in range( self.Nrc ):

         # Omit the corner symbol location 
         if row == exclude :
            continue

         # Check other rows for a match
         if  self.store[row,col] == sym :
            found = row
            break
 
      # Return the result of search
      return found     


   # Establish a start location and perform an initial move when
   # the the algorithm state is proper.
   def  JMproperMove( self ) :

      # Select row and column location at random, skipping current
      # start row and column
      self.rowstart = randix( self.Nrc )
      self.colstart = randix( self.Nrc )
      self.symstart = self.store[self.rowstart, self.colstart]

      # Select different corner symbol
      self.symcorn = randix( self.Nrc - 1 )
      if  self.symcorn == self.symstart :
         self.symcorn = self.Nrc - 1

      # Find the two corner locations that contain this symbol
      self.rowdiag = self.colmatch( self.colstart, self.symcorn )
      self.coldiag = self.rowmatch( self.rowstart, self.symcorn )
      newsym = self.store[self.rowdiag,self.coldiag]

      # Locate and test the diagonal cell for special case 1
      if  newsym == self.symstart :

         # Case 1 - Rotate assignments in corners and stay in proper state
         self.store[ self.rowstart, self.colstart ] = self.symcorn
         self.store[ self.rowdiag, self.coldiag ]   = self.symcorn
         self.store[ self.rowstart, self.coldiag ]  = self.symstart
         self.store[ self.rowdiag, self.colstart ]  = self.symstart
         self.Moves += 1
         self.clearstate( )

      # Otherwise, apply case 2 and the state becomes improper
      else :

         # Case 2 - Reassign symbols at three cells unambiguously 
         self.store[ self.rowstart, self.colstart ] = self.symcorn
         self.store[ self.rowstart, self.coldiag ]  = self.symstart
         self.store[ self.rowdiag, self.colstart ]  = self.symstart
         self.store[ self.rowdiag, self.coldiag ]   = self.empty

         # Preserve conflict information for empty diagonal cell
         self.symdiag = newsym
         self.Moves += 1
         self.proper = False

   # End of proper move 

        
   # Set up the next operating rectangle at the conflicted cell carried
   # from the previous move.  Apply an "improper move" and check the
   # resulting state. 
   #
   def  JMimproperMove( self ) :

      # Locate the cells with conflicting labels.
      rowconf1 = self.rowstart
      rowconf2 = self.colmatch( self.coldiag, self.symstart, self.rowstart )
      colconf1 = self.colstart
      colconf2 = self.rowmatch( self.rowdiag, self.symstart, self.colstart )

      # Randomly pick corner row and column such that at least one of them is new.
      option = randix(3)
      if  option==2 :
         newrow = rowconf2; newcol = colconf2
      if  option==1 :
         newrow = rowconf2; newcol = colconf1     
      if  option==0 :
         newrow = rowconf1; newcol = colconf2
      symnew = self.store[newrow,newcol]

      # Check for the first special case that returns to proper state
      if symnew == self.symdiag :
         # Case 3 applies.  Perform the move operations.
         self.store[self.rowdiag, self.coldiag] = self.symcorn
         self.store[self.rowdiag, newcol] = self.symdiag
         self.store[newrow, self.coldiag] = self.symdiag
         self.store[newrow, newcol] = self.symstart

         # The state becomes proper.
         self.Moves += 1
         self.clearstate( )

      # Check for the second special case that returns to proper state
      elif symnew == self.symcorn :
         # Case 4 applies.  Perform the move operations.
         self.store[self.rowdiag, self.coldiag] = self.symdiag
         self.store[self.rowdiag, newcol] = self.symcorn
         self.store[newrow, self.coldiag] = self.symcorn
         self.store[newrow, newcol] = self.symstart
 
         # The state becomes proper.
         self.Moves += 1
         self.clearstate( )

      # For any other cases, conflicts change but state remains improper.
      else :
         # Case 5 applies.  Perform the move operations.

         # Pick a symbol for resolving the old conflict.
         if  choice() :

            # Resolve with previous diagonal symbol
            self.store[self.rowdiag, self.coldiag] = self.symdiag
            self.store[newrow, self.coldiag] = self.symcorn
            self.store[self.rowdiag, newcol] = self.symcorn

            # Update the conflict information
            self.store[newrow,newcol] = self.empty
            self.symdiag = symnew
            swap = self.symcorn
            self.symcorn = self.symstart
            self.symstart = swap

         else :

            # Resolve with previous corner symbol
            self.store[self.rowdiag, self.coldiag] = self.symcorn
            self.store[newrow, self.coldiag] = self.symdiag
            self.store[self.rowdiag, newcol] = self.symdiag

            # Update the conflict information
            self.store[newrow,newcol] = self.empty
            self.symcorn = self.symstart
            self.symstart = self.symdiag
            self.symdiag = symnew
             
         # Move to new location
         self.rowstart = self.rowdiag
         self.colstart = self.coldiag
         self.rowdiag = newrow
         self.coldiag = newcol

         # The state remains improper after Case 5 move.
         self.Moves += 1

      # End of improper move



   # Apply at least the configured minimum number of shuffling moves to
   # scramble the contents of the Latin Square storage, leaving the result
   # in a proper state. 
   #
   def  shuffle( self ) :

      self.Moves = 0

      while self.Moves < self.Minmoves or not self.proper :

         if  self.proper : 
            self.JMproperMove( )

         else :
            self.JMimproperMove( )

      # End of shuffler method

       
# End of random Latin Square class using J-M iterations