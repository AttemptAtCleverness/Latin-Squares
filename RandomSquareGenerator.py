# 5. "Improved" Jacobson-Matthews Shuffler
# https://www.ridgerat-tech.us/latinsq/jmmethod/rlatin5JMi.py
# https://www.ridgerat-tech.us/latinsq/jmmethod/JMMethod.html#FootE
import numpy as np
from numpy.random import Generator, PCG64
#from numpy.random import Generator, MT19937

'''
This module implements the Jacobson-Mathhews "shuffler" method 
for constructing pseudo-random Latin Squares, following the implementation
of the rlatin5JM module, but with the addition of reverse-indexing 
arrays that enable O(1) lookup to discover locations of indicated 
symbols within the Latin Square storage area.  This reduces algorithm 
complexity at the expense of considerably more storage and more
general overhead, but regardless, might be highly beneficial for 
large Latin Square cases. 

Usage:
  # Construct generator using default configuration options
  import  rlatin5JMi 
  rg = rlatin5JMi.rlatin( 64 )
  # Build the Latin Square
  rg.shuffle()
  # Access the result
  ls = rg.store

This code is provided Sept. 10, 2023 by its author Larry Trammell under 
the CC0 license; that means you can do pretty much anything you want with 
this as long as you cite appropriate credit to the author for providing 
it.  For details see:
        https://creativecommons.org/licenses/by/4.0/ 
'''

rng = Generator(PCG64())
#rg = Generator(MT19937())

# Local randomizing functions

# Make a "flip a coin" binary choice, returning True or False
def choice( ) :
   rv = rng.choice(2, shuffle=False)
   if rv==0 : 
     return False
   return True


# Pick an integer randomly and uniformly in the range of 0 to Nrc-1 
def randix( NN ) :
  rv = rng.choice(NN, shuffle=False)
  return rv

# Apply a specified permutations to 2D storage contents.  This can 
# be a row permutation, a column permutation, or both.  Which permutations
# are applied depends on which permutation index lists are provided. 
# The rowperm and colperm are index lists of length Nrc containing index 
# values from 0 to N-1 inclusive.
# 
# Note that doing this INVALIDATES all of the row and column reverse
# indexing data.  YOU MUST REBUILD that to further use this method!  
#
def  permute(storage, rowperm=None, colperm=None) :

   # Permutation not specified
   if rowperm is None and colperm is None :
      # No operation
      return
    
   # Row permutation
   if rowperm is not None :
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
   def  __init__(self, Nrc, Minmoves=0, source=None, rng=None) :
      # Establish the configuration parameters
      self.Nrc = Nrc
      self.Minmoves = Minmoves
      self.rng = rng or Generator(PCG64())
      
      # Number of moves default heuristic
      if Minmoves==0 :
         self.Minmoves = 3*Nrc

      # Various "flag" constants for content and state
      self.empty = np.short(-1)
      self.filled = np.short(1)
      self.unfilled = np.short(0)

      # Initialize statistics
      self.Moves = 0

      # State
      self.clearstate()

      # Supplemental reverse-indexing arrays to facilitate symbol search
      self.rowx = 0
      self.colx = 0

      # Allocate and initialize storage
      self.store = np.zeros((Nrc,Nrc), dtype=np.short)
      self.fill(source)
      #print("Initial")
      #print("Store\n",self.store)
      #print("rowx \n",self.rowx)
      #print("colx\n",self.colx)


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
   # This invalidates any existing data in the reverse indexing arrays.
   # The rebuild is applied automatically.
   #
   def  fill(self, pattern=None) :
      # Copy from array in memory
      if  isinstance(pattern,np.ndarray):
         # If initializer pattern is specified, copy it
         for row in range(self.Nrc) :
            for col in range(self.Nrc) :
               self.store[row,col] = pattern[row,col]
         return

      # Copy from contents of file         
      if  type(pattern)=="file" :
         # Not yet supported
         return
     
      # Default: fill with an artificial band pattern
      for row in range(self.Nrc) :
         for col in range(self.Nrc) :
            valx = (row+col) % self.Nrc
            self.store[row,col] = valx

      # Optional permutation to break location correlations
      self.rcpermute( )

      # The reverse-indexing arrays must be rebuilt!
      self.buildindex( )



   # Build or rebuild the reverse-lookup indexes for row and column
   # searches.  The contents of the current storage areas MUST BE
   # valid Latin Square data.  
   #
   def  buildindex(self) :
      # Clear the indexing arrays
      self.rowx = np.full((self.Nrc,self.Nrc), self.empty, dtype=np.short)
      self.colx = np.full((self.Nrc,self.Nrc), self.empty, dtype=np.short)

      # Go through storage item by item and build indexing
      for row in range(self.Nrc) :
         for col in range(self.Nrc) :
            sym = self.store[row,col]
            self.rowx[row,sym] = col
            self.colx[col,sym] = row
   


   # Apply a randomized row-column permutations to contents of
   # the Latin Square storage area.  This completely destroys any prior
   # row-column indexing!  Rebuild them to continue using any methods 
   # of this class.
   #
   def  rcpermute(self) :
      rowix = list(range(self.Nrc))
      self.rng.shuffle(rowix)
      colix = list(range(self.Nrc))
      rng.shuffle(colix)
      permute( self.store, rowperm=rowix, colperm=colix )



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
      self.rowdiag = self.colx[ self.colstart, self.symcorn ]
      self.coldiag = self.rowx[ self.rowstart, self.symcorn ]
      newsym = self.store[self.rowdiag,self.coldiag]

      # Locate and test the diagonal cell for special case 1
      if  newsym == self.symstart :

         # Case 1 
         # Assign the new cell values
         self.store[ self.rowstart, self.colstart ] = self.symcorn
         self.store[ self.rowdiag, self.coldiag ]   = self.symcorn
         self.store[ self.rowstart, self.coldiag ]  = self.symstart
         self.store[ self.rowdiag, self.colstart ]  = self.symstart

         # Adjust the index arrays to reflect the new symbol locations
         self.colx[ self.colstart, self.symcorn]  = self.rowstart
         self.rowx[ self.rowstart, self.symcorn]  = self.colstart
         self.colx[ self.coldiag, self.symcorn]   = self.rowdiag
         self.rowx[ self.rowdiag, self.symcorn]   = self.coldiag
         self.colx[ self.coldiag, self.symstart]  = self.rowstart
         self.rowx[ self.rowstart, self.symstart] = self.coldiag
         self.colx[ self.colstart, self.symstart] = self.rowdiag
         self.rowx[ self.rowdiag, self.symstart]  = self.colstart

         # Clear state areas and return to proper state  
         self.Moves += 1
         self.clearstate( )

      # Otherwise, apply case 2 and the state becomes improper
      else :

         # Case 2 
         # Corner cells have conflicts. Clear them.
         self.store[ self.rowstart, self.coldiag ] = self.empty
         self.store[ self.rowdiag, self.colstart ] = self.empty
         self.colx[ self.coldiag,  self.symcorn ]  = self.empty
         self.rowx[ self.rowstart, self.symcorn ]  = self.empty
         self.colx[ self.colstart, self.symcorn ]  = self.empty
         self.rowx[ self.rowdiag,  self.symcorn ]  = self.empty

         # The diagonal cell is the new center of conflicts. Clear it.
         self.symdiag = newsym
         self.store[ self.rowdiag, self.coldiag ] = self.empty
         self.colx[ self.coldiag, newsym ] = self.empty
         self.rowx[ self.rowdiag, newsym ] = self.empty

         # The start cell gets a new assignment unambigouously
         self.store[ self.rowstart, self.colstart] = self.symcorn
         self.colx[ self.colstart, self.symstart ] = self.empty
         self.rowx[ self.rowstart, self.symstart ] = self.empty
         self.colx[ self.colstart, self.symcorn ]  = self.rowstart
         self.rowx[ self.rowstart, self.symcorn ]  = self.colstart

         # Preserve conflict information for empty diagonal cell
         self.Moves += 1
         self.proper = False
      
   # End of proper move 



   # Set up the next operating rectangle at the conflicted cell carried
   # forward from the previous move.  Apply an "improper move" and check the
   # resulting state. 
   #
   def  JMimproperMove( self ) :

      #print("Setting up improper move")
      #print("Store\n",self.store)
      #print("STATE:  symstart ",self.symstart,"  symcorn ", self.symcorn,  
      #   "symdiag ", self.symdiag )
      #print("LOCATION: rowstart ",self.rowstart,"  rowdiag ",self.rowdiag,
      #   "  colstart ", self.colstart, "  coldiag ",self.coldiag)
      #print("rowx \n",self.rowx)
      #print("colx\n",self.colx)
      #quit()

      # Locate cells with conflicting labels, using index arrays.
      rowconf1 = self.rowstart
      rowconf2 = self.colx[self.coldiag,self.symstart]
      colconf1 = self.colstart
      colconf2 = self.rowx[self.rowdiag, self.symstart]

      #print("Candidates for next row: ", rowconf1, rowconf2)
      #print("Candidates for next col: ", colconf1, colconf2)
      
      # Randomly pick corner row and column such that at least one of them is new.
      option = randix(3)
      #print("Option ", option,"  selected, but use 2.")
      option = 2
      if  option==2 :
         # Both the new row and the new column were selected
         newrow = rowconf2 
         newcol = colconf2
         oldrow = rowconf1
         oldcol = colconf1
         symnew = self.store[newrow,newcol]
         #print(" Using new row ", newrow," and new column ", newcol)
         #print(" Restoring row ", oldrow," and old column ", oldcol)

         # Clear the new locations and restore the old
         self.store[self.rowdiag,newcol]  = self.empty
         self.store[newrow, self.coldiag] = self.empty
         self.rowx[newrow,self.symstart]  = self.empty
         self.colx[newcol,self.symstart]  = self.empty
         self.rowx[self.rowdiag,self.symstart]  = self.empty
         self.colx[self.coldiag,self.symstart]  = self.empty
         #print( "After clearing conflicts")
         #print("Store\n",self.store)

         # Restore the previous corner cells 
         self.store[self.rowstart,self.coldiag]  = self.symstart
         self.store[self.rowdiag, self.colstart] = self.symstart
         self.rowx[self.rowstart,self.symstart]   = self.coldiag
         self.colx[self.coldiag,self.symstart]   = self.rowstart
         self.rowx[self.rowdiag,self.symstart]   = self.colstart
         self.colx[self.colstart,self.symstart]  = self.rowdiag
         #print( "After restoring unused locations")
         #print("Store\n",self.store)
         #print("rowx\n",self.rowx)
         #print("colx\n",self.colx)
                   
      if  option==1 :
         # Use the new row and the previous column
         newrow = rowconf2
         newcol = colconf1
         oldrow = rowconf1
         oldcol = colconf2
         symnew = self.store[newrow,newcol]

         # Clear the new row location in the old start column
         self.store[newrow,self.colstart]  = self.empty
         self.rowx[newrow,self.symstart]   = self.empty
         self.colx[self.colstart,self.symstart]  = self.empty

         # Restore the old corner cell in the column
         self.store[oldrow,self.colstart]  = self.symstart
         self.rowx[oldrow,self.symstart]  = self.colstart
         self.colx[self.colstart,self.symstart]   = oldrow
   
      if  option==0 :
         # Use the new column and previous row
         newrow = rowconf1 
         newcol = colconf2
         oldrow = rowconf2
         oldcol = colconf1
         symnew = self.store[newrow,newcol]

         # Clear the new row location in the old start column
         self.store[newrow,self.colstart]  = self.empty
         self.rowx[newrow,self.symstart]   = self.empty
         self.colx[self.colstart,self.symstart]  = self.empty

         # Restore the old corner cell in the column
         self.store[oldrow,self.colstart]  = self.symstart
         self.rowx[oldrow,self.symstart]  = self.colstart
         self.colx[self.colstart,self.symstart]   = oldrow


      # Check for the first special case that returns to proper state
      if symnew == self.symdiag :
         #print("Case 3 move")

         # Case 3 applies.  This unambiguously assigns symbols to all cells 
         # in the rectangle, and returns to proper state.
         self.store[self.rowdiag, self.coldiag] = self.symcorn
         self.store[self.rowdiag, newcol] = self.symdiag
         self.store[newrow, self.coldiag] = self.symdiag        
         self.store[newrow, newcol] = self.symstart

         self.rowx[self.rowdiag,self.symcorn] = self.coldiag
         self.colx[self.coldiag,self.symcorn] = self.rowdiag
         self.rowx[self.rowdiag,self.symdiag] = newcol
         self.colx[newcol,self.symdiag] = self.rowdiag
         self.rowx[newrow,self.symdiag] = self.coldiag
         self.colx[self.coldiag,self.symdiag] = newrow
         self.rowx[newrow,self.symstart] = newcol
         self.colx[newcol,self.symstart] = newrow

         # The state becomes proper.
         self.Moves += 1
         self.clearstate( )
         #print("Store\n",self.store)
         #print("rowx \n",self.rowx)
         #print("colx\n",self.colx)

      # Check for the second special case that returns to proper state
      elif symnew == self.symcorn :
         #print("Case 4 move")

         # Case 4 applies.  Perform the move operation.  This unambiguously
         # values all cells in the rectangle.
         self.store[self.rowdiag, self.coldiag] = self.symdiag
         self.store[self.rowdiag, newcol] = self.symcorn
         self.store[newrow, self.coldiag] = self.symcorn
         self.store[newrow, newcol] = self.symstart

         self.rowx[self.rowdiag,self.symdiag] = self.coldiag
         self.colx[self.coldiag,self.symdiag] = self.rowdiag
         self.rowx[self.rowdiag,self.symcorn] = newcol
         self.colx[newcol,self.symcorn] = self.rowdiag
         self.rowx[newrow,self.symcorn] = self.coldiag
         self.colx[self.coldiag,self.symcorn] = newrow
         self.rowx[newrow,self.symstart] = newcol
         self.colx[newcol,self.symstart] = newrow

         # The state becomes proper.
         self.Moves += 1
         self.clearstate( )
         #print("Store\n",self.store)
         #print("rowx \n",self.rowx)
         #print("colx\n",self.colx)


      # For any other cases, conflicts change but state remains improper.
      # Conflicts are resolved at the diagonal cell, but remain present
      # at the new diagonal and corner cells.
      else :
         # Case 5 applies.  
         # Pick a symbol for resolving the old diagonal cell conflict.
         if  choice() :
            #print("Choice 1, case 5 resolvs with old diagonal")

            # Resolve with previous diagonal symbol
            self.store[self.rowdiag, self.coldiag] = self.symdiag
            self.rowx[self.rowdiag,self.symdiag] = self.coldiag
            self.colx[self.coldiag,self.symdiag] = self.rowdiag

            # Update the conflict state information
            self.symdiag = symnew
            swap = self.symcorn
            self.symcorn = self.symstart
            self.symstart = swap

         else :
            #print("Choice 0, case 5 resolvs with old corner")

            # Resolve with previous corner symbol
            self.store[self.rowdiag, self.coldiag] = self.symcorn
            self.rowx[self.rowdiag,self.symcorn] = self.coldiag
            self.colx[self.coldiag,self.symcorn] = self.rowdiag

            # Update the conflict state information
            self.store[newrow,newcol] = self.empty
            self.symcorn = self.symstart
            self.symstart = self.symdiag
            self.symdiag = symnew

         # The corner cells remain empty.  Clear the new
         # conflict cell.
         self.store[newrow,newcol] = self.empty
         self.rowx[newrow,symnew]  = self.empty
         self.colx[newcol,symnew]  = self.empty
             
         # Move to new location
         self.rowstart = self.rowdiag
         self.colstart = self.coldiag
         self.rowdiag = newrow
         self.coldiag = newcol
         #print("After case 5 move:")
         #print("  Store\n",self.store)
         #print("  rowx \n",self.rowx)
         #print("  colx\n",self.colx)


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
