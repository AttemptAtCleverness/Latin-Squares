import LatinSquares
import unittest

class TestLatinSquaresFunctions(unittest.TestCase):
    
    def test_generate_all_squares(self):
        # Run with low n to save execution time 
        three = LatinSquares.generate_all_squares(3)
        four = LatinSquares.generate_all_squares(4)
        
        self.assertEqual(len(three),12)
        self.assertEqual(len(four),576)

    def test_generate_all_squares_standardized(self):
        four = LatinSquares.generate_all_squares_standardized(4)
        five = LatinSquares.generate_all_squares_standardized(5)
        
        self.assertEqual(len(four),24)
        self.assertEqual(len(five),1344)
        
        return
    
    def test_generate_all_squares_reduced(self):
        four = LatinSquares.generate_all_squares_reduced(4)
        five = LatinSquares.generate_all_squares_reduced(5)
        six = LatinSquares.generate_all_squares_reduced(6)
        
        self.assertEqual(len(four),4)
        self.assertEqual(len(five),56)
        self.assertEqual(len(six),9408)
        return
    
    def test_standardize(self):
        square1 = [
            [3, 1, 2],
            [2, 3, 1],
            [1, 2, 3]
        ]
        square1_standardized = [
            [1, 2, 3],
            [3, 1, 2], 
            [2, 3, 1] 
        ] 
        
        square2 = [
            [4, 1, 3, 2],
            [2, 4, 1, 3],
            [1, 3, 2, 4],
            [3, 2, 4, 1]
        ]
        square2_standardized = [
            [1, 2, 3, 4],
            [4, 1, 2, 3],
            [2, 3, 4, 1],
            [3, 4, 1, 2] 
        ]
        square3 = [
            [1, 2, 3, 4, 5],
            [2, 3, 4, 5, 1],
            [3, 4, 5, 1, 2],
            [4, 5, 1, 2, 3],
            [5, 1, 2, 3, 4]
        ]

        self.assertEqual(LatinSquares.standardize(square1),square1_standardized)
        self.assertEqual(LatinSquares.standardize(square2),square2_standardized)
        self.assertEqual(LatinSquares.standardize(square3), square3)

    def test_reduce(self):
        return
    
    

if __name__ == '__main__':
    unittest.main()
    # running individual function
    ''' 
    suite = unittest.TestSuite()
    suite.addTest(TestLatinSquareFunctions('{FUNCTION NAME HERE}'))
    runner = unittest.TextTestRunner()
    runner.run(suite)
    '''
    