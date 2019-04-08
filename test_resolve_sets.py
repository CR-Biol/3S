import random
import unittest

import resolve_sets

# Testing PlateSeq class
class TestPlateSeq(unittest.TestCase):
    def setUp(self):
        self.plate_small = resolve_sets.PlateSet(2, 3)
        self.plate_default = resolve_sets.PlateSet()
        self.plate_other_default = resolve_sets.PlateSet()
        self.null = resolve_sets.FILL * "x"


    def test_init(self):
        self.assertEqual(self.plate_small.num_rows, 2)
        self.assertEqual(self.plate_small.num_cols, 3)
        self.assertEqual(self.plate_small.num_wells, 6)
        self.assertEqual(self.plate_small._null, self.null)
        self.assertEqual(
            self.plate_small.grid,
            [[self.null, self.null, self.null],
             [self.null, self.null, self.null]]
        )

        self.assertEqual(self.plate_default.num_rows, 16)
        self.assertEqual(self.plate_default.num_cols, 24)
        self.assertEqual(self.plate_default.num_wells, 384)

        with self.assertRaises(TypeError):
            resolve_sets.PlateSet("a", "b")

        with self.assertRaises(TypeError):
            resolve_sets.PlateSet(1.2, 3.2)

        try:
            resolve_sets.PlateSet(1.0, 2.0)
        except:
            self.fail("Failed to initialize PlateSet with integer-like floats.")


    def test_value_of(self):
        self.assertEqual(self.plate_small.value_of(1, 2), self.null)
        self.assertEqual(self.plate_small.value_of((1, 2)), self.null)
        
        with self.assertRaises(ValueError):
            self.plate_small.value_of(-1, 1)

        with self.assertRaises(ValueError):
            self.plate_small.value_of(100, 1)


    def test_insert_and_remove(self):
        self.plate_small.insert("TEST", (1, 1))
        # Check if ONLY position (1, 1) has been changed
        for row in range(self.plate_small.num_rows):
            for col in range(self.plate_small.num_cols):
                if (row, col) == (1, 1):
                    self.assertEqual(self.plate_small.value_of(row, col), "TEST")
                else:
                    self.assertEqual(self.plate_small.value_of(row, col), self.null)
        # Undo the change and check if the value has been reset correctly.
        self.plate_small.remove((1, 1))
        self.assertEqual(self.plate_small.value_of(1, 1), self.null)
    

    def test_get_rows(self):
        # Fill small plate with arbitrary values.
        self.plate_small.insert("A", (0, 0))
        self.plate_small.insert("B", (0, 1))
        self.plate_small.insert("C", (0, 2))
        self.plate_small.insert("D", (1, 0))
        self.plate_small.insert("E", (1, 1))
        self.plate_small.insert("F", (1, 2))
        
        self.assertEqual(
            self.plate_small.get_rows(),
            [["A", "B", "C"], ["D", "E", "F"]]
        )
        # Clear out small plate.
        self.plate_small.clear()
        

    def test_get_cols(self):
        # Fill small plate with arbitrary values.
        self.plate_small.insert("A", (0, 0))
        self.plate_small.insert("B", (0, 1))
        self.plate_small.insert("C", (0, 2))
        self.plate_small.insert("D", (1, 0))
        self.plate_small.insert("E", (1, 1))
        self.plate_small.insert("F", (1, 2))

        self.assertEqual(
            self.plate_small.get_cols(),
            [["A", "D"], ["B", "E"], ["C", "F"]]
        )
        # Clear out small plate.
        self.plate_small.clear()


    def test_count_defaults(self):
        self.assertEqual(
            self.plate_small.count_defaults(), 
            self.plate_small.num_wells
        )

        self.assertEqual(
            self.plate_default.count_defaults(),
            self.plate_default.num_wells
        )

        self.plate_small.insert("TEST", 1, 1)
        self.assertEqual(
            self.plate_small.count_defaults(),
            self.plate_small.num_wells - 1
        )
        self.plate_small.clear()


    def test_fill_randomly_and_clear(self):
        random.seed(1)
        self.plate_default.fill_randomly()
        # Check if plate is not empty anymore.
        for row in range(self.plate_default.num_rows):
            for col in range(self.plate_default.num_cols):
                self.assertNotEqual(self.plate_default.value_of(row, col), self.null)

        self.plate_default.clear()
        # Check if plate is empty again.
        for row in range(self.plate_default.num_rows):
            for col in range(self.plate_default.num_cols):
                self.assertEqual(self.plate_default.value_of(row, col), self.null)


    def test_fill_by_other(self):
        random.seed(1)
        # Fill one plate randomly and fill the other one with its values.
        self.plate_default.fill_randomly()
        self.plate_other_default.fill_by_other(self.plate_default)
        for row in range(self.plate_default.num_rows):
            for col in range(self.plate_default.num_cols):
                self.assertEqual(
                    self.plate_default.value_of(row, col),
                    self.plate_other_default.value_of(row, col)
                ) 
        # Trying to combine plates with incompatible dimensions should raise ValueError.
        with self.assertRaises(ValueError):
            self.plate_small.fill_by_other(self.plate_default)
        
        # Clear all plateas
        self.plate_default.clear()
        self.plate_other_default.clear()
        self.plate_small.clear()


    def test_copy(self):
        random.seed(1)
        self.plate_default.fill_randomly()
        copied_plate = self.plate_default.copy()
        self.assertTrue(self.plate_default.is_equal_to(copied_plate))
        copied_plate.remove(0, 0)
        self.assertFalse(self.plate_default.is_equal_to(copied_plate))
        self.plate_default.clear()


    def test_hpools(self):
        # Fill small plate with arbitrary values.
        self.plate_small.insert("A", (0, 0))
        self.plate_small.insert("B", (0, 1))
        self.plate_small.insert("C", (0, 2))
        self.plate_small.insert("D", (1, 0))
        self.plate_small.insert("DUPLICATE", (1, 1))
        self.plate_small.insert("DUPLICATE", (1, 2))

        self.assertEqual(
            self.plate_small.hpools(),
            {0: {"A", "B", "C"}, 1: {"D", "DUPLICATE"}}
        )
        # Clear plate.
        self.plate_small.clear()


    def test_vpools(self):
        # Fill small plate with arbitrary values.
        self.plate_small.insert("A", (0, 0))
        self.plate_small.insert("B", (0, 1))
        self.plate_small.insert("DUPLICATE", (0, 2))
        self.plate_small.insert("D", (1, 0))
        self.plate_small.insert("E", (1, 1))
        self.plate_small.insert("DUPLICATE", (1, 2))

        self.assertEqual(
            self.plate_small.vpools(),
            {0: {"A", "D"}, 1: {"B", "E"}, 2: {"DUPLICATE"}}
        )
        # Clear plate.
        self.plate_small.clear()


    def test_indices_to_well(self):
        self.assertEqual(
            self.plate_default._indices_to_well(0, 0), 
            "Plate 1 A1")
        self.assertEqual(
            self.plate_default._indices_to_well(0, 0, start_at=2), 
            "Plate 2 A1")
        self.assertEqual(
            self.plate_default._indices_to_well(0, 11), 
            "Plate 1 A12")
        self.assertEqual(
            self.plate_default._indices_to_well(8, 0),
            "Plate 3 A1")
        self.assertEqual(
            self.plate_default._indices_to_well(0, 12),
            "Plate 2 A1")
        self.assertEqual(
            self.plate_default._indices_to_well(0, 12, start_at=3),
            "Plate 4 A1")


# Testing solve() function
class TestSolve(unittest.TestCase):
    def test_solve(self):
        pass


if __name__ == "__main__":
    unittest.main()
