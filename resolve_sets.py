r"""
For Sequencing the promoter library, a pooled amplicon sequencing approach is chosen. 
4 Plates each will be combined into one 'set' of a 2x2 arrangement (24x16 wells total).
From each set, an aliquot of all horizontal rows and all vertical columns will be pooled
(h-pools and v-pools) and PCR amplified. Said amplicon will be sequenced using Illumnia
Nextera Amplicon Seq. 

EXAMPLE OF ONE SET CONSISTING OF FOUR 96 WELL PLATES:
    _______  _______  |
    |  1  |  |  2  |  |
    |_____|  |_____|  |
    _______  _______  |   v-pool sampling 
    |  3  |  |  4  |  |     => corresponds to columns
    |_____|  |_____| \|/
   ------------------>
      h-pool sampling
        => corresponds to rows

Each pool PCR reaction will be performed with a specific and unique barcode pair on each 
primer. Therefore, nested pools can be demultiplexed after the sequencing allowing the 
reconstruction of each  pools and its sequenced promoters. To identify single wells/promoters,
one searches for barcodes / barcode pairs in both an h- and v-pool thus giving an x and a y
coordinate allowing to further demultiplex the reads. The second demultiplexing step is the 
the main goal of the script.
"""

import random
import itertools

__author__ = "Christian Rauch"
__version__ = "0.2a"

FILL = 5

class PlateSet:
    def __init__(self, rows=16, cols=24):
        if (isinstance(rows, float) and isinstance(cols, float)):
            if not (rows == int(rows) and cols == int(cols)):
                raise TypeError(
                    "Tried to initiliaze PlateSeq objects with non float-like data.")
            else:
                rows = int(rows)
                cols = int(cols)
        elif not (isinstance(rows, int) or isinstance(cols, int)):
            # rows and/or cols is not a float and not an integer
            raise TypeError("Tried to initiliaze PlateSeq objects with non-number data.")
        self.num_rows = rows
        self.num_cols = cols
        self.num_wells = self.num_cols * self.num_rows
        self._null = "x" * FILL
        self.grid = [[self._null for _ in range(self.num_cols)] for _ in range(self.num_rows)]
        

    def __repr__(self):
        as_str = ""
        for row in self.grid:
            as_str += str(row) + "\n"
        return as_str


    def value_of(self, *position: "row and col as int"):
        """Get the value of grid position."""
        try:
           row, col = position
        except ValueError:
            row, col = position[0]
        if not (row in range(self.num_rows) and col in range(self.num_cols)):
            raise ValueError("Trying to get value of well outside of grid.")
        return self.grid[row][col]


    def get_rows(self):
        rows = []
        for row in self.grid:
            rows.append(row)
        return rows


    def get_cols(self):
        cols = {i: [] for i in range(self.num_cols)}
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                cols[col_idx].append(self.value_of(row_idx, col_idx))
        cols_as_list = [cols[i] for i in range(self.num_cols)]
        return cols_as_list


    def insert(self, value, *position):
        """Set the value of grid position to desired value."""
        try:
           row, col = position
        except ValueError:
            row, col = position[0]
        if not (row in range(self.num_rows) and col in range(self.num_cols)):
            raise ValueError("Trying to insert a value in a well outside of grid.")
        self.grid[row][col] = value
        return None


    def remove(self, *position):
        """Reset the value of a grid position to default."""
        try:
           row, col = position
        except ValueError:
            row, col = position[0]
        if not (row in range(self.num_rows) and col in range(self.num_cols)):
            raise ValueError("Trying to remove value from a well outside of grid.")
        self.grid[row][col] = self._null
        return None


    def count_defaults(self):
        """Counts the number of default values in the PlateSet.
        Used to compute the number of unsolved wells after the solve() function call."""
        num_default = 0
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                if self.value_of(row_idx, col_idx) == self._null:
                    num_default += 1
        return num_default


    def clear(self):
        """Clears the complete grid and resets its values."""
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                self.remove(row_idx, col_idx)
        return None


    def fill_randomly(self, seed=False):
        """Randomly fills the complete grid with numbers of the alphabet. 
        Randomness can be seeded.
        WARNING: Will override any preexisting content!""" 
        if seed:
            random.seed(seed)
        elements = list(itertools.combinations("ABCDEFGHIJKLMNOPQRSTUVWXYZ", FILL))
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                value = "".join(random.choice(elements))
                self.insert(value, row_idx, col_idx)
        return None


    def fill_by_other(self, other):
        """Fills all wells with the well value of another grid."""
        if self.num_rows != other.num_rows and self.num_cols != other.num_rows:
            raise ValueError(f"Cannot fill {self.num_rows}x{self.num_cols} " \
                + f"plate with {other.num_rows}x{other.num_cols} plate (incompatible dimensions).")

        for row_idx in range(other.num_rows):
            for col_idx in range(other.num_cols):
                self.insert(other.value_of(row_idx, col_idx), row_idx, col_idx)
        return None


    def copy(self):
        """Returns a copy of itself.
        Deep copying algorithm."""
        copy = PlateSet(int(self.num_rows), int(self.num_cols))
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                copy.insert(self.value_of(row_idx, col_idx), row_idx, col_idx)
        return copy


    def is_equal_to(self, other):
        """Returns True if current grid is equal the another grid.
        Returns False otherwise."""
        if not self.num_rows == other.num_rows and self.num_cols == other.num_cols:
            return False
        for row_idx in range(self.num_rows):
            for col_idx in range(self.num_cols):
                if not self.value_of(row_idx, col_idx) == other.value_of(row_idx, col_idx):
                    return False
        return True


    def hpools(self):
        """Retunrs a dictionary containing row numbers as keys and sets of elements representing hpools."""
        hpools = {i: set() for i in range(self.num_rows)}
        for row_idx, row in enumerate(self.get_rows()):
            for value in row:
                hpools[row_idx].add(value)
        return hpools
    

    def vpools(self):
        """Retunrs a dictionary containing col numbers as keys and sets of elements representing vpools."""
        vpools = {i: set() for i in range(self.num_cols)}
        for col_idx, col in enumerate(self.get_cols()):
            for value in col:
                vpools[col_idx].add(value)
        return vpools


    def _indices_to_well(self, *position, start_at=1):
        """Transforms and row_idx, col_idx pair to a plate number with well ID.
        
        WORKS ONLY WITH PLATESETS SMALLER OR EQUAL SIZED AS DEFAULT."""
        try:
           row, col = position
        except ValueError:
            row, col = position[0]
        if row > 16 or col > 24:
            raise ValueError(
                f"Cannot transform index to well number for {self.num_rows}x{self.num_cols} PlateSet"
            )
        plate_row_offset = row // 8
        plate_col_offset = col // 12
        offset_to_num = {
            (0, 0): 0,
            (0, 1): 1,
            (1, 0): 2,
            (1, 1): 3
        } 
        plate_number = start_at + offset_to_num[(plate_row_offset, plate_col_offset)]

        row = row % 8
        col = col % 12
        well = "ABCDEFGH"[row] + str(col + 1)

        return f"Plate {plate_number} {well}"


    def export(self, filename, type="csv", start_at=1):
        """Exports PlateSeq instance into a file for later use.

        filename:      Name of the later file. If no extension is given, type is used for that.
        type:          File format of the output file.
        start_as:      Output file will contain human readible plate number/well ID pairs, 
                       starting with plate number specified as start_at.
        """ 
        if not "." in filename: # No '.' in filename indicates missing file extension
            filename += "." + type
        if type.lower() == "csv":
            to_print = ""
            for row_idx in range(self.num_rows):
                for col_idx in range(self.num_cols):
                    to_print += ";".join([self._indices_to_well(row_idx, col_idx, start_at=start_at),
                                          self.value_of(row_idx, col_idx)])
                    to_print += "\n"
            with open(filename, "w") as outfile:
                outfile.write(to_print)
        else:
            raise NotImplementedError
        

def solve(
    hpools: "dict = {0: {...}, 1: {...}",
    vpools: "dict = {0: {...}, 1: {...}",
    num_rows=16,
    num_cols=24
    ) -> "Solved PlateSet":
    """Takes hpools and vpools as input and returns a solved PlateSet."""
    solved_set = PlateSet(num_rows, num_cols)
    solved_items = {}
    for hpool_idx, hpool_set in hpools.items():
        for item in hpool_set:
            for vpool_idx, vpool_set in vpools.items():
                if item in vpool_set:
                    if item in solved_items:
                        solved_items[item].append((hpool_idx, vpool_idx))
                    else:
                        solved_items[item] = [(hpool_idx, vpool_idx)]
    for item in solved_items:
        if len(solved_items[item]) == 1:
            solved_set.insert(item, solved_items[item][0])
    return solved_set


def fasta_to_pool():
    pass


if __name__ == "__main__":
    number_of_tests = 10
    completely_solved = 0
    unsolved = []

    for _ in range(number_of_tests):
        randomset = PlateSet()
        randomset.fill_randomly()
        solvedset = solve(randomset.hpools(), randomset.vpools())
        if randomset.is_equal_to(solvedset):
            completely_solved += 1
        else:
            unsolved.append(solvedset.count_defaults())


    from statistics import mean
    print(f"{completely_solved} out of {number_of_tests} sets could be solved completely. Else {mean(unsolved)} could not be resolved.")

    # testset = PlateSet()
    # testset.fill_randomly(seed=4)

    from pprint import pprint
    solvedset.export("test1.csv", "csv", 6)
    pprint(solvedset)
    # solvedset.export("test2", "CSV")
    # solvedset.export("test3", "xlsx")

    # solved = solve(testset.hpools(), testset.vpools())
    # print(testset)
    # print()
    # print(solved)
    # print()
    # print(testset.is_equal_to(solved))

    # print(testset)
    # pprint(testset.hpools())
    # pprint(testset.vpools())
