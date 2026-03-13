import h5py
import cooler
from CoolProcessor import CoolProcessor

class MulticoolProcessor:
    """
    A class to process multi-cell .scool files, specifically for 2GB .mcool files 
    like 'nagano_10kb_raw.scool'. It provides functionality to list cells, 
    extract and process individual cell data, and save specific cells to new .cool files.
    
    Attributes:
    -----------
    fn : str
        File path to the .scool file.

    cellNames : list
        List of all cell names in the .scool file, initialized in the constructor.
    """

    def __init__(self, fn="sourceData/nagano_10kb_cell_types.scool"):
        """
        Initializes the MulticoolProcessor by loading the .scool file and extracting 
        the names of all cells in the file.
        
        Parameters:
        -----------
        fn : str
            File path to the .scool file.
        """
        self.fn = fn
        with h5py.File(fn, 'r') as f:
            # List all groups that correspond to cells
            self.cellNames = [name for name in f['cells']]
        
    def getCellNames(self):
        """
        Returns the list of cell names in the .scool file.
        
        Returns:
        --------
        list
            A list containing all the cell names from the .scool file.
        """
        return self.cellNames
    
    def readCell(self, cellName):
        """
        Creates a CoolProcessor object for a specific cell.
        
        Parameters:
        -----------
        cellName : str
            The name of the cell to be processed.
        
        Returns:
        --------
        CoolProcessor
            A CoolProcessor object initialized with the cell's data.
        """
        cool_uri = f"{self.fn}::/cells/{cellName}"
        return CoolProcessor(cool_uri, cellName)

    def getCoolObject(self, cellName="hicBuildMatrix_MATRIX_on_data_16179_and_data_16116"):
        """
        Creates and returns a Cooler object for a specific cell.
        
        Parameters:
        -----------
        cellName : str, optional
            The name of the cell to be loaded as a Cooler object. Default is a specific cell.
        
        Returns:
        --------
        cooler.Cooler
            A Cooler object representing the cell's Hi-C matrix data.
        """
        cool_uri = f"{self.fn}::/cells/{cellName}"
        c = cooler.Cooler(cool_uri)
        return c
    
    def saveCoolToSeperateFile(self, cellName, resultFn):
        """
        Saves the Hi-C matrix data of a specific cell to a separate .cool file.
        
        Parameters:
        -----------
        cellName : str
            The name of the cell whose data is to be saved.
        
        resultFn : str
            The path where the new .cool file should be saved, including the filename.
        
        Returns:
        --------
        None
        """
        cool_uri = f"{self.fn}::/{cellName}"
        cooler.fileops.cp(cool_uri, resultFn)
