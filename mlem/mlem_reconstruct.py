from ctypes import *
import os
import struct
import numpy as np


class MLEMReconstructor:
    """
    Python wrapper class to perform MLEM reconstruction. The class
    provides a method to call the C-based reconstruction using the
    configuration provided by the class variables.

    The reconstructed image is stored in a 1D array of 32-bit (4-byte) floats
    with values corresponding to:

        | ``(x0,y0,z0), (x1,y0,z0) ... (xN,y0,z0)``
        | ``(x0,y1,z0), (x1,y1,z0) ... (xN,y1,z0)``
        | ``...``
        | ``(x0,yN,z0), (x1,yN,z0) ... (xN,yN,z0)``

        | ``(x0,y0,z1), (x1,y0,z1) ... (xN,y0,z1)``
        | ``(x0,y1,z1), (x1,y1,z1) ... (xN,y1,z1)``
        | ``...``
        | ``(x0,yN,z1), (x1,yN,z1) ... (xN,yN,z1)``

        ``...``

        | ``(x0,y0,zN), (x1,y0,zN) ... (xN,y0,zN)``
        | ``(x0,y1,zN), (x1,y1,zN) ... (xN,y1,zN)``
        | ``...``
        | ``(x0,yN,zN), (x1,yN,zN) ... (xN,yN,zN)``

    **NOTE:** the LOR points may need to be sorted for the reconstruction to
    be correct

    :param prefix: the filename prefix for all saved files
    :param niterations: the number of iterations to perform on reconstruction
    :param save_every: save every specified number of iterations
    :param TOF: boolean to enable TOF
    :param TOF_resolution: the TOF resolution in ps
    :param img_size_xy: the image size in the x and y dimensions (in mm)
    :param img_size_z: the image size in the z dimension (in mm)
    :param img_nvoxels_xy: the number of voxels in the x and y dimensions
    :param img_nvoxels_z: the number of voxels in the z dimension
    :param libpath: the path to the C++ reconstruction library
    """

    def __init__(self, prefix: str = "mlem", niterations: int = 1,
                 save_every: int = -1, TOF: bool = True,
                 TOF_resolution: float = 200.,
                 img_size_xy: float = 180.0, img_size_z: float = 180.0,
                 img_nvoxels_xy: int = 60, img_nvoxels_z: int = 60,
                 libpath: str = "lib/libMLEM.so"):

        # Set default values for key variables.
        self.prefix = prefix
        self.niterations = niterations
        self.save_every = save_every
        self.TOF = True
        self.TOF_resolution = TOF_resolution
        self.img_size_xy = img_size_xy
        self.img_size_z = img_size_z
        self.img_nvoxels_xy = img_nvoxels_xy
        self.img_nvoxels_z = img_nvoxels_z

        # Load the C library.
        self.lib = cdll.LoadLibrary(libpath)

    def read_image(self, niter: int) -> np.ndarray:
        """
        Reads a reconstructed image stored as "(prefix)(niter).raw".

        :param niter: the iteration number to be read
        """

        # Construct an empty image.
        xdim = self.img_nvoxels_xy
        ydim = self.img_nvoxels_xy
        zdim = self.img_nvoxels_z
        img_arr = np.zeros([xdim,ydim,zdim])

        # Open the file containing the reconstructed image.
        fpath = "{}{}.raw".format(self.prefix,niter)
        if(os.path.isfile(fpath)):

            fimg = open(fpath,'rb')
            fdata = fimg.read()
            print("Read {} bytes".format(len(fdata)))

            # Read the image of nvoxels 32-bit (4-byte) floats.
            nvoxels = xdim*ydim*zdim
            if(len(fdata) == nvoxels*4):

                # Unpack the floats into a 1D array.
                s_arr = struct.unpack_from('f'*nvoxels, fdata)

                # Fill the 3D image from the array.
                for ivox,w in enumerate(s_arr):
                    i = int(ivox % xdim)
                    j = int(ivox / xdim) % ydim
                    k = int(ivox / (xdim*ydim))
                    img_arr[i,j,k] = w
            else:
                print("ERROR: unexpected data length {}".format(len(fdata)))
        else:
            print("ERROR: file {} not found".format(fpath))

        return img_arr

    def reconstruct(self, lor_x1, lor_y1, lor_z1, lor_t1,
                    lor_x2, lor_y2, lor_z2, lor_t2) -> np.ndarray:
        """
        Performs reconstruction by calling a C-implemented function. The
        reconstruction is list-mode, so the inputs are the coordinates
        :math:`(x,y,z,t)` for each line of response (LOR). These coordinates
        are specified as separate lists.

        :param lor_x1: x-coordinates for the first point in the LOR
        :type lor_x1: list
        :param lor_y1: y-coordinates for the first point in the LOR
        :type lor_y1: list
        :param lor_z1: y-coordinates for the first point in the LOR
        :type lor_z1: list
        :param lor_t1: time coordinates for the first point in the LOR
        :type lor_t1: list
        :param lor_x2: x-coordinates for the second point in the LOR
        :type lor_x2: list
        :param lor_y2: y-coordinates for the second point in the LOR
        :type lor_y2: list
        :param lor_z2: y-coordinates for the second point in the LOR
        :type lor_z2: list
        :param lor_t2: time coordinates for the second point in the LOR
        :type lor_t2: list
        :returns: array of shape [`img_size_xy`, `img_size_xy`,`img_size_z`] containing the reconstructed image
        """

        # Ensure arrays are all the same size.
        ncoinc = len(lor_x1)
        if(len(lor_y1) != ncoinc or len(lor_z1) != ncoinc
           or len(lor_t1) != ncoinc or len(lor_x2) != ncoinc
           or len(lor_y2) != ncoinc or len(lor_z2) != ncoinc
           or len(lor_t2) != ncoinc):
            print("ERROR: all LOR arrays must contain the same # of values")

        # ---------------------------------------------------------
        # Call the C-function for reconstruction.
        # ---------------------------------------------------------

        # Create C-arrays
        lor_x1 = (c_float * ncoinc)(*lor_x1)
        lor_y1 = (c_float * ncoinc)(*lor_y1)
        lor_z1 = (c_float * ncoinc)(*lor_z1)
        lor_t1 = (c_float * ncoinc)(*lor_t1)
        lor_x2 = (c_float * ncoinc)(*lor_x2)
        lor_y2 = (c_float * ncoinc)(*lor_y2)
        lor_z2 = (c_float * ncoinc)(*lor_z2)
        lor_t2 = (c_float * ncoinc)(*lor_t2)
        prefix = c_char_p(self.prefix.encode('utf-8'))

        # Set the argument types and return type.
        self.lib.MLEM_TOF_Reco.argtypes = (c_int, c_bool, c_float,
         c_float, c_float, c_int, c_int, c_int,
         POINTER(c_float), POINTER(c_float), POINTER(c_float), POINTER(c_float),
         POINTER(c_float), POINTER(c_float), POINTER(c_float), POINTER(c_float),
         c_char_p, c_int)

        self.lib.MLEM_TOF_Reco.restype = POINTER(c_float)

        # Call the function.
        img = self.lib.MLEM_TOF_Reco(self.niterations, self.TOF,
                                self.TOF_resolution, self.img_size_xy,
                                self.img_size_z, self.img_nvoxels_xy,
                                self.img_nvoxels_z, ncoinc,
                                lor_x1, lor_y1, lor_z1, lor_t1,
                                lor_x2, lor_y2, lor_z2, lor_t2,
                                prefix, self.save_every)

        # --------------------------------------------------
        # Prepare the numpy image from the C float pointer.
        # --------------------------------------------------

        # Construct some quantities for ease of use later.
        xdim = ydim = self.img_nvoxels_xy
        zdim = self.img_nvoxels_z
        nvoxels = xdim*ydim*zdim

        # Create a 3D numpy array of the correct dimensions.
        img_arr = np.zeros([xdim,ydim,zdim])

        # Extract the information from the C array into a 1D numpy array.
        rimg = np.array([img[i] for i in range(nvoxels)])

        # Fill the final 3D image, extracting the voxel indices from the
        #  single index in the 1D array.
        for ivox,w in enumerate(rimg):
            i = int(ivox % xdim)
            j = int(ivox / xdim) % ydim
            k = int(ivox / (xdim*ydim))
            img_arr[i,j,k] = w

        return img_arr
