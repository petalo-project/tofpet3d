Reconstruction
********************

Here we consider a 3D PET MLEM reconstruction algorithm, the derivation of which
comes from `these slides, based on the IAEA Nuclear Medicine Physics handbook
<https://humanhealth.iaea.org/HHW/MedicalPhysics/e-learning/Nuclear_Medicine_Handbook_slides/Chapter_13._Image_Reconstruction.pdf>`_.

We consider list-mode reconstruction, meaning that the input will be the
measured lines of response (LORs), each defined by two interaction points,
:math:`(x_1,y_1,z_t,t_1)` and :math:`(x_2,y_2,z_2,t_2)`, and we seek to
reconstruct a 3D image in a volume divided into :math:`N` 3D voxels. To derive the MLEM
procedure, we will consider several quantities. Note that in the following
definitions, index :math:`i` will run over the lines of response and index
:math:`j` will run over the voxels of the reconstructed image:

- :math:`\lambda_j`: the value of the reconstructed image at voxel :math:`j`
- :math:`A_{ij}`: the probability that the emission of electron-positron annihilation
  radiation within voxel :math:`j` will lead to detection in LOR :math:`i` (referred to as the "system matrix")
- :math:`y_i`: the measured number of counts in LOR :math:`i`

In practice, since LORs are not binned in any way (we consider a list of all
the LORs observed), :math:`y_i = 1` for each measured LOR and all other
LORs (theoretically an infinite number) will have :math:`y_i = 0`.

We now consider a theoretical construction :math:`x_{ij}`, which contains the
number of measured counts in LOR :math:`i` that were emitted from voxel :math:`j`.
We model these as Poisson distributed about a mean value :math:`\bar{x}_{ij} = A_{ij}\lambda_j`.
The `log-likelihood for observation <https://en.wikipedia.org/wiki/Poisson_distribution>`_
of the image with voxels :math:`\lambda_j` given the measured :math:`x_{ij}` is

.. math::
  L(\lambda) = \sum_{ij}x_{ij}\ln(A_{ij}\lambda_j) - A_{ij}\lambda_j.

To maximize the likelihood, we look for :math:`\lambda_{j'}` such that :math:`\partial L/\partial\lambda_{j'} = 0`.

.. math::
  \frac{\partial L}{\partial\lambda_{j'}} = \sum_{ij}\Bigl(\frac{x_{ij}}{\lambda_j}\delta_{jj'} - A_{ij}\delta_{jj'}\Bigr)
  = \frac{\sum_{i}x_{ij'}}{\lambda_j'} - \sum_{i}A_{ij'} = 0.

Therefore we require that for each :math:`\lambda_j`

.. math::
  \lambda_j = \frac{\sum_{i}x_{ij}}{\sum_{i}A_{ij}}.

This will be solved iteratively, denoting the voxel values at step :math:`k` as
:math:`\lambda^{(k)}_j` starting with an initial guess :math:`\lambda^{(0)}_j`
and solving for each successive :math:`\lambda^{(k+1)}_j` given the measurements
:math:`x^{(k)}_{ij}` at step :math:`k`,

.. math::
  \lambda^{(k+1)}_j = \frac{\sum_{i}x^{(k)}_{ij}}{\sum_{i}A_{ij}}.

Since we don't have all the information in :math:`x_{ij}`, we must now express it
in terms of actual measurements that have been made. If we just use the expected
value :math:`\langle x^{(k)}_{ij}\rangle = A_{ij}\lambda^{(k)}_j`, we will get
:math:`\lambda^{(k+1)}_j = \lambda^{(k)}_j`. To actually incorporate the
measured LORs, we write the measured number of counts in each LOR :math:`y_{i}`
as

.. math::
  y_{i} = \sum_{j}A_{ij}\lambda^{(k)}_{j} + b_{i},

where the term :math:`b_{i}` contains additional contributions to the LOR, such
as counts due to random scatters. Then we divide :math:`y_i` by the entire right
side of the equation above, writing

.. math::
  \frac{y_{i}}{\sum_{j}A_{ij}\lambda^{(k)}_{j} + b_{i}} = 1,

and we multiply :math:`\langle x^{(k)}_{ij}\rangle` by 1,

.. math::
  \langle x^{(k)}_{ij}\rangle \rightarrow \langle x^{(k)}_{ij}\rangle\times 1 = A_{ij}\lambda^{(k)}_j \times \frac{y_{i}}{\sum_{j}A_{ij}\lambda^{(k)}_{j} + b_{i}}.

Now we can write the iterative equation for :math:`\lambda^{(k+1)}_j` with the
above expression for the expected value :math:`\langle x^{(k)}_{ij}\rangle`,

.. math::
  \lambda^{(k+1)}_j = \frac{\lambda^{(k)}_j}{\sum_{i}A_{ij}}\times \sum_{i}\Biggl(\frac{y_{i}A_{ij}}{\sum_{j}A_{ij}\lambda^{(k)}_{j} + b_{i}}\Biggr)

Since we will consider each LOR individually, we have for measured lines of response
:math:`i_{m}`, :math:`y_{i_m} = 1` and for all others :math:`y_{i} = 0`, so we can
rewrite the sum over :math:`i` as the sum over measured LORs :math:`i_m`,

.. math::
  \lambda^{(k+1)}_j = \frac{\lambda^{(k)}_j}{\sum_{i}A_{ij}}\times \sum_{i_{m}}\Biggl(\frac{A_{i_{m}j}}{\sum_{j}A_{i_{m}j}\lambda^{(k)}_{j} + b_{i_{m}}}\Biggr).

**This is the final equation that must be implemented in the code.** Some notes
on the content of this equation:

- :math:`A_{i_{m}j}` will be calculated on-the-fly for each measured LOR :math:`i_{m}`. Note that this
  is done by considering the *radiological path* of the line of response through the active volume, that is,
  the set of active voxels encountered by the LOR from :math:`(x_1,y_1,z_t)` to :math:`(x_2,y_2,z_2)`. These
  voxels can be determined using the Siddon algorithm
  [`R. Siddon. Med. Phys. 12 (2), 252 (1985) <https://aapm.onlinelibrary.wiley.com/doi/abs/10.1118/1.595715>`_].
  Furthermore, the voxels can be weighted using the time of flight information,
  such that rather than considering the path as a uniform line in which the
  probability that the electron-positron emission came from any given voxel on
  the line is equal, the distribution can be a gaussian centered on the voxel
  corresponding to the recorded time-of-flight.
- :math:`P_{i_{m}} \equiv \sum_{j}A_{i_{m}j}\lambda^{(k)}_{j} + b_{i_{m}}` can then be
  calculated using the present image voxels :math:`\lambda^{(k)}_{j}` and a
  given model for the additive contributions :math:`b_{i_{m}}`. Because
  :math:`P_{i_{m}}` is a sum over contributions to LOR :math:`i_{m}` from all
  voxels, it is called the *projection* (or *forward projection*) onto :math:`i_{m}`.
- The sum :math:`BP_{j} = \sum_{i_{m}}(A_{i_{m}j}/P_{i_{m}})`
  distributes the measured LORs back onto the image according to the probabilities
  :math:`A_{i_{m}j}` that the observed LOR was due to emission from voxel :math:`j`.
  This is called the *back projection*.
- The sum :math:`S_{j} \equiv \sum_{i}A_{ij}` is the total probability of observing any LOR
  emitted from voxel :math:`j`. This is called the *sensitivity matrix*, as it
  describes how "sensitive" the detector is to emission from a specific location
  in the active region.

One iteration of the algorithm is then performed as follows.

1. For each LOR :math:`i_{m}`:

     a. Compute the projection :math:`P_{i_{m}}` of the current image.
     b. Compute the contribution to the backprojection :math:`A_{i_{m}j}/P_{i_{m}}`
        and add it to the backprojection :math:`BP_{j}` of voxel :math:`j`.

2. Compute the :math:`\lambda^{(k+1)}_{j}` by multiplying each value
   :math:`BP_{j}` by the current voxel divided by the corresponding value in
   the sensitivity matrix, :math:`\lambda^{(k)}_{j}/S_{j}`.

A helper class ``MLEMReconstructor`` has been created to call the reconstruction
algorithm, implemented in C, from Python.

.. automodule:: mlem.mlem_reconstruct
   :members:
