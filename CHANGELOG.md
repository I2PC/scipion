
## v1.0.1 (2016-06-30) 

- Several protocol fixes:
   - Fixed bug when creating the output for Frealign (in some cases some information from input particles was not properly propagated)
   - Fixed some bugs in movie alignment protocols (summovie and unblur) and tests added
   - Some minor bugs fixed in Relion protocols
   - Bugs fixed in Resmap protocol when using two half volumes
- Fixed several bugs in Spider protocols: 
   - converting input particles with alignment 
   - wrong regular expression for replacing some variables in script template
   - parsing of the resulting dendrogram
   - some additional validations and removed unused code
- Bugfixes and inprovements in Xmipp protocols:
   - Protocols screen-classes merged into one: compare-reprojections
   - Complete refactoring of operate-particles and operate-volumes protocols (previously called 2D and 3D calculator). Tests added
- Picking and Viewer: 
   - Warning if particles are picked in a temporary folder and the SetOfParticles was not created
   - Improved implementation of assign-tiltpairs protocol in Xmipp and some refactoring of picking methods
   - Fixed bug that caused GUI to freeze sometimes
   - Some bugs fixed when displaying and exporting particles
   - Sorting arrows displayed after sorting by a column. Hourglass displayed while sorting.
   - Some bug fixed when creating subset from classes
- Other fixes or improvements:
   - ImageHandler's methods convert and writeStack now accepts alignment parameters
   - Fixed bug when displaying Movies summary (sqlite files were not closed)
   - Fixed bug when spawning Eman process to write particles
   - Added REMOTE_MESA_LIB environment var for using OpenGL in remote desktops
   - Created a LegacyProtocol class to read deprecated protocols
   - Cleanup in some tests and added new ones for core classes or functions

## v1.0.0 (2016-02-20) Augusto

- Allows to combine several EM software packages (~ 100 protocols):
  - All protocols from Xmipp 
  - Most of protocols from Relion
  - MDA protocols from Spider
  - Some protocols from Eman2/Sparx
  - From Grigorieff lab: CTFFIND, FREALIGN, unblur and summovie.
  - A few tools from Bsoft
  - ResMap, gEMpicker, dogpicker, motioncorr
- Full tracking and reproducibility:
  - Display runs as a list or a tree.
  - Inspect the parameters of a previous run
  - Repeat one or several runs
  - Export/Import a workflow template
- Data analysis:
  - Visualization and operation with Sets. (Particles, Micrographs, CTFs, etc)
  - Visualization of Volumes
  - Resolution and angular distribution plots
