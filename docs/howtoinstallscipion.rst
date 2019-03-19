.. figure:: https://github.com/I2PC/scipion/wiki/images/scipion_logo.gif
   :alt: http://scipion.cnb.csic.es

.. _howtoinstallscipion:
Installing Scipion 2.0 from sources
===================================

Step1: Download
---------------

You can install Scipion anywhere, as long as you have write permissions.
If you want to share Scipion installation among different users you can
use ``/usr/local`` or a similar path. If it is only for you, you can use
your home directory as well. If you have a previous installation, we
recommend starting this installation fresh in a new folder.

From GitHub
~~~~~~~~~~~

Clone Scipion repository (install **git** if not present in your
system):

::

    git clone https://github.com/I2PC/scipion.git
    cd scipion 
    git checkout release-2.0.0

Git will create a ``scipion`` directory under your current path; you do
not need to create it manually.

Step2: Dependencies
-------------------

To install Scipion from source, some development libraries are required.
You can install them with (this example is for Debian/Ubuntu distros.

-  For ubuntu 16:

``sudo apt-get install gcc-5 g++-5 cmake openjdk-8-jdk libxft-dev libssl-dev libxext-dev libxml2-dev libreadline6 libquadmath0 libxslt1-dev libopenmpi-dev openmpi-bin  libxss-dev libgsl0-dev libx11-dev gfortran libfreetype6-dev scons libfftw3-dev libopencv-dev curl``

-  For ubuntu 18:

``sudo apt-get install gcc-5 g++-5 cmake openjdk-8-jdk libxft-dev libssl-dev libxext-dev libxml2-dev libreadline6-dev libquadmath0 libxslt1-dev libopenmpi-dev openmpi-bin  libxss-dev libgsl0-dev libx11-dev gfortran libfreetype6-dev scons libfftw3-dev libopencv-dev curl``

-  For debian 9.7:

``sudo apt-get install gcc-5 g++-5 cmake openjdk-8-jdk libxft-dev libssl-dev libxext-dev libxml2-dev libreadline7 libquadmath0 libxslt1-dev libopenmpi-dev openmpi-bin  libxss-dev libgsl-dev libx11-dev gfortran libfreetype6-dev scons libfftw3-dev libopencv-dev curl``

-  **Note 1** - gcc and g++

You might already have g++ and gcc installed, but it is important to
have version 5 for compatibility with certain plugins like Xmipp3. Apart
from installing them, you'll need to set them as defaults (e.g. updating
the links resulting of typing ``which gcc`` and ``which g++`` to point
to ``gcc-5`` and ``g++-5`` respectively, or using update alternatives as
suggested
`here <https://askubuntu.com/questions/1087150/install-gcc-5-on-ubuntu-18-04>`__
).

-  **Note 2** - CUDA recommendations:

The following table contains the current compatibility of the plugins
that use CUDA and have some restrictions with the CUDA versions. To
maximize compatibility with all plugins we recommend to install
**CUDA-8.0**, but it's up to you to choose a different version depending
on the plugins you'd like to use (depending on this you might not need
CUDA at all).

+-----------------+------------+------------+------------+------------+------------+------------+------------+------------+
| Cuda version -> |      5.0   |      6.0   |      6.5   |      7.0   |      7.5   |      8.0   |      9.1   |      9.2   |
+=================+============+============+============+============+============+============+============+============+
| motioncorr      |            |            |            |            |            | ☑          | ☑          | ☑          |
+-----------------+------------+------------+------------+------------+------------+------------+------------+------------+
| gctf            | ☑          | ☑          | ☑          | ☑          | ☑          | ☑          |            |            |
+-----------------+------------+------------+------------+------------+------------+------------+------------+------------+
| gautomatch      |            |            |            | ☑          | ☑          | ☑          |            |            |
+-----------------+------------+------------+------------+------------+------------+------------+------------+------------+

Step3: Configure and Install
----------------------------

**Configure**

After installing the `dependencies <#step2-dependencies>`__, you can
proceed to generate configuration files. If you had a previous Scipion
installation, it is a good idea to make a copy of your current config
files ``~/.config/scipion/scipion.conf`` and
``<your_scipion_home>/config/scipion.conf``. Now run:

::

    cd scipion
    ./scipion config

You will be asked to share **scipion usage only** data. Sharing `usage
data <https://github.com/I2PC/scipion/wiki/Collecting-Usage-Statistics-for-Scipion>`__
will help to make Scipion better.

This command will generate the configuration files (if not present) and
will try to automatically find configuration paths.

If everything is OK (all green in the output) you can proceed to the
next step. If there is a problem (red colored output), you will need to
edit ``config/scipion.conf`` file in your preferred text editor and run
``./scipion config`` again.

One known change for Ubuntu 18 are the MPI paths in
``<your_scipion_home>/config/scipion.conf``:

::

    MPI_LIBDIR = /usr/lib/x86_64-linux-gnu/openmpi/lib
    MPI_INCLUDE = /usr/lib/x86_64-linux-gnu/openmpi/include/

Read more about `editing the configuration
file <Scipion-Configuration>`__.

The file ``config/hosts.conf`` contains some properties of the execution
machine. This configuration file is particularly important for clusters
that use a Queue System. If you are installing Scipion on a cluster, you
probably will want to check `how to configure an execution
host <Host-Configuration>`__.

**Compile**

To compile and install Scipion, just run:

::

    ./scipion install -j 5

``-j 5`` tells the Scipion installer to use 5 processors (cores) for
compilation. You should adjust this value according to your system.

If you have problems compiling Scipion, see
`Troubleshooting <https://github.com/I2PC/scipion/wiki/Troubleshooting>`__
page.

Step 4: Installing Xmipp3 and other EM Plugins
----------------------------------------------

Scipion can use many EM plugins. It is almost **mandatory to install
``scipion-em-xmipp``** (i.e. Scipion will run without it but with very
limited functionality). \* **For developers**: Developers might want to
build xmipp from the latest development version, please head
`here <https://github.com/I2PC/xmipp/wiki/Migrating-branches-from-nonPluginized-Scipion-to-the-new-Scipion-Xmipp-structure#xmipp>`__
if this is your case. You might also want to check how to [[install
plugins from the command line\|Developer guide: Installing plugins from
the command line]]. \* **For users**: To list and install plugins
including Xmipp, you can use the Plugin manager as shown below
(recommended) or alternatively, use the [[command line tool\|Developer
guide: Installing plugins from the command line]] mentioned for
developers. \* Run Scipion ``./scipion`` Because we haven't installed
xmipp yet, you'll see a message saying something like this in the
terminal:

::

    ```
    Scipion v2.0 (2019-03-12) Diocletian (release-2.0.0-fixes 50b9908)

    >>>>> python  /home/yaiza/Desktop/scipion/pyworkflow/apps/pw_manager.py 

     >>> WARNING: Xmipp binaries not found. Ghost active.....BOOOOOO!
      > Please install Xmipp to get full functionality. 
    (Configuration->Plugins->scipion-em-xmipp in Scipion manager window)

    ```

    * Open Plugin Manager (it will take a minute to load)

    [[images/guis/scipion_config_menu.png|Scipion projects manager]]

    * Select Xmipp to install it by clicking on the empty checkbox on the left. 

    [[images/guis/plugin_manager_install_xmipp.png|install xmipp]]

    * Add the number of processors you'd like to use (the more, the merrier!). Then click on the install button on the operations tab

    [[images/guis/plugin_manager_install_xmipp_install_button.png|install xmipp]]

    * Now we can check the progress on the Output log tab (or go make some coffee, Xmipp installation will take a bit!). 
    You might have to refresh the logs by clicking on the refresh symbol on the right. Please note that messages might not appear in order if we are using more than 1 processor.

    [[images/guis/plugin_manager_xmipp_install_logs.png|install xmipp logs]]

    * When the operation gets a green check, it's done!

    [[images/guis/plugin_manager_xmipp_done.png|install xmipp logs]]

    **Note**: if xmipp installation fails, you might have to uninstall it with the plugin manager:

    [[images/guis/plugin_manager_xmipp_uninstall.png|uninstall xmipp]]

    And manually remove leftover elements:

    ```
    rm -rf software/em/xmipp*
    ```

    * Now when we close and re-launch Scipion, we should get no messages.

    ```` ./scipion

    Scipion v2.0 (2019-03-12) Diocletian (release-2.0.0-fixes 50b9908)

    >>>>> python  /home/yaiza/Desktop/scipion/pyworkflow/apps/pw_manager.py 
    ````

Please refer to the `Plugin manager guide <Plugin-Manager>`__ to get
more details about plugin installation options.

Step 5: Cleaning up (Optional)
------------------------------

After Scipion is installed and properly working (see how to run tests in
the next section) one could clean some temporary files to free some disk
space after installation.

Remove the files under ``software/tmp`` folder:

::

    rm -rf sofware/tmp/*

The downloaded .tgz files of the EM packages can also be removed:

::

    rm -rf sofware/em/*.tgz

Next Steps
----------

-  Test your installation by running at least the *Small* and *Medium*
   tests mentioned in `running tests page <Running-Tests>`__.
-  Complete some of the `Scipion
   Tutorials <User-Documentation#Tutorials>`__.
