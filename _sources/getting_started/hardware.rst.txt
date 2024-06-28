.. _hardware:

Hardware recommendations
================================================

Goo has been developed with two main machines in hands: 

1. M2 Max MacBook Pro 32GB Unified Memory
2. RTX 3080 GPU, Intel i9 CPU and 64GB RAM

Goo exhibits a sub-quadratic time complexity on the number of cells, typically around :math:`O(N^{1.5})` where :math:`N` is the number of cells in the simulations. This is achieved because Blender restrict the computation of interations to a maximum distance. In other words, a force that is 3 cell diameter away from a certain cell will not affect it, and the calculation of the interaction is ignored by Blender engines. 


GPU
-----

It is very complex to accellerate physics-based simulations with GPUs: they are inherently serial in their computation. 

Rendering simulations really benefit from GPUs. 

Running Goo on computing clusters
-----------------------------------

Goo will be supported to run in background mode on computing clusters in the future. Building Blender from source on computing clusters is the limiting factor at the moment. We are exploring docker containers. 

O2 is Harvard Medical School's computing cluster. If you don't have an O2 account, follow RC's `instructions <https://harvardmed.atlassian.net/wiki/spaces/O2/pages/1918304257/How+to+request+or+retain+an+O2+account>`__.

1. Log in to O2

.. code-block:: bash

    ssh <HMS-ID>@o2.hms.harvard.edu


Work in progress. Stay tuned for updates.