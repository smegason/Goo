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

Goo will be supported to run in background mode on computing clusters in the future. Building Blender from source on computing clusters is the limiting factor at the moment. 

conda config --add channels conda-forge

libxi-devel-cos6-x86_64 gcc_linux-64 gxx_linux-64 subversion make cmake mesa-libgl-devel-cos6-x86_64 mesa-libegl-devel-cos6-x86_64 libx11-devel-cos6-x86_64 libxxf86vm-devel-cos7-x86_64 libxi-devel-cos6-x86_64 xorg-libxcursor xorg-libxrandr xorg-libxinerama libstdcxx-ng

conda install 

gcc
gcc-c++
git-lfs
subversion 
make 
cmake 
mesa-libGL-devel 
mesa-libEGL-devel 
libX11-devel 
libXxf86vm-devel 
libXi-devel 
libXcursor-devel 
libXrandr-devel 
libXinerama-devel 
libstdc++-static

qt-wayland wayland-protocols libxkbcommon dbus kernel-headers_linux-64

wayland-devel 
wayland-protocols-devel 
libxkbcommon-devel 
dbus-devel 
kernel-headers


Additional package asked when running `make`
- vulkan-tools

Prompted to install the following packages:

  libexpat           conda-forge/linux-64::libexpat-2.6.2-h59595ed_0 
  libvulkan-loader   conda-forge/linux-64::libvulkan-loader-1.3.250.0-h1fe2b44_0 
  libxcb             conda-forge/linux-64::libxcb-1.15-h0b41bf4_0 
  pthread-stubs      conda-forge/linux-64::pthread-stubs-0.4-h36c2ea0_1001 
  vulkan-tools       conda-forge/linux-64::vulkan-tools-1.3.250-hac7e632_1 
  wayland            conda-forge/linux-64::wayland-1.23.0-h5291e77_0 
  xorg-kbproto       conda-forge/linux-64::xorg-kbproto-1.0.7-h7f98852_1002 
  xorg-libx11        conda-forge/linux-64::xorg-libx11-1.8.9-h8ee46fc_0 
  xorg-libxau        conda-forge/linux-64::xorg-libxau-1.0.11-hd590300_0 
  xorg-libxdmcp      conda-forge/linux-64::xorg-libxdmcp-1.1.3-h7f98852_0 
  xorg-xextproto     conda-forge/linux-64::xorg-xextproto-7.3.0-h0b41bf4_1003 
  xorg-xproto        conda-forge/linux-64::xorg-xproto-7.0.31-h7f98852_1007 


- shaderc