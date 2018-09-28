In order to get this working in Fedora, we need to run two commands.

First we need to run:
    yum install /usr/include/X11/xlib.h
    yum install libXtst-devel
    yum install libXrandr-devel
    yum install mesa-libGL-devel
    yum groupinstall "X Software Development"
