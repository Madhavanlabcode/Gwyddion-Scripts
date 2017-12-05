# Gwyddion-Scripts

This repository is meant to host python scripts written and maintained by the Madhavan group to be used with the software Gwyddion,.

## What is Gwyddion? 

From the website, "Gwyddion is a modular program for SPM (scanning probe microscopy) data visualization and analysis". Gwyddion is a powerful, modular open source sofware that is very convenient for much of the data analysis we do. However, Gwyddion does not have all of the functionality we need for data analysis. Fortunately, Gwyddion allows for python scripting to add on additional functionality. That's where this repository come in.

## How do I install Gwyddion to use these scripts?

Enabling python scripting or Pygwy as it called in Gwyddion is not trivial for several reasons having to do with python versions, pygwy dependencies, and different versions of Gwyddion. However here are the steps that are recommended from the Gwyddion website found ![here](http://gwyddion.net/documentation/user-guide-en/installation-ms-windows.html#installation-ms-windows-pygwy).

1. ![Download the 32-bit version of Gwyddion](http://gwyddion.net/download.php#stable-windows)
2. ![Download Python 2.7.13](https://www.python.org/downloads/release/python-2713/)
3. ![Download three Python packages Gwyddion needs](https://sourceforge.net/projects/gwyddion/files/pygtk-win32/)

Once you have downloaded and installed all of the things in the list, Gwyddion should be able to recognize and utilize your scripts. Upon starting, gwyddion recognizes python files placed in a specific folder on your computer. On Windows this folder is ~\gwyddion\pygwy and on Unix OS's its ~\.gwyddion\pygwy. This is where you should clone this repository if you want to use these scripts.

After getting everything set up, I recommend reading the Gwyddion tutorial about how to use pygwy and python scripting with their software, found ![here](http://gwyddion.net/documentation/user-guide-en/pygwy.html). It will go through basic file format and give you an idea of what is capabable with pygwy/python scripting.

If you decide to make your own scripts, you will undoubtably look for Gwyddion's ![python API](http://gwyddion.net/documentation/head/pygwy/), which will allow you to access the data from the loaded images/files.

## List of Necessary Python Libraries

In order to have the full capabiltiy of all the scripts here, you need have some extra python modules downloaded. As you write scripts and want to include libraries, think if you can use ones that are already required or are robust to updates. If you do add a necessary library, place it here.

* NumPy (1.13)
