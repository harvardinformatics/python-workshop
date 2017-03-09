# python-workshop

Harvard FAS Informatics materials for the Practical Python on Odyssey

### Prerequisites
Workshop participants must have basic Unix proficiency, some operational knowledge of Python, and an Odyssey account.

Python basics will be covered, but not in depth.


### Production branch contains the answer
A checkout of the production branch will provide a hisnhers.py script that works.

To try out the functional script follow these steps:

* Checkout the production branch of the project

    ```
    $ git clone https://github.com/harvardinformatics/python-workshop.git
    $ cd python-workshop
    $ git checkout production
    ```

* Setup PYTHONPATH

    ```
    $ export PYTHONPATH=`pwd`:$PYTHONPATH
    ```

* Install the lookkool annotation package; make sure gcc is 4.9 or better

    ```
    $ module load gcc/4.9.3-fasrc01
    $ pip install git+https://github.com/harvardinformatics/lookkool.git
    ```

* Set LD_LIBRARY_PATH to find the liblookkool.so.  It should be in site-packages.

    ```
    $ export LD_LIBRARY_PATH=/path/to/site-packages:$LD_LIBRARY_PATH
    ```

