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

* Setup PYTHONPATH (assuming you're in the pythong-workshop directory)

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

* Run hisnhers.py against the example file in the data subir

   ```
    $ cd data
    $ ../bin/hisnhers.py example.fq
    Length 0: 40
    Base counts- A: 8   T: 9    C: 10   G: 13   
    Length 1: 40
    Base counts- A: 11  T: 5    C: 16   G: 8    
    Length 2: 40
    Base counts- A: 8   T: 6    C: 9    G: 17   
    Length 3: 40
    Base counts- A: 9   T: 8    C: 12   G: 11   
    Length 4: 40
    Base counts- A: 11  T: 9    C: 10   G: 10   
    Length 5: 40
    Base counts- A: 9   T: 6    C: 15   G: 10   
    Length 6: 40
    Base counts- A: 9   T: 2    C: 12   G: 17   
    Writing to example.fa
    Elapsed assembly time 5 seconds
    Elapsed annotation time 7 seconds
    Elapsed annotation time 14 seconds
   ```

* Checkout the results in example.fa.annotations

    ```
    $ head -20 example.fa.annotations
    {
        "contig2": [
            {
                "start": 2, 
                "seqid": "contig2", 
                "end": 5, 
                "key": "palindrome"
            }, 
            {
                "start": 3, 
                "seqid": "contig2", 
                "end": 5, 
                "key": "start_codon"
            }, 
            {
                "start": 31, 
                "seqid": "contig2", 
                "end": 34, 
                "key": "palindrome"
            }, 
    ```
