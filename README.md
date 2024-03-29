# BOPS: <br><small>small protein complex prediction algorithm **B**ased **O**n **P**rotein-protein interaction network **S**egmentation</small>

## Table of contents
* [Features](#Features)
* [Algorithm](#Process)
* [Method](#Method)
* [Experiment](#Experiment)
* [Email](#Email)
* [Acknowledgments](#Acknowledgments)

## Features
The BOPS algorithm can find small protein complexes from the protein protein interaction network to help biological experiments.

## Algorithm
First, the BOPS algorithm calculates the balanced weights, and replaces the original weights with balanced weights. Second, the BOPS algorithm divides the graphs larger than MAXP until the original PPIN is divided into small networks and the details of segmentation is described in (3.4.2). Third, the BOPS algorithm enumerates the connected subset of each small network, calculates the cohesion of each connected subset, identifies potential protein complexes based on cohesion and removes those that are similar.

## Method

#### Compiling
The BOPS algorithm is coded in C++. You can use C++11 or or higher version of the compiler to get the program.

Windows：
```
g++ BOPS.cpp -o BOPS.exe -std=c++11
```
Linux:
```
g++ BOPS.cpp -o BOPS -std=c++11
```
#### Runing
Run the BOPS algorithm to get the result as follow:

Windows:
```
BOPS.exe PPI_file result_file
```
Linux:
```
./BOPS PPI_file result_file
```

Specifically, you can input balanced_index on your own, the default value of it is 1.5.

Windows:
```
BOPS.exe PPI_file result_file balanced_index
```
Linux:
```
./BOPS PPI_file result_file balanced_index
```

## Experiment

The "evaluation" folder includes the methods used to evaluate the method and the resulting predicted complexes based on:

- the number of predicted complexes
- Fscore: the harmonic mean of Recall and Precision
- Precision: the rate of predicted protein complexes that match at least one reference complex 
- Recall: the rate of reference protein complexes that match at least one predicted complex
- ACC: the geometric accuracy
- Sn: the rate of the maximum-sum number of matched proteins to the total number of proteins in the set of reference protein complex.
- PPV: the rate of the maximumsum number of matched proteins to the total matched number of proteins in the set of predicted protein complex. 

To evaluate the "Predicted Complexes.txt" file, it runs as follows:
The match is coded in python. You must use python 2.7 of the interpreter to get the program.
			
  			python match.py -n ppi_name.txt reference_name.txt result_name.txt
            

The parameters of these methods are set as the recommended values as mentioned in their original papers. For our method, we set the balance index β to 1.5 as to recommended values.

| Datasets         | \#predicated | F\-score | ACC    |
|------------------|--------------|----------|--------|
| Krogan\-core     | 704          | 0\.558   | 0\.528 |
| Krogan\-extended | 778          | 0\.538   | 0\.482 |
| Gavin            | 832          | 0\.668   | 0\.560 |
| Collins          | 794          | 0\.614   | 0\.550 |

The experimental performance proves that the BOPS algorithm can obtain the best results on F-score and ACC when identifying small protein complexes. And the performance of BOPS is better than most of algorithms with respect to the whole protein complexes.

## Email
If you have any question, you can contact with jiaqinglv@foxmail.com,22151303@zju.edu.cn,liangbing@dlut.edu.cn and zhyj@dmu.edu.cn.
## Acknowledgments
The project is supported in part by College Student Innovation and Entrepreneurship Training Program Support Project of China (2019101411600010093,2020101411600010129). We are very grateful for the research done by our predecessors.

At last, thanks very much for using BOPS.
