# LAWS-HiC

LAWS-HiC is a computational tool designed to enhance the power of standard Hi-C peak callers by post-processing their initial results. Using a *Locally Adaptive Weighting and Screening* (LAWS) approach [1], LAWS-HiC adjusts p-values from standard Hi-C peak callers by accounting for the local spatial dependency structure of data, alleviating the common assumption of statistical independency among interactions made by most standard peak callers. This method enables more accurate identification of significant chromatin interactions.

## Installation
LAWS-HiC now is available as an R package. This package can be installed from GitHub:
```r
install.packages("devtools")
devtools::install_github("lb-zhou/LAWSHiC")
```

## Usage

This package contains three functions:
* bh.func(): this function is used to tune $\tau$. This function is based on the LAWS described by Cai et al. (2021) [1];

* laws_hic(): adjusts 3D peak calling results;

* laws_hic_paral(): adjusts 3D peak calling results with parallelized job processing.

### Input data requirement and format

The required input of LAWS-HiC are:

* results from a standard 3D peak caller (e.g., FitHiC2): only bin pair coordinators and p-values are needed for LAWS-HiC;

* domain list, in this example domain list is the list of topologically associated domains (TADs).

Here, we provide an example of 3D peak caller results and an example of domain list.

#### 3D peak caller results:

```{r}

      chr fragmentMid1 fragmentMid2    p_value
    <int>        <int>        <int>      <num>
 1:    21     16505000     15605000 1.00000000
 2:    21     16505000     15615000 0.01300292
 3:    21     16505000     15625000 1.00000000
 4:    21     16505000     15635000 0.01406848
 5:    21     16505000     15645000 0.39551740
 6:    21     16505000     15655000 1.00000000
 7:    21     16505000     15665000 0.40508290
 8:    21     16505000     15675000 1.00000000
 9:    21     16505000     15685000 0.41507940
10:    21     16505000     15695000 0.01804282
11:    21     16505000     15705000 0.01883900
12:    21     16505000     15715000 0.43088820
13:    21     16505000     15725000 0.43637470
14:    21     16505000     15735000 1.00000000
15:    21     16505000     15745000 0.44767060

```

Please check `example_data.tsv` for the full example data. The four columns of the standard 3D peak caller results should be:

* chr: chromosome number;

* fragmentMid1, fragmentMid2: coordinators of bin pairs (middle position of the fragment in this example);

* p_value: p-value from the standard peak caller.

#### Domain list:

```{r}

      chr       x1       x2
    <int>    <int>    <int>
 1:    21 18755000 18965000
 2:    21 16235000 16440000
 3:    21 15405000 15490000
 4:    21 16445000 16880000
 5:    21 15805000 16145000
 6:    21 17095000 17275000
 7:    21 16360000 16440000
 8:    21 15630000 15755000
 9:    21 17100000 17200000
10:    21 15895000 16060000

```

Please check `example_domain_list.csv` for the full example domain list data. The three columns of the domain list should be:

* chr: chromosome number;

* x1: start position of the domain;

* x2: end position of the domain.

### LAWS-HiC

For `laws_hic()` and `laws_hic_paral()` function, the input parameters are:

* input: data of standard 3D peak caller results with format previously shown;

* domain_input: data of domain list with format previously shown;

* chr: chromosome number;

* resolution: resolution of the original Hi-C data in kb (e.g., 10 for 10kb resolution);

* gc_frequency: only for `hic_laws_paral()`. Garbage collection frequency, default is 20.

`laws_hic()` is the function to adjust the 3D peak calling results with LAWS procedure:

```{r}
library(LAWSHiC)
LAWSHiC::laws_hic(input=example_data, domain_input=example_domain_list, 
                  chr=21, resolution=10)
```

For long chromosomes or large datasets, parallel processing is recommended. The `laws_hic_paral()` function enables parallel execution to enhance computational effciency.

```{r}
library(LAWSHiC)
library(doParallel)

num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)
registerDoParallel(cl)

LAWSHiC::laws_hic(input=example_data, domain_input=example_domain_list,  
                  chr = 21, resolution = 10)
```

Both functions will output the tuned $\tau$ and $h$ for each domain:

```{r}

Processed TAD boundary  1 / 12  with tau =  0.2688476  and h =  1.237808 
Processed TAD boundary  2 / 12  with tau =  0.5472209  and h =  1.210957 
Processed TAD boundary  3 / 12  with tau =  0.4205185  and h =  0.8722585 
Processed TAD boundary  4 / 12  with tau =  0.2951233  and h =  1.6203 
Processed TAD boundary  5 / 12  with tau =  0.8130685  and h =  1.479749 
Processed TAD boundary  6 / 12  with tau =  0.3582338  and h =  1.161979 
Processed TAD boundary  7 / 12  with tau =  0.8927248  and h =  0.864941 
Processed TAD boundary  8 / 12  with tau =  0.1364388  and h =  1.004806 
Processed TAD boundary  9 / 12  with tau =  0.6925968  and h =  0.9334673 
Processed TAD boundary  10 / 12  with tau =  0.595248  and h =  1.117409 
Processed TAD boundary  11 / 12  with tau =  0.8927248  and h =  0.8722585 
Processed TAD boundary  12 / 12  with tau =  0.8088843  and h =  0.864941

```

### Output format

```{r}

     chr fragmentMid1 fragmentMid2      p_value     i     j       p_laws        pi
    <int>        <int>        <int>        <num> <num> <num>        <num>     <num>
 1:    21     16135000     15835000 6.425272e-03  1614  1584 4.148520e-04 0.9393502
 2:    21     16135000     15845000 4.992228e-01  1614  1585 8.467593e-03 0.9833213
 3:    21     16135000     15855000 4.418422e-04  1614  1586 4.939890e-06 0.9889434
 4:    21     16135000     15865000 1.083029e-01  1614  1587 1.331206e-03 0.9878577
 5:    21     16135000     15875000 1.425163e-04  1614  1588 1.417433e-06 0.9901522
 6:    21     16135000     15885000 3.812398e-03  1614  1589 2.428471e-05 0.9936704
 7:    21     16135000     15895000 4.659940e-03  1614  1590 2.968562e-05 0.9936699
 8:    21     16135000     15905000 1.066683e-10  1614  1591 1.071007e-12 0.9900593
 9:    21     16135000     15915000 3.680861e-01  1614  1592 4.661362e-03 0.9874946
10:    21     16135000     15925000 2.239241e-03  1614  1593 2.409945e-05 0.9893523
11:    21     16135000     15935000 6.948450e-04  1614  1594 4.387884e-06 0.9937247
12:    21     16135000     15945000 1.187583e-02  1614  1595 3.253732e-05 0.9972677
13:    21     16135000     15955000 2.432552e-01  1614  1596 2.482740e-04 0.9989804
14:    21     16135000     15965000 4.833000e-01  1614  1597 7.484818e-04 0.9984537
15:    21     16135000     15975000 5.125350e-01  1614  1598 5.584812e-03 0.9892210

```

The columns in the output data are: 

* chr: chromosome number;

* fragmentMid1, fragmentMid2: coordinators of bin pairs (same as input);

* p_value: p-value from the standard peak caller;

* i,j: identifiers of bin pairs;

* p_laws: LAWS-HiC adjusted p-value;

* pi: LAWS-HiC estimated spatial level;

## References

[1] Cai, T. Tony, Wenguang Sun, and Yin Xia. “LAWS: A Locally Adaptive Weighting and Screening Approach to Spatial Multiple Testing.” Journal of the American Statistical Association 117, no. 539 (2021): 1370–83. doi:10.1080/01621459.2020.1859379.
