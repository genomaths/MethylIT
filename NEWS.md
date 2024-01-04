# MethylIT 0.3.2.7
* Current version 
* Actively used in methylome analysis  
* In development         
   * 06-23-24
      - New functions added to identify DMRs and others accumulative
        improvements.
      
# MethylIT 0.3.2.6
* Current version 
* Actively used in methylome analysis  
* In development         
   * 01-15-23:
      - Update function 'uniqueGRfilterByCov' adding parameter 'and.min.cov'.
        To get the results with previous default parameter set:
        and.min.cov = FALSE    

   * 02-07-23:      
      - Add Random Forest algorithm to functions 'evaluateDIMPclass' and
        'estimateCutPoint'.

# MethylIT 0.3.2.5
* Current version 
* Actively used in methylome analysis
* In development     
   * 05-12-22:
      - Add a new functions getDMGs wrapping several steps of DMG estimations.
      - Replace function 'data.table' with functions from dplyr package.

# MethylIT 0.3.2.3
* Current version 
* Actively used in methylome analysis
* In development     
   * 04-06-22:
      - Fixed issues with parallel computation.
   * 11-04-21:    
         1) Upgrading corresponding to R version 4.1.1.     
              - Improvement of parallel computation
              - Improvement on some non-linear algorithms.

# MethylIT 0.3.2.2
* Actively used in methylome analysis
* In development     
   * 07-14-21:    
         1) Upgrading corresponding to R version 4.1.0.     
              - Prevent potential crash originated by changes in some 
                R packages internally used by Methyl-IT.   
              - *J-Divergence* is now available in *estimateCutPoint* function

# MethylIT 0.3.2.1
* Actively used in methylome analysis
* In development     
   * 06-09-20:    
         1) Upgrading corresponding to R version 4.0.0.     
              - Prevent potential crash originated by changes in some 
                R packages internally used by Methyl-IT.

# MethylIT 0.3.2
   * 04-17-20:    
         1) A new vignette added:      
              - DMR detection with Methyl-IT
   * 03-08-20:  
         1) Two new vignette added:    
              - Principal Components and Linear Discriminant Analyses with Methyl-IT   
              - Association Between Gene Expression and Cytosine DNA Methylation at gene-body  
   * 03-23-20:  
         1) Update documentation and code style reformatted

## New features

* Improvement of the '_countTest_' function applied for DMG identification. The
  new version was named '_countTest2_'
  
* Added function 'gofReport' search for the best fitted model between the set
  of models requested by the user.
  
* The functions using machine-learning algorithms were improved

* The documentation was notably improved

# MethylIT 0.3.1

* Available at https://git.psu.edu/genomath/MethylIT

# MethylIT release

November 17, 2018:
* Initial development.    
