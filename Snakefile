with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()



rule all:
         input:
            expand("{sample}_preprocessed.rds", sample = samples),
            expand("{sample}_filtered.rds", sample = samples), 
            expand("{all}.rds", all =config['ALL']),
            expand("{all}_analysed.rds", all =config['ALL']),
            expand("{all}_annotated.rds", all =config['ALL']),

rule preprocess: 
       input:  
            "{sample}_filtered_feature_bc_matrix.h5", 
       params: 
            "{sample}" 
       output: 
          "{sample}_preprocessed.rds",
       shell: 
           """
           Rscript preprocess.R {params} 
           """ 

rule filter: 
     input: 
         "{sample}_preprocessed.rds",
     output: 
         "{sample}_filtered.rds"
     params: 
         nCount1 = config['nCount1'],
         nCount2 = config['nCount2'], 
         nfeatures1 = config['nfeatures1'],
         nfeatures2 = config['nfeatures2'], 
         mt = config['mt'] 
     shell: 
        """ 
        Rscript filter.R {input} {params[0]} {params[1]} {params[2]} {params[3]} {params[4]} 
        """ 


rule merge: 
     input: 
        expand("{sample}_filtered.rds", sample = samples)
     params: 
        all =config['ALL']   
     output: 
       expand("{all}.rds", all =config['ALL'])
     shell: 
        """
         Rscript merge.R {params.all} {input} 
        """ 

rule analyse: 
    input: 
        expand("{all}.rds", all =config['ALL'])
    output: 
        expand("{all}_analysed.rds", all =config['ALL'])
    shell: 
       """ 
       Rscript analyse.R {input}
       """

rule annotate: 
    input: 
      expand("{all}_analysed.rds", all =config['ALL']),
      annotations = config['ANNOTATIONFILE'] 
    output: 
       expand("{all}_annotated.rds", all =config['ALL'])
    shell: 
       """
       Rscript annotate.R {input[0]}  {input[1]}   
       """ 


