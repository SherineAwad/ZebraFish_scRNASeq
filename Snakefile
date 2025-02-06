with open(config['SAMPLES']) as fp:
    samples = fp.read().splitlines()



rule all:
         input:
            expand("{sample}_pre.rds", sample = samples),
            expand("{sample}_filter.rds", sample = samples), 
            expand("{all}.rds", all =config['ALL']),
            expand("{all}_analysed.rds", all =config['ALL']),
            expand("{all}_annotated.rds", all =config['ALL']),
            #expand("{selection}_analysed.rds", selection = config['SNAME']), 
            #expand("{selection}_annotated.rds", selection = config['SNAME']),

rule preprocess: 
       input:  
            "{sample}_filtered_feature_bc_matrix.h5", 
       params: 
            "{sample}" 
       output: 
          "{sample}_pre.rds",
       shell: 
           """
           Rscript preprocess.R {params} 
           """ 

rule filter: 
     input: 
         "{sample}_pre.rds",
     output: 
         "{sample}_filter.rds"
     shell: 
        """ 
        Rscript filter.R {input} 
        """ 


rule merge: 
     input: 
        expand("{sample}_filter.rds", sample = samples)
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
      expand("{all}_analysed.rds", all =config['ALL']) 
    output: 
       expand("{all}_annotated.rds", all =config['ALL'])
    shell: 
       """
       Rscript annotateClusters.R {input}   
       """ 


