#before this make sure to'pip install PyWGCNA'
#note, I had to install rust and visual C++ prerequisites
#before installing PyWGCNA (not enough space tho)

#new code using PyWGCNA as a package
import os
import PyWGCNA
import pandas as pd

#chatgpt specific
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt


class coexpressionanalysis:

  def __init__(self):
    pass

  def WGCNA(self, geneExppaths, metadatapaths, variables, correlationtype, outputpath):
    #note the geneExpfilepaths, metadatapaths and variables are lists
    #of all the paths for the different files. One file per variable

    #First prepare the input to be in the right format:
    #rows/obs are samples/sale info and cols/var are genes/gene info
    
    #create a list to store the individual dataframes
    dfs = []
    for i in geneExppaths:
      #create database from csv file
      df = pd.read_csv(i, index_col=0)
      #in this dataframe, each row is a gene and each column is a sample 
      #append them to make one big dataframe
      dfs.append(df)
    #concatenate dataframes, but first check that indexes match
    if all(df.index.equals(dfs[0].index) for df in dfs):
      combined_df = pd.concat(dfs, axis=1)
    else:
      print("Row indices do not match, cannot concatenate.")
    #transpose it to get it in the right format
    geneexp = combined_df.T
    #print('transformed df')     #troubleshooting
    #print(geneexp.head())   #troubleshooting
      
    ####This section is not necessary ***
    '''
    #making it a csv file to be used as an input for PyWGCNA
    combined_df.to_csv(path_or_buf=os.path.join(outputpath,'geneexpression.csv'))
    filepath = os.path.join(outputpath,'geneexpression.csv')
    '''
    
    #initialyze WGCNA object
    wgcna = PyWGCNA.WGCNA(geneExp=geneexp, outputPath=outputpath, save=True)

    #preprocessing and module building
    wgcna.runWGCNA()
    
    #inserting metadata *****************Test this section*********
    #the variable order corresponds to the geneExppaths and metadatapaths order
    #giving each condition a color and separating by variables for later analysis
    #initiating empty list to store metadataframes
    metadfs = []
    #initializing color mapping dictionary
    color_mapping = {}
    for i in enumerate(metadatapaths):
      #create database from csv file
      metadf = pd.read_csv(metadatapaths[i], index_col=0)
      metadfs.append(metadf)
      unique_conditions = metadf['condtion'].unique()
      for j in enumerate(unique_conditions):
        #creating dictionary of colors 
        color_mapping.update({unique_conditions[j]:wgcna.metadataColors()[j]})
      wgcna.setMetadataColor(variables[i], color_mapping)
          
    if all(metadf.index.equals(metadfs[0].index) for metadf in metadfs):
      combined_metadf = pd.concat(metadfs, axis=1)
    else:
      print("Metadata row indices do not match, cannot concatenate.")

    wgcna.updateSampleInfo(sampleInfo=combined_metadf, sep=',')
    print('sample info')              #troubleshooting
    print(wgcna.geneExpr.obs.head(5)) #troubleshooting

    '''
    #### we don't need this for our data but we may want to include it
    #updating gene info
    geneList = PyWGCNA.getGeneList(dataset='mmusculus_gene_ensembl',
                                  attributes=['ensembl_gene_id', 
                                              'external_gene_name', 
                                              'gene_biotype'],
                                  maps=['gene_id', 'gene_name', 'gene_biotype'])
    wgcna.updateGeneInfo(geneList)
    '''
    #analyzing results
    wgcna.analyseWGCNA()

    #Get candidate genes for modules with high correlations to the variables of interest
    ##### if we are not assuming that all the data they give us is for variables of interest we need to add
    ##### an input that is a list of the variables of interest only

    #get the module trait correlations
    module_trait_correlation = wgcna.moduleTraitCor
    #this part finds the name of the module with the highest absolute correlation coefficients for each trait
    ###### maybe we would want something that is not absolute 
    ####### if we know low light is better than high light 
    ###### in this case we should name low light -1, control 0 and high light 1
    ###### then look for modules that are highly correlated in the way we are interested
    
    #initialize a correlated genes dataframe
    correlatedgenesdf = pd.DataFrame(columns =['Variable', 'Genes'])

    #make variables list into a dataframe for manipulation
    variablesdf = pd.DataFrame(variables)
    #look for modules that are correlated to 
    for i in enumerate(variables):
      if correlationtype[i]=='positive':
        #this will extract the module with the highest positive correlation only
        moduleofI = module_trait_correlation[variablesdf[i].nlargest(1)]
      elif correlationtype[i]=='negative':
        #this will extract the module with the highest negative correlation only
        moduleofI = module_trait_correlation[variablesdf[i].nsmallest(1)]
      correlatedgenes = wgcna.getGeneModule(moduleofI)
      #create dataframe with the genes from the modules
      correlatedgenesdf = correlatedgenesdf.append(
        {'Variable': variables[i],
         'Genes':', '.join(correlatedgenes)}, 
        ignore_index=True
      )
    #now save data frame to csv file
    correlatedgenesdf.to_csv('candidate genes from highest variable correlation')

    #for correlation visualization
    wgcna.module_trait_relationships_heatmap
    ##### maybe add network visualization, not needed tho

    #save object as a .p file
    wgcna.saveWGCNA()

    
    #important for the pathway co-expression
    #asking first for other geneIDs 
    ##***make sure this is in the right format**
    pathwaygenes = input('Please input a list of geneIDs of genes you know are involved in the pathway of interest: ')
    # Find the intersection of the requested genes and the genes in the dataset
    valid_genes = set(pathwaygenes).intersection(set(wgcna.datExpr.var_names))
    #make a list of the modules that other pathway genes belong to
    pathwaymodules = []
    for i in valid_genes:
      pathwaymodules.append(wgcna.getGeneModule(i))
    pathwaymodulesdf = pd.DataFrame(pathwaymodules)
    pathwaymodule_counts = pathwaymodulesdf.value_counts()
    print(pathwaymodule_counts) #troubleshooting
    common_module = pathwaymodule_counts.idxmax()
    potentialgenes = wgcna.getModulesGene(common_module)
    #making list into dataframe to easily make into a csv file
    potentialgenes = pd.DataFrame(potentialgenes)
    #this is one of the final products
    potentialgenes.to_csv('candidate genes pathway coexpression')
    
    pass

#In the main file, execution would be something like this

directory = r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data'
geneexpressionpaths = [r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\vst_counts.csv']
metadatapaths = [r'C:\Users\palmy\Documents\McGill\U4\BIEN 470\Workspace Justin\capstoneenv\test data\metadata.csv']
variables = []
correlationtype = []

cea = coexpressionanalysis()
cea.WGCNA(geneexpressionpaths, metadatapaths,variables, correlationtype, directory)





